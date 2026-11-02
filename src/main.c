#define AETHER_IMPLEMENTATION
#include "aether/aether.h"

#include "laplace/laplace.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*-------- F L O W   F I E L D ----------------------------------------------*/

/*
 * Read-only view of a laplace solution file. Field pointers point into the
 * memory-mapped file. Coordinates are in file units (length_scale metres per
 * unit); all other fields are SI.
 */
typedef struct FlowField
{
    i64 ni, nj, nk;
    f64 x0, x1;          /* streamwise domain bounds, file units */
    f64 length_scale;    /* metres per coordinate unit */
    f64 gamma, R_gas;    /* -, J/(kg K) */

    const f64 *x;
    const f64 *y;
    const f64 *z;
    const f64 *u;
    const f64 *v;
    const f64 *w;
    const f64 *p;
    const f64 *T;
    const f64 *rho;
} FlowField;

/* Air state at a point, SI */
typedef struct FlowSample { f64 u, v, w, rho, p, T; } FlowSample;

typedef u8 Location;
enum Location_
{
    Location_Invalid = 0,       /* NaN position, or beyond a face a droplet cannot cross */
    Location_Inside,
    Location_WallImpact,
    Location_Exited,            /* beyond an inlet or outlet face */
};

/* Field array stored under `name`, or NULL if the file has no such field. Arrays start at hdr->header_size. */
static const f64* flow_field_find(const LaplaceHeader* hdr, view raw, i64 n, const char* name)
{
    const f64* base = (const f64*)(raw.data + hdr->header_size);
    for (i64 q = 0; q < hdr->n_fields; ++q)
        if (strncmp(hdr->field_names[q], name, sizeof(hdr->field_names[q])) == 0)
            return base + q * n;
    return NULL;
}

/*
 * Checks the header of a mapped laplace file and points f's fields into the
 * mapping. On failure returns false and describes the first problem found in
 * *out_err (allocated from arena); *f is left untouched.
 */
static b8 flow_field_open(view raw, Arena* arena, LaplaceHeader* hdr, FlowField* f, str8* out_err)
{
#define FAIL(...) do { *out_err = str8_push_fmt(arena, __VA_ARGS__); return false; } while (0)

    if (raw.size < sizeof(*hdr)) FAIL("file is smaller than the %llu-byte header", (unsigned long long)sizeof(*hdr));
    memcpy(hdr, raw.data, sizeof(*hdr));

    if (hdr->magic != LAPLACE_MAGIC)
        FAIL("not a laplace file (bad magic number)");
    if (hdr->version != LAPLACE_VERSION)
        FAIL("file format version %llu, expected %llu", (unsigned long long)hdr->version, (unsigned long long)LAPLACE_VERSION);
    if (hdr->header_size < sizeof(*hdr) || hdr->header_size > raw.size || hdr->header_size % sizeof(f64) != 0)
        FAIL("header_size %llu must be a multiple of 8, at least %llu, and within the file",
             (unsigned long long)hdr->header_size, (unsigned long long)sizeof(*hdr));
    if (hdr->ni < 2 || hdr->nj < 2 || hdr->nk < 2)
        FAIL("grid %lld x %lld x %lld needs at least 2 nodes per direction", (long long)hdr->ni, (long long)hdr->nj, (long long)hdr->nk);
    if (hdr->n_fields < 1 || hdr->n_fields > LAPLACE_MAX_FIELDS)
        FAIL("field count %lld outside 1..%d", (long long)hdr->n_fields, LAPLACE_MAX_FIELDS);

    /* bound each dimension by the file size first, so ni*nj*nk cannot overflow */
    u64 cap = raw.size / sizeof(f64);
    u64 ni = (u64)hdr->ni, nj = (u64)hdr->nj, nk = (u64)hdr->nk;
    if (ni > cap || nj > cap / ni || nk > cap / (ni * nj))
        FAIL("grid %llu x %llu x %llu is larger than the file", (unsigned long long)ni, (unsigned long long)nj, (unsigned long long)nk);

    i64 n = hdr->ni * hdr->nj * hdr->nk;
    u64 expected = hdr->header_size + (u64)hdr->n_fields * (u64)n * sizeof(f64);
    if (raw.size != expected)
        FAIL("file is %llu bytes, header implies %llu", (unsigned long long)raw.size, (unsigned long long)expected);

    if (!(hdr->gamma > 1.0))        FAIL("gamma = %g, must be > 1", hdr->gamma);
    if (!(hdr->R_gas > 0.0))        FAIL("R_gas = %g, must be > 0", hdr->R_gas);
    if (!(hdr->length_scale > 0.0)) FAIL("length_scale = %g, must be > 0", hdr->length_scale);

    FlowField out    = {0};
    out.ni           = hdr->ni;
    out.nj           = hdr->nj;
    out.nk           = hdr->nk;
    out.length_scale = hdr->length_scale;
    out.gamma        = hdr->gamma;
    out.R_gas        = hdr->R_gas;

    struct { const char* name; const f64** dst; } required[] =
    {
        { "x",   &out.x   }, { "y", &out.y }, { "z", &out.z },
        { "u",   &out.u   }, { "v", &out.v }, { "w", &out.w },
        { "rho", &out.rho }, { "p", &out.p }, { "T", &out.T },
    };
    for (u64 q = 0; q < ARRAY_COUNT(required); ++q)
    {
        *required[q].dst = flow_field_find(hdr, raw, n, required[q].name);
        if (!*required[q].dst) FAIL("required field \"%s\" is not in the file", required[q].name);
    }

    out.x0 = out.x1 = out.x[0];
    for (i64 i = 1; i < out.ni; ++i)
    {
        f64 xi = out.x[laplace_idx(i, 0, 0, out.nj, out.nk)];
        out.x0 = MIN(out.x0, xi);
        out.x1 = MAX(out.x1, xi);
    }

    *f = out;
    return true;

#undef FAIL
}

/*-------- F L O W   L O O K U P --------------------------------------------*/

/*
 * Lookup tables for a separable grid, in metres. Every i-plane lies at one
 * x, and on station i
 *
 *     y(i, j) = y_lo[i] + (y_hi[i] - y_lo[i]) * s[j]
 *     z(i, k) = z_lo[i] + (z_hi[i] - z_lo[i]) * t[k]
 *
 * with the same distributions s and t (0 at the first node, 1 at the last)
 * on every station. Within a cell the trilinear map of the corner nodes then
 * separates, so a point's cell and local coordinates follow exactly from
 * three 1-D searches, with no iteration.
 */
typedef struct FlowGrid
{
    i64 ni, nj, nk;
    const f64 *xs;              /* station x (ni), strictly increasing */
    const f64 *y_lo, *y_hi;     /* y at j = 0 and j = nj-1 on each station (ni) */
    const f64 *z_lo, *z_hi;     /* z at k = 0 and k = nk-1 on each station (ni) */
    const f64 *s;               /* distribution in j (nj), strictly increasing from 0 to 1 */
    const f64 *t;               /* distribution in k (nk), strictly increasing from 0 to 1 */
    LaplaceBC  face_bcs[6];     /* i-min, i-max, j-min, j-max, k-min, k-max */
} FlowGrid;

static const char* const face_names[6] = { "i-min", "i-max", "j-min", "j-max", "k-min", "k-max" };

/*
 * Builds g's tables from f (allocated from arena) and checks that every node
 * lies where the tables put it, to 1e-9 of the largest coordinate; that each
 * symmetry face is a plane of constant y or z; and that no i-face is a
 * symmetry face. On failure returns false and describes the first problem in
 * *out_err; *g is left untouched.
 */
static b8 flow_grid_build(const FlowField* f, const LaplaceHeader* hdr, Arena* arena, FlowGrid* g, str8* out_err)
{
#define FAIL(...) do { *out_err = str8_push_fmt(arena, __VA_ARGS__); return false; } while (0)

    const i64 ni = f->ni, nj = f->nj, nk = f->nk;
    const f64 L  = f->length_scale;

    f64* xs   = arena_push_array(arena, f64, ni);
    f64* y_lo = arena_push_array(arena, f64, ni);
    f64* y_hi = arena_push_array(arena, f64, ni);
    f64* z_lo = arena_push_array(arena, f64, ni);
    f64* z_hi = arena_push_array(arena, f64, ni);
    f64* s    = arena_push_array(arena, f64, nj);
    f64* t    = arena_push_array(arena, f64, nk);

    /* station tables from the edges of each i-plane */
    f64 extent = 0;
    for (i64 i = 0; i < ni; ++i)
    {
        xs[i]   = f->x[laplace_idx(i, 0,      0,      nj, nk)] * L;
        y_lo[i] = f->y[laplace_idx(i, 0,      0,      nj, nk)] * L;
        y_hi[i] = f->y[laplace_idx(i, nj - 1, 0,      nj, nk)] * L;
        z_lo[i] = f->z[laplace_idx(i, 0,      0,      nj, nk)] * L;
        z_hi[i] = f->z[laplace_idx(i, 0,      nk - 1, nj, nk)] * L;

        if (!(y_hi[i] > y_lo[i]) || !(z_hi[i] > z_lo[i]))
            FAIL("station i = %lld: y or z does not increase from the first node to the last", (long long)i);
        if (i > 0 && !(xs[i] > xs[i - 1]))
            FAIL("station x is not strictly increasing at i = %lld", (long long)i);

        extent = MAX(extent, MAX(fabs(xs[i]), MAX(MAX(fabs(y_lo[i]), fabs(y_hi[i])), MAX(fabs(z_lo[i]), fabs(z_hi[i])))));
    }

    /* cross-stream distributions, from station 0 */
    for (i64 j = 0; j < nj; ++j) s[j] = (f->y[laplace_idx(0, j, 0, nj, nk)] * L - y_lo[0]) / (y_hi[0] - y_lo[0]);
    for (i64 k = 0; k < nk; ++k) t[k] = (f->z[laplace_idx(0, 0, k, nj, nk)] * L - z_lo[0]) / (z_hi[0] - z_lo[0]);
    for (i64 j = 1; j < nj; ++j) if (!(s[j] > s[j - 1])) FAIL("y is not strictly increasing in j at j = %lld", (long long)j);
    for (i64 k = 1; k < nk; ++k) if (!(t[k] > t[k - 1])) FAIL("z is not strictly increasing in k at k = %lld", (long long)k);

    /* separability: every node must sit where the tables put it */
    const f64 tol = 1e-9 * extent;
    for (i64 i = 0; i < ni; ++i)
    {
        for (i64 k = 0; k < nk; ++k)
        {
            for (i64 j = 0; j < nj; ++j)
            {
                i64 m  = laplace_idx(i, j, k, nj, nk);
                f64 dx = f->x[m] * L - xs[i];
                f64 dy = f->y[m] * L - (y_lo[i] + (y_hi[i] - y_lo[i]) * s[j]);
                f64 dz = f->z[m] * L - (z_lo[i] + (z_hi[i] - z_lo[i]) * t[k]);
                if (!(fabs(dx) <= tol && fabs(dy) <= tol && fabs(dz) <= tol))
                    FAIL("grid is not separable: node (i,j,k) = (%lld, %lld, %lld) is %.3g m from its separable position",
                         (long long)i, (long long)j, (long long)k, MAX(fabs(dx), MAX(fabs(dy), fabs(dz))));
            }
        }
    }

    /* a symmetry face must be a plane of constant y or z, so reflection across it is a coordinate flip */
    const f64* face_coord[6] = { NULL, NULL, y_lo, y_hi, z_lo, z_hi };
    for (i32 q = 0; q < 6; ++q)
    {
        if (hdr->face_bcs[q] != LaplaceBC_Symmetry) continue;
        if (!face_coord[q]) FAIL("a symmetry condition on the %s face is not supported", face_names[q]);
        for (i64 i = 1; i < ni; ++i)
            if (!(fabs(face_coord[q][i] - face_coord[q][0]) <= tol))
                FAIL("symmetry face %s is not a plane: it moves %.3g m between i = 0 and i = %lld",
                     face_names[q], fabs(face_coord[q][i] - face_coord[q][0]), (long long)i);
    }

    FlowGrid out = { .ni = ni, .nj = nj, .nk = nk, .xs = xs, .y_lo = y_lo, .y_hi = y_hi,
                     .z_lo = z_lo, .z_hi = z_hi, .s = s, .t = t };
    memcpy(out.face_bcs, hdr->face_bcs, sizeof(out.face_bcs));
    *g = out;
    return true;

#undef FAIL
}

/*
 * Cell index c in [0, n-2] with tab[c] <= v <= tab[c+1], for strictly
 * increasing tab and tab[0] <= v <= tab[n-1]. When hint is a valid cell the
 * search starts there and gallops outward (1, 2, 4, ... cells) before
 * bisecting, so a value in the hinted cell costs two comparisons and one m
 * cells away costs O(log m). Without a hint it bisects the whole table.
 */
static i64 search_cell(const f64* tab, i64 n, f64 v, i64 hint)
{
    i64 lo = 0, hi = n - 1;                         /* invariant: tab[lo] <= v <= tab[hi] */
    if (hint >= 0 && hint <= n - 2)
    {
        if (v < tab[hint])
        {
            hi = hint;
            for (i64 step = 1; ; step *= 2)
            {
                i64 probe = hint - step;
                if (probe <= 0)       { lo = 0; break; }
                if (tab[probe] <= v)  { lo = probe; break; }
                hi = probe;
            }
        }
        else if (v > tab[hint + 1])
        {
            lo = hint + 1;
            for (i64 step = 1; ; step *= 2)
            {
                i64 probe = hint + 1 + step;
                if (probe >= n - 1)   { hi = n - 1; break; }
                if (tab[probe] >= v)  { hi = probe; break; }
                lo = probe;
            }
        }
        else
        {
            return hint;
        }
    }
    while (hi - lo > 1)
    {
        i64 mid = lo + (hi - lo) / 2;
        if (tab[mid] <= v) lo = mid; else hi = mid;
    }
    return lo;
}

/* exact at both ends: w = 0 gives a, w = 1 gives b */
static f64 lerp(f64 a, f64 b, f64 w) { return a * (1.0 - w) + b * w; }

/* Cell found by a droplet's previous lookup; the next search starts there. i < 0 means none. */
typedef struct CellHint { i64 i, j, k; } CellHint;

#define CELL_HINT_NONE ((CellHint){ -1, -1, -1 })

/* Where a point lies in the grid */
typedef struct CellCoord
{
    i64 i, j, k;            /* cell, named by its lowest-index corner node */
    f64 xi, eta, zeta;      /* local coordinates in the cell, each in [0, 1] */
    f64 sy, sz;             /* -1 if the point was reflected across a symmetry plane of constant y (z), else +1 */
    i32 face;               /* when not inside: the face the point is beyond, 0..5 = i-min..k-max; else -1 */
} CellCoord;

/* what lying beyond a face with this boundary condition means for a droplet */
static Location location_beyond(LaplaceBC bc)
{
    switch (bc)
    {
        case LaplaceBC_Wall:   return Location_WallImpact;
        case LaplaceBC_Inlet:
        case LaplaceBC_Outlet: return Location_Exited;
        default:               return Location_Invalid;
    }
}

/*
 * Locates the point (x, y, z), in metres. A point beyond a symmetry plane is
 * reflected into the computed domain first. Inside the domain, fills the
 * cell and local coordinates and returns Location_Inside; otherwise sets
 * out->face to the face the point is beyond and returns what that face's
 * boundary condition means. When hint is non-NULL the search starts at the
 * hinted cell, and on success the hint is updated to the cell found.
 */
static Location flow_locate(const FlowGrid* g, f64 x, f64 y, f64 z, CellHint* hint, CellCoord* out)
{
    CellCoord c = { .sy = 1.0, .sz = 1.0, .face = -1 };
    *out = c;
    if (isnan(x) || isnan(y) || isnan(z)) return Location_Invalid;

    /* symmetry planes are constant in x (checked at build), so station 0 holds their position */
    if      (g->face_bcs[2] == LaplaceBC_Symmetry && y < g->y_lo[0]) { y = 2.0 * g->y_lo[0] - y; c.sy = -1.0; }
    else if (g->face_bcs[3] == LaplaceBC_Symmetry && y > g->y_hi[0]) { y = 2.0 * g->y_hi[0] - y; c.sy = -1.0; }
    if      (g->face_bcs[4] == LaplaceBC_Symmetry && z < g->z_lo[0]) { z = 2.0 * g->z_lo[0] - z; c.sz = -1.0; }
    else if (g->face_bcs[5] == LaplaceBC_Symmetry && z > g->z_hi[0]) { z = 2.0 * g->z_hi[0] - z; c.sz = -1.0; }

    if      (x < g->xs[0])         c.face = 0;
    else if (x > g->xs[g->ni - 1]) c.face = 1;
    if (c.face >= 0) { *out = c; return location_beyond(g->face_bcs[c.face]); }

    c.i  = search_cell(g->xs, g->ni, x, hint ? hint->i : -1);
    c.xi = (x - g->xs[c.i]) / (g->xs[c.i + 1] - g->xs[c.i]);

    /* cross-section bounds at x, then the position as a fraction of them: 0 on the first node of j (k), 1 on the last */
    f64 y0 = lerp(g->y_lo[c.i], g->y_lo[c.i + 1], c.xi), y1 = lerp(g->y_hi[c.i], g->y_hi[c.i + 1], c.xi);
    f64 z0 = lerp(g->z_lo[c.i], g->z_lo[c.i + 1], c.xi), z1 = lerp(g->z_hi[c.i], g->z_hi[c.i + 1], c.xi);
    f64 sv = (y - y0) / (y1 - y0);
    f64 tv = (z - z0) / (z1 - z0);

    if      (sv < 0.0) c.face = 2;
    else if (sv > 1.0) c.face = 3;
    else if (tv < 0.0) c.face = 4;
    else if (tv > 1.0) c.face = 5;
    if (c.face >= 0) { *out = c; return location_beyond(g->face_bcs[c.face]); }

    c.j    = search_cell(g->s, g->nj, sv, hint ? hint->j : -1);
    c.k    = search_cell(g->t, g->nk, tv, hint ? hint->k : -1);
    c.eta  = (sv - g->s[c.j]) / (g->s[c.j + 1] - g->s[c.j]);
    c.zeta = (tv - g->t[c.k]) / (g->t[c.k + 1] - g->t[c.k]);

    if (hint) *hint = (CellHint){ c.i, c.j, c.k };
    *out = c;
    return Location_Inside;
}

/* The 8 corner nodes of a cell and their trilinear weights. Corner q is offset by (q & 1, q >> 1 & 1, q >> 2 & 1) in (i, j, k). */
typedef struct CellStencil { i64 m[8]; f64 w[8]; } CellStencil;

static CellStencil cell_stencil(const FlowField* f, const CellCoord* c)
{
    CellStencil st;
    for (i32 q = 0; q < 8; ++q)
    {
        i32 di = q & 1, dj = (q >> 1) & 1, dk = (q >> 2) & 1;
        st.m[q] = laplace_idx(c->i + di, c->j + dj, c->k + dk, f->nj, f->nk);
        st.w[q] = (di ? c->xi : 1.0 - c->xi) * (dj ? c->eta : 1.0 - c->eta) * (dk ? c->zeta : 1.0 - c->zeta);
    }
    return st;
}

static f64 stencil_apply(const CellStencil* st, const f64* field)
{
    f64 sum = 0;
    for (i32 q = 0; q < 8; ++q) sum += st->w[q] * field[st->m[q]];
    return sum;
}

/* Air state at a located point. Velocity normal to a symmetry plane the point was reflected across is flipped back. */
static FlowSample flow_sample(const FlowField* f, const CellCoord* c)
{
    CellStencil st = cell_stencil(f, c);
    FlowSample  s;
    s.u   = stencil_apply(&st, f->u);
    s.v   = stencil_apply(&st, f->v) * c->sy;
    s.w   = stencil_apply(&st, f->w) * c->sz;
    s.rho = stencil_apply(&st, f->rho);
    s.p   = stencil_apply(&st, f->p);
    s.T   = stencil_apply(&st, f->T);
    return s;
}

/* flow_locate, then flow_sample when the point is inside; *out is untouched otherwise */
Location locate_and_sample(const FlowField* f, const FlowGrid* g, f64 x, f64 y, f64 z, CellHint* hint, FlowSample* out)
{
    CellCoord c;
    Location  loc = flow_locate(g, x, y, z, hint, &c);
    if (loc == Location_Inside) *out = flow_sample(f, &c);
    return loc;
}

/*-------- D R O P L E T S --------------------------------------------------*/

typedef u8 DropletStatus;
enum DropletStatus_
{
    DropletStatus_Active = 0,
    DropletStatus_Impacted,
    DropletStatus_Evaporated,
    DropletStatus_Exited,
};

static const char* droplet_status_name(DropletStatus s)
{
    switch (s)
    {
        case DropletStatus_Impacted:   return "impacted";
        case DropletStatus_Evaporated: return "evaporated";
        case DropletStatus_Exited:     return "exited";
        default: return "max_steps";
    }
}

typedef struct Droplet
{
    f64 x, y, z;
    f64 u, v, w;
    f64 T_mean;
    f64 d, d0;
    i64 last_cell_idx;
    DropletStatus status;
} Droplet;

static str8 format_bytes(Arena* arena, u64 n)
{
    if (n >= 1024ull * 1024ull * 1024ull) { return str8_push_fmt(arena, "%.2f GB", (f64)n / (1024.0 * 1024 * 1024)); }
    if (n >= 1024ull * 1024ull )          { return str8_push_fmt(arena, "%.2f MB", (f64)n / (1024.0 * 1024)); }
    if (n >= 1024ull )                    { return str8_push_fmt(arena, "%.2f kB", (f64)n / (1024.0)); }
    return str8_push_fmt(arena, "%llu bytes", (unsigned long long)n);
}


int main(int argc, char** argv)
{

    int exit_code = EXIT_FAILURE;
    u64 program_start = time_mark();

    if (args_has_flag(argc, argv, STR("--help"), STR("-h")))
        printf("lagrange.exe\n  -- help text comming soon\n");

    Arena* memory  = arena_alloc(GB(1));
    Arena* scratch = arena_alloc(MB(4));

    const char* path = "data/flow.bin";
    view raw = file_map(path);
    if (!raw.data)
    {
        fprintf(stderr, "lagrange: cannot open flow file '%s'\n", path);
        goto cleanup_arena;
    }

    LaplaceHeader hdr;
    FlowField     f;
    str8          err;
    if (!flow_field_open(raw, scratch, &hdr, &f, &err))
    {
        fprintf(stderr, "lagrange: %s: %.*s\n", path, (int)err.size, (const char*)err.data);
        goto cleanup_files;
    }

    FlowGrid grid;
    if (!flow_grid_build(&f, &hdr, memory, &grid, &err))
    {
        fprintf(stderr, "lagrange: %s: %.*s\n", path, (int)err.size, (const char*)err.data);
        goto cleanup_files;
    }

    str8 size = format_bytes(scratch, raw.size);
    fprintf(stdout, "flow file: %s  (%lld x %lld x %lld nodes, %.*s)\n", path,
            (long long)f.ni, (long long)f.nj, (long long)f.nk, (int)size.size, (const char*)size.data);


    {

    }


    /* output results */
    {

    }


    f64 elapsed = time_elapsed_sec(program_start, time_mark());
    fprintf(stdout, "Program executed successfully --- %.4f s\n", elapsed);

    exit_code = EXIT_SUCCESS;

cleanup_files:
    file_unmap(raw);
cleanup_arena:
    arena_release(scratch);
    arena_release(memory);
exit:
    return exit_code;
}
