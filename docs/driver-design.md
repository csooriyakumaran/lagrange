# Trajectory Driver: Design and Implementation Plan

`src/main.c` is a droplet-trajectory driver for an icing wind tunnel. It reads
a converged potential-flow solution from `laplace`, releases droplets at the
spray-bar plane, tracks each one through the contraction to the end of the
domain, and reports the droplet state at chosen sampling planes and where
each trajectory ends. The droplet physics (`droplet_step`) is described in
`docs/model-cleanup-and-port.md`; this document covers everything around it,
plus the full implementation order.

## Status

- **Done — flow-file loading.** `flow_field_open` memory-maps the file,
  validates the header (magic, version, `header_size`, grid dimensions
  against the file size, gas properties), and finds the required fields by
  name. On failure it reports the first problem and the program exits. The
  program currently prints a one-line summary of the loaded file.
- **Done — flow lookup.** `flow_grid_build` builds the lookup tables and
  checks the grid at load; `flow_locate`, `flow_sample` and
  `locate_and_sample` find a point's cell and interpolate the air state
  ("Flow lookup" below).
- **Separate tool — flow-file inspection.** Mesh and domain statistics, flow
  at stations, stagnation conditions and conservation checks are printed by
  `laplace-inspect` in the `laplace` repository
  (`application/c/inspect.c`). Run it on a flow file before using the file
  here.
- **Not written:** `droplet_step`, the sweep loop, time-step control, plane
  records, the run log.

## Pipeline

```
data/flow.bin --> FlowField (mmap'd, read-only) + FlowGrid (tables built and checked at load)
                     |
                     v
     release grid at the spray-bar plane: diameter x (y, z) offset --> Droplet[]
                     |
                     v
     per droplet, to completion:
         repeat { locate_and_sample; droplet_step; record sampling-plane crossings }
                     |
                     v
     PlaneRecord[n_droplets x n_planes] + outcome per droplet
                     |
                     v
     fixed-width run log
```

## Flow field

`data/flow.bin` is the node-centred output of the `laplace` potential-flow
solver on a body-fitted structured grid with flat index
`j + nj*(k + nk*i)`. Axes: x streamwise (index i), y lateral (index j), z
vertical, positive up (index k). Coordinates are stored in file units
(`length_scale` metres per unit, 1e-3 for mm); all other fields are SI. The
header records the boundary condition on each of the six faces.

```c
typedef struct FlowField {
    i64 ni, nj, nk;
    f64 x0, x1;              /* streamwise domain bounds, file units */
    f64 length_scale;        /* metres per coordinate unit */
    f64 gamma, R_gas;        /* -, J/(kg K) */
    const f64 *x, *y, *z;
    const f64 *u, *v, *w;
    const f64 *p, *T, *rho;
} FlowField;
```

## Flow lookup

### Grid structure

Grids from `laplace_grid` are separable:

- every i-plane lies at a single x;
- on station i, `y(i, j) = y_lo(i) + (y_hi(i) - y_lo(i)) * s_j` and
  `z(i, k) = z_lo(i) + (z_hi(i) - z_lo(i)) * t_k`;
- the distributions s and t (0 at the first node, 1 at the last) are the
  same on every station.

Storing both bounds per station, rather than assuming a symmetry plane at 0,
also covers full-domain grids and sections whose opposite walls follow
different contours.

`flow_grid_build` builds these tables in metres from the mapped coordinates
and refuses the file, naming the first failure, if:

- station x is not strictly increasing, or y or z does not increase from the
  first node to the last;
- any node lies further than 1e-9 of the largest coordinate from where the
  tables put it;
- a symmetry face is not a plane of constant y or z, or an i-face is a
  symmetry face.

The check reads every coordinate once, about 15–25 ms for the current
2.1 million nodes.

```c
typedef struct FlowGrid {
    i64 ni, nj, nk;
    const f64 *xs;              /* station x (ni), strictly increasing */
    const f64 *y_lo, *y_hi;     /* y at j = 0 and j = nj-1 on each station (ni) */
    const f64 *z_lo, *z_hi;     /* z at k = 0 and k = nk-1 on each station (ni) */
    const f64 *s, *t;           /* distributions in j (nj) and k (nk), 0 to 1 */
    LaplaceBC  face_bcs[6];     /* i-min, i-max, j-min, j-max, k-min, k-max */
} FlowGrid;
```

### Locating and sampling

```c
Location   flow_locate(const FlowGrid* g, f64 x, f64 y, f64 z, CellHint* hint, CellCoord* out);
FlowSample flow_sample(const FlowField* f, const CellCoord* c);
Location   locate_and_sample(const FlowField* f, const FlowGrid* g, f64 x, f64 y, f64 z,
                             CellHint* hint, FlowSample* out);
```

Positions are in metres. `flow_locate`:

1. returns `Invalid` for a NaN coordinate;
2. reflects the point across any symmetry plane it lies beyond, recording
   the reflection in `CellCoord.sy` / `sz`;
3. for x outside the stations, reports the i-face it is beyond;
4. finds the station cell i and local coordinate ξ, interpolates the
   cross-section bounds linearly in ξ, and expresses y and z as fractions of
   them; a fraction below 0 or above 1 means the point is beyond a j- or
   k-face;
5. finds j and k in the s and t tables, and the local coordinates η and ζ.

Beyond a face, `CellCoord.face` names it and the result follows that face's
boundary condition: wall → `WallImpact`; inlet or outlet → `Exited`; anything
else → `Invalid`. A point exactly on a boundary is inside. The droplet
centre is the test point; droplet radii of at most a millimetre are
negligible against walls metres away.

Within a cell, `y = y_lo(ξ) + (y_hi(ξ) - y_lo(ξ)) * s(η)` with each factor
linear in its local coordinate, which is exactly the trilinear map of the
cell's corner nodes. The lookup is therefore the exact inverse of that map,
with no iteration.

`flow_sample` interpolates u, v, w, rho, p and T trilinearly, computing the 8
corner weights once per lookup (`cell_stencil`), and flips v (w) back for a
point that was reflected in y (z).

**Cell hint.** `CellHint` is optional and in/out (NULL for none; start each
droplet from `CELL_HINT_NONE`). It holds the cell of the droplet's previous
lookup. Each 1-D search starts at the hinted cell and gallops outward
(1, 2, 4, … cells) before bisecting: two comparisons when the point is still
in the hinted cell, O(log m) when it is m cells away. A stale or invalid hint
costs only time; results are identical with and without one. Measured on the
current file: 113 ns per locate-and-sample along a trajectory with the hint,
144 ns without, and about 440 ns for scattered random points, where memory
access dominates. The hint gains little on a separable grid; it is what keeps
the search local on a general grid.

**Alternative for non-separable grids** (non-rectangular sections, corner
fillets, stations not at constant x): invert each cell's trilinear map by
Newton iteration and walk to neighbouring cells starting from the hinted
cell. Roughly three times the code, with robustness work in skewed cells. The
interpolation is the same, so switching does not change results, and the
interface above stays the same.

### Symmetry and gravity

The flow is solved on a quarter of the cross-section (symmetry on j-min and
k-min). The flow is symmetric; droplet trajectories are not. Gravity acts
along −z and breaks the up–down symmetry across the k-min plane, and
supercooled large drops sag on the order of a metre through the contraction
(`docs/model-cleanup-and-port.md`, "Operating regime"). Droplets therefore
live in the full cross-section, and only the flow lookup is folded.

Gravity does not break the lateral symmetry, so with a symmetric j-min face
the release grid needs only y ≥ 0; results for y < 0 are mirror images. The
release grid must span the full height. With a full-domain flow file the
fold does nothing and nothing else changes.

## Droplets

```c
typedef struct Droplet {
    f64 x, y, z;                  /* m, full cross-section */
    f64 u, v, w;                  /* m/s */
    f64 d, d0;                    /* m */
    f64 T_shell[MAX_SHELLS];      /* K, centre to surface */
    f64 t;                        /* s since release */
    f64 We_max, Re_max;           /* diagnostics */
    DropletStatus status;         /* Active, Impacted, Evaporated, Exited, MaxSteps */
} Droplet;
```

### Time step

The exponential trajectory integrator and implicit shell solve are stable at
any step, so dt is set by accuracy only:

```
dt = min( C_cell * min(dx_i / |u_d|, dy_j / |v_d|, dz_k / |w_d|),   C_cell ≈ 0.2
          C_mass * m / |mdot| )                                       C_mass ≈ 0.01
```

where dx_i, dy_j, dz_k are the sizes of the droplet's current cell:
`dx_i = xs[i+1] - xs[i]`, `dy_j = (y_hi - y_lo) * (s[j+1] - s[j])` and
`dz_k = (z_hi - z_lo) * (t[k+1] - t[k])` at the droplet's x. The first term
keeps the droplet within a fraction of a cell per step, so the frozen air
state is accurate; the second bounds the mass change per step. With 100 mm
streamwise cells this is roughly 1.6e-4 s in the test section and a few
thousand steps per droplet for the whole tunnel. A step limit guards against
stalled trajectories and ends them as `MaxSteps`.

### Cost

About 10⁴ droplets × a few thousand steps × (a ~0.1 µs flow lookup + a
15-shell tridiagonal solve) is a few times 10⁷ steps — seconds on one core.
Performance does not shape the design:

- **Array of structs `Droplet[]`.** It mirrors the prototype's per-droplet
  state, which keeps the port easy to check against it. Struct-of-arrays and
  vectorisation are unnecessary at this cost.
- **Each droplet runs to completion before the next starts.** Trajectories end
  at different step counts and only per-droplet results are needed. The flow
  field is read-only and each droplet's scratch is private, so the loop can be
  threaded across droplets later without changing its shape.

### Memory

- `memory` arena, sized for the run: grid tables, `Droplet[]`,
  `PlaneRecord[]`, outcomes.
- `scratch` arena, reset per droplet with `arena_begin_temp` /
  `arena_end_temp`: the tridiagonal workspace inside `droplet_step`.
- Storage grows with n_droplets × n_planes, not with the number of steps.

## Sweep loop

```c
/* release grid: diameters x offsets in the spray-bar plane */
droplets = arena_push_array_zero(memory, Droplet,     n_d * n_y * n_z);
records  = arena_push_array_zero(memory, PlaneRecord, n_d * n_y * n_z * n_planes);

for each (diameter, y offset, z offset):
    initialize a Droplet at x_release: velocity, temperature and diameter per the release conditions

for each droplet:
    temp = arena_begin_temp(scratch)
    hint = CELL_HINT_NONE
    for step in 0 .. max_steps while droplet.status == Active:
        loc = locate_and_sample(&flow, &grid, droplet position, &hint, &sample)
        if loc == WallImpact: droplet.status = Impacted; record impact position; break
        if loc != Inside:     droplet.status = Exited;   break
        prev = droplet
        droplet.status = droplet_step(&droplet, &sample, &params, dt, scratch)
        for each sampling plane between prev.x and droplet.x: record crossing
        dt = next step size
    arena_end_temp(temp)

write the run log
```

## Outputs

**Plane records.** When a droplet crosses a sampling plane between two steps,
its state is interpolated linearly in time to the exact plane x. Each record
holds: y, z, u, v, w, d, d / d0, surface and mean droplet temperature, local
air speed and temperature, slip, We_max and Re_max so far, and time since
release. Default planes: the test-section entrance (x = 24.2 m in the current
file) and the model station.

**Outcome per droplet:** reached the outlet, impacted a wall (with impact
position and face), evaporated, or hit the step limit.

**Optional trajectory dump:** a full per-step history for a few selected
droplets, written as CSV for plotting. It is additive — one buffer sized for
a single droplet — and does not change the shape above.

**Run log:** human-readable and not parsed by any program, in the style of
Fortran formatted output — 132 columns, numbers in Ew.d format, sections
headed by banners. Sections:

1. Inputs: command line and configuration, program version, run start time.
2. Flow file: dimensions, extents, boundary conditions, inlet and outlet
   conditions.
3. Release conditions and model parameters.
4. Per-diameter summary at each sampling plane: fraction reaching the plane;
   mean, min and max of droplet temperature, d / d0, droplet-to-air speed
   ratio and vertical centroid shift; maximum Weber number.
5. Outcomes: fraction impacted and where, fraction evaporated, steps taken.
6. Run time.

The fixed-width record writer in `laplace-inspect` (fields placed at fixed
columns, Fortran A/I/E formatting, '*'-fill on overflow, output to several
streams) is the model for this log. Including the program version requires
the generated `version.h` directory on the target's include path.

## Open decisions (driver)

- **D1 — What the program reports.** Recommended: plane records and outcomes
  as above. Collection efficiency on a body applies only once a model body is
  in the flow solution; the current flow file is the empty tunnel. If a body
  is added, starting offsets will need clustering near the stagnation
  streamline.
- **D5 — Release conditions.** Needed: spray-bar x position; release grid
  (full height; y ≥ 0 with a symmetric j-min face); initial droplet speed
  (absolute, or as a ratio to the local air speed); initial droplet
  temperature; list of diameters. The prototype's start point (x = 11.45 m,
  air 60 m/s) does not correspond to a location in the current tunnel.

D2 (size range and breakup), D3 (freezing) and D4 (humidity input) are in
`docs/model-cleanup-and-port.md`.

## Validation (coupling and driver)

Physics tests are in `docs/model-cleanup-and-port.md`.

| Test | Checks | Reference |
|---|---|---|
| Sample at every grid node | interpolation | node values, to round-off |
| Interpolate the coordinate fields at a located point | exact inverse of the cell map | the point itself (after folding) |
| Sample at cell centres | interpolation | same cell at local (0.5, 0.5, 0.5); mean of the 8 corners |
| Sample at mirrored points | symmetry fold | identical scalars; v (w) negated for reflection in y (z) |
| Lookups with moving, stale and invalid hints | cell hint | identical to unhinted lookups |
| Points beyond each face, on boundaries, NaN | boundary outcomes | face and `Location` per boundary condition |
| Node moved off the separable position; symmetry face bent | grid check | refused at load |
| 1 µm tracer | lookup and integrator together | slip → 0 along the path |
| Releases at (y, z) and (−y, z) | symmetric handling | mirror-image trajectories (and in z, with gravity off) |
| Uniform air in a straight duct | plane-crossing interpolation | analytic trajectory |
| Centreline release | whole pipeline | 1-D centreline estimate in `docs/model-cleanup-and-port.md` |

The first seven rows pass against the current flow file (round-trip error
6e-16 of the domain size, node values to 3e-16); the test program is not yet
part of the repository.

## Implementation plan

Each phase ends when its tests pass.

0. **Decisions.** Settle D1–D5 and record them in these two documents.
1. **Corrected Python prototype.** Physics and numerics fixes, then the
   physics tests in Python; save reference outputs.
   (`docs/model-cleanup-and-port.md`, "Sequencing")
2. **C physics kernel.** Closures and `droplet_step`; match the stored Python
   outputs; physics tests in C.
3. **Flow coupling — done.** Grid tables, separability check,
   `flow_locate` / `flow_sample` with symmetry folding and the cell hint.
   Remaining: move the lookup tests into the repository.
4. **Driver.** Release grid, sweep loop, adaptive dt, plane records,
   outcomes; the tracer, mirror, straight-duct and centreline tests.
5. **Output.** The fixed-width run log and the optional trajectory dump.
6. **Studies.** Diameter sweep and spray-bar position sensitivity; then
   decide D2 from the Weber numbers of 3-D trajectories.

## Upstream items in `laplace`

These affect the flow file this program reads. None blocks phases 1–2.

- **Axis convention.** This program takes y (j) as lateral and z (k) as
  vertical. `laplace_grid` builds y from the wall half-height and z from the
  half-width (`src/laplace_core.f90`), and `laplace.h`, the Fortran docs and
  the README call j vertical and k spanwise, while the C application's
  options and comments call j lateral and k vertical. Square sections hide
  the difference; for rectangular ones, height and width arrive swapped.
  Proposed fix: make the generator axis-neutral — it takes the wall
  half-extent in y and in z instead of a half-height and half-width — and
  state the convention once in the README and `laplace.h`; the application
  passes the half-width as the y extent and the half-height as the z extent.
- **Velocity across symmetry planes is not zero.** Node values of v on j-min
  and w on k-min reach 0.18 m/s in the current file. After folding, droplets
  near those planes see a small spurious push across them. Setting v = 0 on
  j-min and w = 0 on k-min in `laplace_output` removes it.
- **Mass flow is not conserved along the contraction.** The area-integrated
  mass flow deviates from the inlet value by up to −1.3% (i = 547,
  x = 18.7 m), although the solve is converged (max residual 5.6e-7 against a
  tolerance of 1e-6). Candidates: the node-velocity reconstruction in
  `laplace_output`, or discretisation error. A mass-flow check from the
  solver's own face fluxes and a grid-refinement study distinguish the two.
  For droplets this is a local air-speed error of about 1%.
- **Format version.** The header layout now carries `header_size`, but
  `LAPLACE_VERSION` is still 1, so readers cannot tell old and new files
  apart. Bump it to 2.
- **Flow-file reading is implemented twice**, in `flow_field_open` here and
  in `laplace-inspect`, and the writer in the `laplace` C application keeps
  `n_fields`, the field names and the write calls as three lists that must
  agree by hand; a wrong write order produces a valid-looking file with
  swapped fields. Proposed:
  - Split `laplace.h` into the kernel ABI (`laplace.h`; users link
    `laplacelib`) and the file format (`laplace_file.h`: header struct,
    constants, boundary-condition enum, `laplace_idx`; `stdint.h` and
    `string.h` only, never linked).
  - Add pure functions on the header to `laplace_file.h`, with no file I/O
    (callers keep their own mapping and streaming):
    `laplace_header_init` fills magic, version, `header_size`, `n_fields` and
    the names from one `{name, array}` table, which the writer also uses for
    its write loop; `laplace_header_check(h, file_size)` returns 0 or an
    error code; `laplace_field_offset(h, name)` returns a byte offset or −1.
  - Fortran cannot call C `static inline` functions. If a Fortran program
    ever writes flow files, move the bodies into a small compiled C library
    with the same declarations.

## Build notes

- `FetchContent` with `SOURCE_SUBDIR` (used to take only the `laplace`
  headers) requires CMake 3.18; `cmake_minimum_required` currently says 3.14.
- `version.h` is generated into the build's `generated` directory, which is
  not on the target's include path yet.
