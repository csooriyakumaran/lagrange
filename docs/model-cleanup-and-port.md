# Droplet Model: Review and Porting Plan

Source reviewed: `prototype/Airbus droplet model 2.py`, a single-droplet
Lagrangian trajectory, heat-transfer and mass-transfer model: drag, gravity,
buoyancy, convective heating with radial shell conduction inside the droplet,
and Sherwood-number evaporation. Target: a C kernel, `droplet_step`, called by
the trajectory driver in `docs/driver-design.md`, which supplies local air
conditions from the `laplace` potential-flow solution of the icing wind
tunnel.

Open decisions are labelled D1–D5. D2–D4 concern the physics and are listed
at the end of this document; D1 and D5 concern the driver and are in
`docs/driver-design.md`, which also holds the full implementation order.

## Operating regime in the target tunnel

The current flow file (`data/flow.bin`) is a large closed-section icing wind
tunnel with a square cross-section:

| Region | x (m) | Cross-section | Centreline air speed | Static T |
|---|---|---|---|---|
| Settling chamber | −36 to ≈ 0 | 9.0 m × 9.0 m | 11.7 m/s | −10.1 °C |
| Contraction, 10.2 : 1 | ≈ 0 to ≈ 24 | narrows to 2.83 m × 2.83 m | 11.7 → 124 m/s | → −17.6 °C |
| Test section | ≈ 24 to 47.6 | 2.82 m × 2.82 m | 128 m/s (M 0.40) | −18.2 °C |

Stagnation conditions are 101 325 Pa and 263.15 K; static pressure falls from
101.2 kPa in the settling chamber to 90.7 kPa in the test section.

### Size dependence: 1-D centreline estimate

Method: one droplet released on the tunnel axis at x = 0 m with the local air
speed, tracked to the test-section entrance at x = 24.2 m. Air speed, density
and temperature are interpolated along the axis nodes of the flow file;
viscosity from Sutherland's law; sphere drag (Schiller–Naumann, Cd = 0.44
above Re = 1000) applied to the full slip vector; gravity acts while the axis
air has no vertical velocity; no evaporation. Weber number
We = ρ_air |slip|² d / σ with σ = 0.0757 N/m. Transit time is 1.2–1.4 s for
all sizes.

| d (µm) | Stokes τ (s) | max slip (m/s) | max Re | max We | u_drop / u_air at x = 24.2 m | sag at x = 24.2 m |
|---:|---:|---:|---:|---:|---:|---:|
| 20   | 0.0013 | 1.7  | 3    | 0.0  | 0.99 | 1.5 cm |
| 100  | 0.033  | 11.4 | 88   | 0.2  | 0.93 | 28 cm  |
| 300  | 0.30   | 25.9 | 600  | 3.3  | 0.81 | 1.2 m  |
| 600  | 1.2    | 36.5 | 1688 | 13.3 | 0.72 | 2.4 m  |
| 1000 | 3.3    | 43.6 | 3357 | 31.5 | 0.66 | 3.7 m  |

Maximum slip occurs at x ≈ 22–23 m, where the air accelerates fastest. The
sag column is an upper bound: away from the axis, the converging flow in the
contraction has a vertical component that pushes droplets back toward the
centreline.

What this means for the model:

- **Icing-cloud sizes (tens of µm, 14 CFR 25 Appendix C).** We < 1, droplets
  stay within about 7% of the air speed, and gravity drift is centimetres to
  tens of centimetres. A rigid-sphere model is valid.
- **Supercooled large drops (hundreds of µm, Appendix O).** We reaches 3 at
  300 µm, 13 at 600 µm and 32 at 1 mm, against an onset of bag breakup near
  We ≈ 12. Drops flatten and their drag rises before they break, so
  rigid-sphere drag is outside its range from a few hundred µm up in this
  contraction. These drops also lag the air by 20–35% at the test section and
  drift downward on the order of a metre, comparable to the test-section
  half-height of 1.41 m.
- **Gravity breaks the up–down symmetry of the flow solution.** Trajectories
  must be computed over the full height of the tunnel even when the flow is
  solved on a quarter domain (see "Symmetry" in `docs/driver-design.md`).

### The prototype's test case

The prototype runs one 595.6 µm droplet starting at x = 11.45 m in uniform air
at 60 m/s, −13 °C, 101 325 Pa and 95% relative humidity, with the droplet at
51 m/s and −3.5 °C; 15 shells, dt = 1e-5 s, 30 000 steps. These conditions
do not correspond to a location in the target tunnel (the centreline air
speed at x = 11.45 m is about 18 m/s); release conditions for coupled runs
are decision D5.

## What the model gets right

- Drag correlation: Schiller–Naumann with Cd floored at 0.44, continuous at
  Re = 1000.
- Property correlations: Kell equation for water density, Sutherland's law
  for air viscosity, Ranz–Marshall for the Sherwood number, Magnus formula
  for saturation vapour pressure over (supercooled) liquid water — all
  standard forms, correctly transcribed.
- State is threaded explicitly through the time loop with no hidden globals.
  The environment queries (air velocity, pressure, temperature, humidity) are
  functions of position, which is the right shape for coupling to an external
  flow field.
- One-way coupling (droplets do not perturb the air) is correct for dilute
  droplet clouds and should be kept.

## Physics issues to resolve before porting

### 1. Drag is applied per component

Each axis uses `0.5 * rho * |v_rel_x| * v_rel_x * A * Cd`, so the cross-stream
force scales with the cross-stream slip alone. The correct force is along the
slip vector with magnitude set by the total slip:

```c
f64 s = sqrt(sx*sx + sy*sy + sz*sz);              /* |V_air - V_drop| */
f64 k = 0.5 * rho_air * s * A * Cd / m;            /* one factor for all three axes */
ax = k * sx;   ay = k * sy + g_eff;   az = k * sz;
```

The per-component form makes cross-stream drag far too weak whenever the
slip is mostly streamwise, which is always the case in a contraction. With
it, the centreline estimate above gives a 600 µm droplet a vertical speed of
3.3 m/s at the test-section entrance instead of 0.56 m/s, and 26% more sag
(58% more at 100 µm). Converging streamlines in the contraction are followed
less closely, so wall impacts are overpredicted. The Stokes branch (Re < 1)
is linear and unaffected.

### 2. Evaporation is decoupled from the energy balance and can change sign

Mass loss is computed two ways — from the humidity difference, and from a
heat-limited rate `h * pi * d^2 * (T_air - T_drop) / L` — and the one with the
smaller magnitude is kept. The latent heat of the evaporated mass is never
removed from the droplet's energy.

The heat-limited rate is negative whenever the droplet is warmer than the
air. In the target tunnel that is the normal case (droplet −3.5 °C, air
−10 °C to −18 °C), so the smaller-magnitude rule picks the negative value
and the droplet gains mass while it should be evaporating.

Fix: drop the heat-limited estimate. The evaporation rate comes from mass
transfer alone — Ranz–Marshall Sherwood number and the vapour-density
difference between the surface (saturated at the surface temperature) and
the local air — and `mdot * L_v(T_surface)` is a heat sink in the outer
shell's energy balance.

A Spalding-number film model (Abramzon–Sirignano) is not needed at these
conditions. The vapour mass fraction at the surface is about 3e-3 and in the
air about 2e-3, so the mass-transfer number B_M ≈ 1e-3 and
ln(1 + B_M) = B_M to within 0.06%: the classical film model above is already
accurate. Revisit only if warm or low-pressure cases are added.

### 3. Coordinate convention

x is streamwise, y (grid index j) is lateral, and z (grid index k) is
vertical, positive up. The prototype applies gravity along −z and buoyancy
along +z, which matches, so this carries over unchanged. Flow files must use
the same convention (`docs/driver-design.md`, "Upstream items in
`laplace`").

### 4. Relative humidity is not conserved through the contraction

The prototype holds relative humidity fixed at each position. Through the
contraction no water condenses or is added, so the conserved quantity is the
humidity ratio

```
w   = 0.622 * p_v / (p - p_v)             (set once at a reference station)
p_v = w * p / (0.622 + w)                 (local vapour pressure, from local static p)
```

while the saturation pressure falls with the static temperature. In the
target tunnel the relative humidity at the test section is 1.74 times the
settling-chamber value: 95% in the settling chamber becomes 166% (supersaturated
over liquid water) at the test section. Humidity must therefore be specified
at a named station and carried as a humidity ratio (decision D4).

### 5. The Nusselt correlation has no identified source

The Sherwood number is textbook Ranz–Marshall
(`Sh = 2 + 0.6 Re^(1/2) Sc^(1/3)`), but the Nusselt number uses an unfamiliar
expression. It tends to Nu = 2 as Re → 0 and to about 0.58 (Re Pr)^(1/2) at
large Re Pr (for air, Pr ≈ 0.71), so it is plausible, but it has no
identified source. Use
Ranz–Marshall for both (`Nu = 2 + 0.6 Re^(1/2) Pr^(1/3)`), which keeps heat
and mass transfer consistent by the heat/mass analogy; keep the old
expression only as a comparison.

### 6. No freezing physics

The droplet starts supercooled but only liquid-phase heat and mass transfer
are modelled; there is no nucleation, latent heat of fusion or ice phase.
This is standard if the model's job ends at the test section and accretion
physics belongs elsewhere (decision D3).

### 7. No deformation or breakup

Drag assumes a rigid sphere at all Weber numbers. Valid for icing-cloud
sizes; not for supercooled large drops in this contraction (see the regime
table). Decision D2.

### 8. Smaller issues

- The Stokes branch below Re = 1 jumps the drag by 15% at Re = 1 (the
  Schiller–Naumann factor there is 1.15), and `Cd = 24 / Re` is evaluated
  before the branch, which divides by zero at zero slip. Writing drag as
  `F = 3 pi mu d f(Re) * V_rel` with `f(Re) = Cd Re / 24`
  (`1 + 0.15 Re^0.687` below Re = 1000, `0.44 Re / 24` above) is continuous
  and finite at Re = 0, and removes the branch.
- Latent heat is the 100 °C value, 2.26e6 J/kg. Near −10 °C it is about
  2.52e6 J/kg (11% higher); use `L_v(T) ≈ 2.501e6 - 2370 T[°C]` J/kg.
- The heat-transfer branch is chosen from the Biot number at the first step,
  but the shell model is still advanced every step even when its result is
  discarded; with an explicit scheme it can go unstable unseen. The implicit
  shell model below removes the branch entirely.
- The prototype computes a humid-air density for Re and drag. In the coupled
  model the air density comes from the flow file (dry air); at these
  humidity ratios the difference is about 0.1%.
- Comment fixes, so they are not "corrected" into bugs during the port: the
  gas constant of air is labelled `kJ/kg*K` but the value (287.05) is in
  J/(kg K), which is correct; the viscosity docstring says kg/m3, the value
  is in Pa s.

## Numerical method

### Trajectory: exponential integrator

Over one step, hold the air velocity, the drag relaxation time
`tau = rho_p d^2 / (18 mu f(Re))` and the gravity-plus-buoyancy acceleration
`g_eff` constant. The equation of motion is then linear and integrates
exactly:

```c
f64 e = exp(-dt / tau);
v_new = v_air + g_eff*tau + (v - v_air - g_eff*tau) * e;
x_new = x + (v_air + g_eff*tau)*dt + (v - v_air - g_eff*tau) * tau * (1 - e);
```

This is stable for any dt, so the drag stiffness of small droplets (τ ≈ 1 ms
at 20 µm) never limits the step. The only error is from holding the air
state and τ fixed over the step, which the driver controls by limiting how
far the droplet moves per step (`docs/driver-design.md`, "Time step"). The
prototype's explicit Euler update with dt = 1e-5 s needs about 30 000 steps
for 0.3 s of flight; the exponential update needs a few thousand for the
whole tunnel.

### Internal conduction: implicit shells

The prototype advances shell conduction explicitly. Explicit diffusion on
equal-thickness spherical shells is stable only for `dt <~ dr^2 / (3 alpha)`
(the centre shell is the binding one), with
`alpha = k_water / (rho_water cp_water) ≈ 1.3e-7 m^2/s`:

| d | shells | explicit limit |
|---|---|---|
| 600 µm | 15 | ≈ 1 ms |
| 20 µm | 15 | ≈ 1 µs |

The prototype's dt = 1e-5 s happens to be stable for its 600 µm droplet and
is unstable for icing-cloud sizes.

Use backward Euler with a tridiagonal (Thomas) solve: unconditionally stable,
O(N) per step, no oscillation when a sudden surface flux meets a large dt
(Crank–Nicolson rings there). First-order time accuracy is sufficient
because dt stays far below the droplet's thermal time,
`rho_p cp_p d^2 / (6 Nu k_air)` — about 6 ms at 20 µm and about 1 s at 600 µm.

With an implicit solve, a single shell is the lumped model, so the Biot-number
branch disappears: small-Biot droplets simply come out nearly uniform. The
shell count becomes a resolution parameter, checked by convergence.

### Heat and mass coupling within a step

Evaluate the evaporation rate at the start of the step from the current
surface temperature, then apply convection and the latent-heat sink as the
outer boundary condition of the implicit shell solve. Lagging the rate by one
step is first-order consistent and adequate because dt is far below the
thermal time; the wet-bulb test below checks the coupling.

### Diameter update

Update mass from the evaporation rate, then diameter from mass and the
mean-temperature water density. Shells are re-spaced to the new radius each
step with temperatures carried by index; the per-step change in radius is
tiny, so the remapping error is negligible.

## Kernel shape

- The ~30 scalar closures (densities, viscosity, conductivities, Re, f(Re),
  Nu, Sh, saturation pressure, L_v) are pure functions with no state. They
  port almost verbatim.
- One update per step:

  ```c
  DropletStatus droplet_step(Droplet* d, const FlowSample* air, const DropletParams* prm,
                             f64 dt, Arena* scratch);
  ```

  `FlowSample` is the interpolated air state (u, v, w, rho, p, T) from the
  driver. `DropletParams` holds run constants: humidity ratio, gravity
  vector, shell count, surface tension. `scratch` provides the tridiagonal
  workspace.
- Shell temperatures are droplet state that persists across steps, so they
  live in `Droplet` (a fixed-capacity array keeps `Droplet` a plain value).
- The kernel also updates per-droplet diagnostics the driver reports: the
  maximum Weber and Reynolds numbers seen so far.
- Array of structs, one droplet at a time (`docs/driver-design.md`, "Cost").

## Validation

Each test has an independent reference; run them in Python on the corrected
prototype first, then in C.

| Test | Checks | Reference |
|---|---|---|
| Closures against the corrected Python prototype | the port | stored prototype outputs for a set of states |
| Relaxation in uniform air, no gravity | trajectory integrator | exact exponential (the integrator should match to round-off) |
| Terminal velocity in still air | drag and gravity | Gunn & Kinzer (1949), 20 °C and 1013 hPa: 0.27, 2.06, 4.03 m/s at 0.1, 0.5, 1.0 mm |
| Sphere with a convective surface, no evaporation | shell solver; convergence order in dr and dt | Carslaw & Jaeger series solution at Bi = 0.1, 1, 10 |
| Single shell | lumped limit | exponential cooling with constant h |
| Evaporation in still air at constant temperature | mass transfer | d² law: d² falls linearly with time |
| Steady droplet temperature in uniform air | heat/mass coupling (issue 2) | psychrometric wet-bulb temperature |

## Open decisions (physics)

- **D2 — Droplet size range.** Recommended: build the rigid-sphere model
  first and report the maximum Weber number per trajectory. Add
  deformation-dependent drag and a breakup criterion only if supercooled
  large drop runs are needed, once 3-D trajectories give actual Weber
  numbers.
- **D3 — Freezing.** Recommended: out of scope; report droplet temperature at
  each sampling plane so supercooling is visible.
- **D4 — Humidity input.** Recommended: relative humidity (or humidity ratio)
  at a named reference station, converted once to a humidity ratio and held
  constant.

## Sequencing

Fix the physics and numerics in the Python prototype first, where iteration
is cheap; porting a numerically fragile model only moves the fragility into
a language where it is harder to debug. These are phases 1 and 2 of the
implementation plan in `docs/driver-design.md`.

**Phase 1 — corrected prototype.**

1. Vector drag with `f(Re)`; gravity and buoyancy stay on z.
2. Evaporation from mass transfer only, with the latent-heat sink;
   temperature-dependent `L_v`.
3. Ranz–Marshall Nusselt number.
4. Implicit shell conduction; one shell as the lumped model; remove the Biot
   branch.
5. Exponential trajectory integrator.
6. Humidity carried as a humidity ratio.

Done when every test in the validation table passes in Python, and reference
outputs for the closures and a few full trajectories are saved for the C
port.

**Phase 2 — C kernel.** Closures and `droplet_step` as above. Done when the
closures match the stored Python outputs to near round-off and the
validation table passes in C.
