# FABLE Plan — `fecore`: a Compatible Finite-Element Dynamical Core
## Implementation Edition

| | |
|---|---|
| **Version** | v1.0 — implementation edition (supersedes v0.1 roadmap; same mission, now specified to hand-off depth) |
| **Date** | 2026-07-12 |
| **Author** | Claude Fable 5 (high-level planning), with David Fillmore |
| **Audience** | Implementing coding agents (Sonnet/Haiku-class) executing one task per session, plus David as gatekeeper |
| **Status** | Sub-project 1 COMPLETE (all gates GO). SP2–SP7 + XC specified below. **No task begins without David's explicit GO.** |
| **Scope note** | Tracked in this fork by David's decision (2026-07-12) so implementing agents on any host see the plan of record. Task briefs/reports remain in `.git/sdd/` (host-local, never committed). |

---

## 0. How to use this plan (instructions to implementing agents)

You are implementing **one task** (`SPx-Ty` or `XC-y`) from §7. Protocol:

1. Read §0–§6 of this file, your task's spec in §7, and Appendices A–D.
   Read `fecore/operator/poisson.py`'s module docstring (dolfinx API notes)
   and Appendix C before writing any dolfinx code.
2. `conda activate fecore`. Run `python -m pytest fecore/tests/ -q` FIRST —
   the suite must be green (plus the recorded xfails) before you start.
   If it isn't, STOP and report.
3. Write a brief to `.git/sdd/task-SPx-Ty-brief.md`: your reading of the task,
   assumptions, API points you will verify. (Task-1 precedent: briefs get API
   details wrong; that's fine — record the corrections in the report.)
4. Implement. Rules:
   - Touch only files your task owns (§7 lists them). If you must touch a
     shared file, say so in the report.
   - Every public function gets a docstring in the style of
     `fecore/operator/poisson.py` (NumPy-style params, returns, notes).
   - New knowledge about dolfinx quirks goes into the relevant module
     docstring (the SP1 pattern) AND your report.
   - At most ONE environment change per task; if you add a package, pin it in
     `fecore/environment.yml` in the same commit and state why.
   - **Logged, never silent**: any skipped case, fallback, or capped
     computation must emit a `log.warning` saying what was dropped and why.
5. Tests: implement the named tests in your task spec (you may add more).
   Gate thresholds are law — you may NOT weaken a threshold to pass. If a
   gate fails and you can characterize the failure, mark the test `xfail`
   with a reason string and record a **NO-GO** verdict; if you cannot
   characterize it, record **BLOCKED**. Never delete a failing gate test.
6. Run the FULL suite. Record pass/xfail counts.
7. Write `.git/sdd/task-SPx-Ty-report.md` using the template in Appendix D
   (goal, what was done, API surprises, gate numbers, verdict GO/NO-GO/BLOCKED,
   any corrections to this plan → list under "Plan feedback").
8. Commit: one task ≈ one commit, message `fecore SPx-Ty: <title>` with body
   summarizing gate results, trailer `Co-Authored-By: <your model name>
   <noreply@anthropic.com>`. Do NOT push.
9. **STOP.** Report the verdict. The next task starts only on David's GO.

**Gate classes** used throughout §7:
- **[S] Structural** — exact/mathematical: conservation to 1e-12 (relative),
  conformity, mode counts, convergence slopes (≥1.9 second-order claims,
  ≥1.7 for nonlinear end-to-end, ≥0.9 lowest-order). Hard thresholds, fixed.
- **[B] Benchmark** — vs pinned reference data or published values. Values
  marked **CALIBRATE** are set by the first GO run (recorded in the report
  and then backfilled into the test as a regression band ±20% unless stated).
  Values marked **⚠verify** are from the planning model's memory of the
  literature: re-derive them from the cited paper before hard-coding.

**Numbers discipline:** any physical constant or case parameter comes from
`fecore/constants.py` / the case module (App. A/B) — never inline literals.

---

## 1. Mission and context (condensed)

Build **`fecore`**: a compatible (mimetic) finite-element dynamical-core
prototype in Python/FEniCSx on MPAS SCVT-derived meshes, culminating in a
working 3-D moist nonhydrostatic dycore verified against MPAS-A on its own
idealized cases. Motivation: the MPAS FV Poisson operator converges at L²
slope **1.380** on the official SCVT bundles with error pinned to pentagon
rings (`docs/design/poisson/verification-plan/tier-a-mms.md`); SP1 proved FEM
restores clean 2nd order **on the same meshes**. Compatible FE (the
LFRic/Gusto/FEEC lineage: Cotter & Shipton 2012) gives TRiSK's mimetic
properties — local mass conservation, no spurious pressure modes, steady
geostrophic modes — **without grid orthogonality**, which is exactly what the
pentagons break.

**Definition of Done (dycore v1.0):** dry JW baroclinic wave at ~120 km and
the reduced-radius moist supercell running stably in MPI (≥16 ranks), all [S]
gates green, benchmark cases within pinned bands, side-by-side FE-vs-FV
reports for the three `~/Data/MPAS` cases, RTD design narrative current,
tagged release + Zenodo DOI. (Production port is explicitly out of scope —
SP7-T4 is a decision memo.)

Non-goals: MPAS Fortran integration, full physics suite, DA, performance
parity (we *measure* honestly instead).

---

## 2. Current state — inherited assets and invariants

Sub-project 1 (Tasks 1–8, branch `feature/fecore`) delivered, all gates GO,
**30 passed / 1 xfail**:

| Asset | File | Notes |
|---|---|---|
| Icosa surface + shell meshes | `fecore/mesh/icosa.py` | `icosa_surface_mesh(refine)`, `icosa_shell_mesh(refine, n_layers, H)`; unit sphere, outward-oriented |
| SCVT-dual surface + shell | `fecore/mesh/scvt_dual.py` | from MPAS `grid.nc` (`cellsOnVertex` triangles on cell centres); `sphere_radius(path)` |
| Extrusion utils | `fecore/mesh/_shell_utils.py` | `_extrude_to_tets`, `_fix_tet_orientations` (Freudenthal–Kuhn, conforming) |
| Elliptic solvers | `fecore/operator/poisson.py` | `solve_surface` (gauge-pinned), `solve_shell` (ground-Dirichlet/top-Neumann); P1/P2; consistent+lumped RHS; **module docstring = dolfinx 0.10 API notes** |
| MMS | `fecore/verify/mms.py` | Y₄² surface; Y₄²·sin(πζ/2) shell |
| Convergence harness | `fecore/verify/convergence.py` | finest-two log-log slope; **coupled vertical/horizontal refinement**; SCVT bundle registry at `~/Data/MPAS/poisson_tier_A2_scvt/meshes/{480,240,120,60}km/grid.nc` |
| Environment | `fecore/environment.yml` | dolfinx 0.10.0, petsc4py 3.24, mpi4py 4.1, MPICH, numpy 2.4, netcdf4 1.7, pytest 9.1, Python 3.14. **No scipy yet.** |

Key SP1 numbers (finest-two L² slopes): surface icosa/SCVT P1-consistent
1.9936/1.9994, P2 2.0083/2.0000, P1-lumped ≈2.00 both; shell icosa/SCVT P1
1.9699/1.9743, P2 2.0187/2.0546. FV reference: 1.380 (L∞ 0.253).

Inherited invariants (do not break):
- The xfail anchor `test_scvt_proper_fem_beats_the_defect`'s companion
  (lumped-FEM-also-converges) stays xfail — it is a scientific record.
- Degree-1 geometry caps ALL elements at O(h²) (Dziuk 1988) — do not chase
  3rd-order slopes before XC-A.
- Physical thin shell (H/R≈3e-3) is ill-conditioned with the SP1 approach —
  normalized geometry (R=1, H=0.25) is the standard test frame until XC-B.
- MPAS `grid.nc` conventions: 1-based indices (subtract 1), `cellsOnVertex`
  rows with any 0 ⇒ incomplete triangle (drop), pentagons = cells with
  `nEdgesOnCell == 5` (exactly 12 on quasi-uniform meshes).

---

## 3. Mathematical formulation compendium

Implementers: this section is the single source for the continuous and weak
forms. Task specs reference equations as (3.x). UFL-style pseudocode uses
`inner, div, grad, curl, cross, jump, avg, dx, dS`; `n_f` = facet normal,
`n` = outward unit sphere normal ≈ `x/R` (use the **analytic**
`SpatialCoordinate(mesh)/R`, not `CellNormal` — orientation-safe, O(h²)
consistent with flat cells).

### 3.0 Notation and spaces

Discrete de Rham complexes (lowest order shown; degrees rise together):

- **LANE-T** (triangles, SCVT-dual + icosa surfaces):
  `Vq = CG2  --n×grad-->  V2 = BDFM1  --div-->  V3 = DG1`
- **LANE-Q** (quads from SCVT-primal/icosa splitting):
  `Vq = Q1   --n×grad-->  V2 = RTCF1  --div-->  V3 = DG0`
- **3-D (LANE-Q hexes)**: `W0=Q1 → W1=NCE1 → W2=NCF1 → W3=DG0`, plus
  `Vθ` = scalar, horizontally-DG0, vertically-CG1 (Charney–Phillips: θ DOFs
  co-located with vertical-velocity DOFs on horizontal facets).
  (basix exposes RTCF/NCF/NCE as the RT / N1curl families on quad/hex cells —
  confirm exact enum names in K2.)
- Slice (2-D x–z quads): `V2 = RTCF1, V3 = DG0, Vθ` as above.

`u` velocity (V2/W2), `D` SWE depth (V3), `ρ` density (W3/V3), `θ` potential
temperature (Vθ), `Π` Exner (same space as ρ), `q` SWE PV (Vq), `b` bathymetry
/ orography, `f = 2Ω z/R` Coriolis parameter on the sphere.

`perp(u) := cross(n, u)` (tangential rotate-by-90°). Structural identity to
unit-test: `perp(perp(u)) = −u` for tangential u.

### 3.1 Elliptic (SP1 recap + mixed form)

Primal (done): `eps*inner(grad(u),grad(v))*dx = f*v*dx`.

Mixed (SP2-T1): find `(σ, φ) ∈ V2 × V3`:

```
(3.1a)  inner(σ, τ)*dx − φ*div(τ)*dx = 0                ∀ τ ∈ V2
(3.1b)  div(σ)*v*dx                  = f*v*dx           ∀ v ∈ V3
```

Closed surface ⇒ compatibility `∫f = 0` (Y₄² satisfies it) and φ defined up
to a constant: attach a PETSc nullspace (constant on the φ block) or pin one
φ DOF; align means before error norms (SP1 pattern). Local conservation
identity to test: for every cell K, `∫_K div σ_h = ∫_K f` holds to solver
tolerance because V3 test functions contain indicator functions.

### 3.2 Rotating shallow water

Continuous, vector-invariant:

```
(3.2a)  ∂u/∂t + q F⊥ + grad(g(D+b) + ½|u|²) = 0
(3.2b)  ∂D/∂t + div(F) = 0,   F = D u,   q = (ζ + f)/D,   ζ = curl_n(u)
```

Weak forms (Shipton–Gibson–Cotter 2018 lineage):

**PV diagnosis** — find `q ∈ Vq`:
```
(3.2c)  q*γ*D*dx = −inner(cross(n, grad(γ)), u)*dx + f*γ*dx     ∀ γ ∈ Vq
```
Sign convention MUST be fixed by the structural test: for solid-body rotation
`u = Ωs cross(ẑ_vec, x)` (tangential on the sphere), the diagnosed relative
vorticity is `2 Ωs z/R`.

**Mass flux projection** — `F ∈ V2`: `inner(F, w)*dx = inner(D*u, w)*dx`.

**Momentum** — ∀ `w ∈ V2`:
```
(3.2d)  inner(w, ∂u/∂t)*dx + q̃ * inner(w, perp(F))*dx
        − div(w)*(g*(D+b) + 0.5*inner(u,u))*dx = 0
```
with APVM stabilized PV `q̃ = q − τ_q * inner(u, grad(q))`, `τ_q = 0.5*Δt`.

**Continuity** — upwind DG, ∀ `φ ∈ V3` (no exterior boundary on the sphere):
```
(3.2e)  φ*∂D/∂t*dx − inner(grad(φ), u)*D*dx
        + jump(φ)*(un('+')*D('+') − un('−')*D('−'))*dS = 0
        where un = 0.5*(dot(u, n_f) + |dot(u, n_f)|)
```
(`dot(u,n_f)` is single-valued for H(div) u — this is what makes the upwind
flux well-defined.) LANE-Q DG0 note: the `grad(φ)` volume term vanishes; the
facet terms carry everything.

**Linearized SWE** (SP2-T5), about rest depth H:
```
(3.2f)  ∂u/∂t + f perp(u) + g grad(h) = 0 ;  ∂h/∂t + H div(u) = 0
```
weak: `inner(w,u_t) + f*inner(w, perp(u)) − g*h*div(w) = 0`;
`φ*h_t + H*φ*div(u) = 0`. Time: implicit midpoint (energy-conserving).
Discrete linear geostrophic balance solve for tests: given divergence-free
`u_h = Π_{V2} perp(grad ψ)`, find `h_h`: `g*h*div(w)*dx = f*inner(w, perp(u_h))*dx`
(+ mean-zero constraint).

### 3.3 Compressible Euler (slice and 3-D)

Variables `(u, ρ, θ)`, Exner from state:
```
(3.3a)  Π = (R_d ρ θ / p00)^(R_d/c_v)
(3.3b)  ∂u/∂t + (vector-invariant advection) + c_p θ grad(Π) + g ẑ + 2Ω×u = 0
(3.3c)  ∂ρ/∂t + div(ρ u) = 0
(3.3d)  ∂θ/∂t + u·grad(θ) = 0
```
(slice: no Coriolis; ẑ = Cartesian vertical. 3-D shell: ẑ → n, gravity g n,
geometry-exact "deep" — note O(H/R) vs shallow-atmosphere MPAS-A when
comparing.)

**Pressure-gradient weak form** (θ ∉ H1 horizontally — facet correction on
vertical/interior facets):
```
(3.3e)  PG(w; θ, Π) = − c_p*div(θ*w)*Π*dx
                      + c_p*jump(θ*w, n_f)*avg(Π)*dS
```
Boundary terms vanish on slip walls (`dot(w,n)=0` essential on V2 top/bottom;
lateral periodic or sponge). Structural test: with hydrostatically balanced
(θ*, Π*) (§3.5), `PG(w;θ*,Π*) + ⟨w, g ẑ⟩` residual is small and h-converging.

**Velocity advection** (vector-invariant, slice): diagnose scalar
`ζ ∈ CG1`: `ζ*γ*dx = −inner(perp2(grad(γ)), u)*dx` with
`perp2((a,b)) = (−b,a)`; then advection term
`ζ*inner(w, perp2(u))*dx − div(w)*0.5*inner(u,u)*dx`.
3-D: diagnose `ω ∈ W1`: `inner(ω,η)*dx = inner(u, curl(η))*dx` (+BC terms),
advection `inner(w, cross(ω, u))*dx − div(w)*K*dx`. Upwinded/SUPG variants
are allowed if the plain form is noisy — record in report (plan feedback).

**θ transport** (Vθ, vertically CG / horizontally DG): SUPG in the vertical,
upwind jumps horizontally:
```
(3.3f)  (χ + τ_θ inner(u, grad(χ))) * (θ_t + inner(u, grad(θ))) * dx
        + horizontal-facet upwind jump terms for the DG direction
```
`τ_θ = 0.5*Δt` default. **ρ transport**: (3.2e) with D→ρ.

**Semi-implicit quasi-Newton timestep** (LFRic/Gusto pattern; SP4-T2, SP5-T2):
θ-method with off-centering α (default 0.55 ⚠tunable). Each step: ~4 Picard
iterations solving the linearized system about time-lagged reference state
`(u*, ρ*, θ*)` for increments `(Δu, Δρ, Δθ)`:

1. **Δθ elimination (pointwise in columns):** Δθ = −α Δt (Δu·ẑ) ∂θ*/∂z at
   the shared Vθ/vertical-velocity DOF sites (this is WHY Vθ is built
   vertically-CG1 horizontally-DG0 — same DOF locations as W2's horizontal
   facets; verify co-location in K3 and exploit nodally).
2. **Δρ elimination:** V3/W3 mass matrix is block-diagonal (DG) — invert
   cellwise, substitute into the momentum block.
3. Remaining **Helmholtz-type mixed system in (Δu, ΔΠ)**: solve via
   (i) hybridization [K4]: break V2 facet continuity, introduce trace
   multiplier λ on facets, static-condense cell-locally (numpy dense per
   cell), CG+AMG on the SPD-ish trace system; or
   (ii) fallback: PETSc `fieldsplit` Schur (registry entry
   `si_fieldsplit`, App. A) on the monolithic block system.
4. Nonlinear residuals recomputed each Picard iteration; converged when
   ‖residual‖ drops 1e-6 relative or 4 iterations, whichever first (log).

### 3.4 Tracers and moisture (SP6)

Flux-form DG transport of mixing ratios `m ∈ DG1` with dycore mass flux;
upwind facet flux as (3.2e); SSPRK3; TVB minmod slope limiter (Cockburn–Shu)
+ positivity clip with mass fixer (log every activation). Consistency:
transporting `m ≡ 1` with the dycore's (ρ, F) must return `m ≡ 1` to 1e-12
(free-stream preservation).

Kessler microphysics (SP6-T3): column-split after dynamics; variables
(θ, qv, qc, qr); saturation adjustment with Tetens
`es(T) = 610.78 exp(17.27 (T−273.15)/(T−35.86)) Pa` ⚠verify constants against
the MPAS Kessler (`src/core_atmosphere/physics/**/module_mp_kessler*`);
autoconversion/accretion/evaporation/sedimentation constants ported from that
file (single source of truth so FE-vs-FV comparisons share microphysics).

### 3.5 Balance solves

**Hydrostatic balance init** (SP4-T1, SP5): given θ₀ profile and top Π_top,
solve the column mixed system for (w-part of V2 trial, Π): weak form of
`c_p θ₀ ∂Π/∂z = −g` with `dot(u,n)=0` at ground; recover
`ρ₀ = p00 Π^(c_v/R_d) / (R_d θ₀)`. Implement as a vertical-only mixed solve
on the full mesh (horizontal terms dropped); Gusto's
`compressible_hydrostatic_balance` is the design reference. Gate: at-rest
hold (§7 SP4-T1).

### 3.6 Structural identities (used across tests)

1. Complex exactness: `‖div(Π_{V2} (n×grad γ))‖ ≤ 1e-10 ‖γ‖` ∀ γ ∈ Vq basis
   sample (curl-then-div = 0 discretely).
2. `perp(perp(u)) = −u` (tangential).
3. H(div) conformity: normal-trace jump of any V2 function = 0 across
   interior facets (sample-based).
4. Local mass balance (3.1)/(3.2e): cellwise flux residual ≤ 1e-10·scale.
5. Solid-body vorticity: q-diagnosis of `u = Ωs ẑ×x` gives `(2Ωs z/R + f)/D`.

---

## 4. Software architecture

### 4.1 Module layout (grows left→right; create only what your task owns)

```
fecore/
  constants.py          # App. A — physical constants, EarthConfig dataclass
  mesh/                 # existing: icosa.py, scvt_dual.py, _shell_utils.py
    quadsplit.py        # SP2-T2: SCVT-primal & icosa → all-quad surface
    columns.py          # SP5-T1: quad→hex / tri→prism columnar shells, terrain
    slice2d.py          # SP4-T0: structured x–z quad meshes + terrain map
    registry.py         # XC-D: named mesh bundles + checksums
  element/
    zoo.py              # SP2: pair factory: pair(name, mesh) → (V2, V3, Vq)
    bdfm.py             # SP2-T3 (LANE-T only): custom BDFM1 via basix
    thetaspace.py       # SP4-T0/K3: Vθ tensor-product scalar space
  operator/
    poisson.py          # existing
    mixed_poisson.py    # SP2-T1
  model/
    swe_linear.py       # SP2-T5
    swe.py              # SP3-T1
    slice_euler.py      # SP4
    euler3d.py          # SP5
    tracer.py           # SP6-T1
  physics/
    moist.py, kessler.py# SP6-T3
  timestep/
    ssprk.py            # SSPRK3 (Shu-Osher form), stage hooks
    si.py               # SP4-T2/SP5-T2 quasi-Newton driver
  solver/
    options.py          # named PETSc option dicts (App. A.3) — single source
    hybrid.py           # K4 hybridization (optional path)
  init/
    balanced.py         # §3.5 balance solves
  verify/
    mms.py, convergence.py   # existing
    norms.py            # relative L2/L∞ via assembled forms (SP1 pattern)
    conservation.py     # time-series trackers + assert helpers
    suites/
      williamson.py galewsky.py slice_cases.py dcmip.py jw.py supercell.py
    references.py       # SP3-T0: pinned reference data loader (~/Data/MPAS/fecore_refs)
  io/
    output.py           # VTX/XDMF writers + scalar time series (SP2-T6)
    latlon.py           # sample fields to lat-lon grids for plots/comparisons
    mpas_compare.py     # SP3-T6/SP5-T5: load MPAS output, common-grid errors
  tests/                # mirrors: test_<module>.py; gates live here
```

### 4.2 API style rules

- Solvers/models are **functions or small classes taking (mesh, config) and
  returning result objects** — follow `solve_surface`'s shape: explicit
  params, NumPy docstrings, `Returns` section. No global state.
- Every model exposes `step(state, dt) -> state` and `run(case, T, dt,
  callbacks) -> Result`; `Result` carries fields + `diagnostics` time series
  (dict of numpy arrays) + provenance (mesh id, pair, options names).
- PETSc options ONLY via `solver/options.py` registry names.
- MPI-clean from SP3-T1 on: no `.gather` of full fields except in io/latlon
  behind `if comm.size == 1` guards or proper parallel gathers; tests must
  pass under `mpirun -n 1` and `-n 4` where marked.
- Case modules (verify/suites) hold ALL case parameters as module-level
  dataclasses (App. B) — tests and models import from there.

### 4.3 Dependency policy

Additions expected (one per task, pinned when first used): `scipy` (SP2-T4
eigensolves), `dolfinx_mpc` (SP4-T0, K5 — may be dropped if K5 fails),
`matplotlib` (SP2-T6, Agg backend only), `pytest-xdist` (optional, SP3+).
NOT allowed without a plan amendment: Firedrake (that's PC-2/PC-11 territory),
pyamg, custom C extensions.

---

## 5. Design decisions and risk register

Decisions: **DECIDED** (fixed), **DEFAULT** (spike may overturn — fallback
named), **OPEN** (David).

| ID | Decision | Status |
|---|---|---|
| D1 | Compatible mixed FE (FEEC) discretization | DECIDED |
| D2 | Horizontal pair: LANE-T `BDFM1×DG1` vs LANE-Q `RTCF1×DG0` (both 2:1 u:h DOF ratio, the Cotter–Shipton criterion) | DEFAULT: decided by K2 spike; planning bias LANE-Q (all-native elements, unified 2-D/3-D) |
| D3 | Vertical/3-D cells: columnar tensor-product (hex if LANE-Q, prism if LANE-T); tets only as last resort | DEFAULT (K2 + SP5-T1) |
| D4 | Equations: vector-invariant SWE → (u,ρ,θ,Π) compressible Euler, geometry-exact deep | DECIDED |
| D5 | Time: SSPRK3 (SWE) → SI quasi-Newton (slice/3-D); HEVI-IMEX fallback | DECIDED (fallback PC-6) |
| D6 | Geometry: affine/bilinear now; isoparametric degree-2 in XC-A | DECIDED |
| D7 | Toolchain: dolfinx 0.10.0 pinned; version bump only as plan amendment (PC-9) | DECIDED |
| D8 | Python prototype; production decision memo SP7-T4 | DECIDED |
| D9 | Transport at lowest order via upwind DG + (if needed) recovered-space upgrade | DEFAULT (PC-7) |
| D10 | Thin-shell conditioning: columnar DOF ordering + column-blocked PC | DEFAULT (K6) |

### Spike register (K-spikes) — where the plan expects revision

Each spike has a **decision table**; outcomes route to the Plan-Change
register (§9). Spikes are ordinary tasks with extra reporting duty: the
report MUST fill in the decision table and flag §9 items.

| Spike | Inside task | Question | Decision table |
|---|---|---|---|
| **K1** | SP2-T1 | Does dolfinx 0.10 assemble/solve Piola-mapped H(div) mixed forms on immersed manifolds (tdim2/gdim3)? Do interior-facet (dS) integrals work there (needed for DG transport)? | Both work → proceed. H(div) fails → PC-2. dS fails → PC-3. |
| **K2** | SP2-T2 | Quad-split SCVT meshes: mesh quality near pentagons; RTCF1×DG0 mixed Poisson rates; basix RT-on-quad availability. | Quality+rates good → **LANE-Q adopted (PC-1)**: skip SP2-T3, retarget SP3/SP5 pairs. Bad → LANE-T: do SP2-T3 (BDFM1), SP5 needs prism elements (raises R2). |
| **K3** | SP4-T0 | Can Vθ (horiz-DG0 ⊗ vert-CG1) be built as a basix custom element with facet DOFs on structured quads, with vertical-only continuity verified on the actual mesh? Are Vθ DOFs co-located with W2 horizontal-facet DOFs? | Works → proceed. Fails → PC-4 (θ ∈ Q1 CG fallback; SI elimination becomes a lumped approximation — note accuracy consequence). |
| **K4** | SP4-T2 | Hand-rolled hybridization/static condensation of the SI Helmholtz block: correct + faster than fieldsplit at slice sizes? | Yes → default solver. No → `si_fieldsplit` permanent (PC-5), revisit at SP7. |
| **K5** | SP4-T0 | dolfinx_mpc periodic constraints on V2 (H(div)) and Vθ at 0.10 compat? | Yes → periodic slice domains. No → walls+sponges variants of the cases (Straka/Schär unaffected; SK94 uses wide domain + side sponges). |
| **K6** | XC-B-T0 | Physical thin shell (H/R≈3e-3, anisotropy 10²–10³): does column-ordered DOF numbering + column-block / line-smoother PC give mesh-independent-ish CG iterations? | Yes → XC-B + real-atmosphere 3-D viable as planned. No → operator rescaling / vertical Schur ladder (report options; PC-12). |

### Risk register (beyond spikes)

| # | Risk | Mitigation | Fallback |
|---|---|---|---|
| R3 | BDFM1 not in basix (LANE-T only) | port FIAT dual basis (SP2-T3 spec) | BDM1×DG1 (Staniforth et al. 2013) |
| R6 | Reference data irreproducible | SP3-T0 pins everything with checksums + provenance | regenerate with pinned namelists |
| R7 | Python/serial perf wall before SP7 | MPI early (SP3-T1), reduced-radius 3-D configs | pull SP7-T3 forward (PC-8) |
| R8 | SCVT-dual/quad-split cell quality near pentagons inflates FE constants | K2 measures; pentagon-ring error decomposition is a standing diagnostic | local smoothing study; the OTHER lane |
| R9 | `core_sw` Williamson support ≠ memory | SP3-T0 reads `src/core_sw/Registry.xml` to confirm `config_test_case` semantics before promising cases | drive cases via init files instead |
| R10 | Sign/orientation errors in perp/curl/PV on manifolds | structural tests (§3.6) are mandatory FIRST tests of SP2-T5/SP3-T1 | — |

---

## 6. Verification ladder and standing gates

| Tier | Content | Lands in |
|---|---|---|
| A — exact/MMS | mixed-Poisson MMS; sphere gravity-mode spectra; TC2; DCMIP small-planet rates; transport rates; XC-A/XC-B reruns | SP2, SP3-T2, SP5-T3, SP6-T1 |
| B — idealized | Williamson 2/5/6, Galewsky; SK94, Straka, Schär, hydrostatic IGW; JW; BF bubble; supercell | SP3–SP6 |
| C — mesh sensitivity | icosa vs SCVT-dual vs quad-split vs var-res; pentagon-ring decomposition on every sphere case | SP3-T6/T7, SP5-T5 |
| D — scaling | strong/weak MPI; cost breakdowns | SP7 |
| E — showcase | the three `~/Data/MPAS` cases FE-vs-FV + electricity | SP3-T6, SP5-T5, SP6-T4, XC-B |

Standing conservation gates (implemented once in `verify/conservation.py`,
asserted per case):

| Quantity | Statement | Gate [S] |
|---|---|---|
| Mass (∫D or ∫ρ) | exact by (3.2e) | rel. drift ≤ 1e-12 per benchmark window (accumulation-adjusted: ≤ 1e-14/step) |
| Tracer mass; free-stream | flux-form DG | ≤ 1e-12; m≡1 preserved to 1e-12 |
| Energy | not exact (time truncation) | report always; drift ≤ 1% per window [B, CALIBRATE] |
| Enstrophy (SWE) | APVM dissipates | monotone decay after spinup [S] |
| Local flux balance | per-cell residuals | ≤ 1e-10·scale [S] |

---

## 7. Task specifications

Sizes: S ≈ one focused session, M ≈ 1–2, L ≈ 2–3. ⏸ = hard stop for David.
Every task: protocol §0, report per App. D.

---

### SP2 — Compatible spaces on the sphere

#### SP2-T1 (M) — Mixed Poisson on the manifold [SPIKE K1]
**Goal:** stand up H(div)×L² mixed solves on both surface families; gate R1.
**Owns:** `operator/mixed_poisson.py`, `element/zoo.py` (initial),
`tests/test_mixed_poisson.py`.
**Method:** implement (3.1) with `pair="RT1xDG0"`:
`basix.ufl.element("RT", cell, 1)` on the manifold mesh + `("DG", 0)`;
`basix.ufl.mixed_element`. Solve monolithic with nullspace handling (§3.1);
solver registry `mixed_lu` (MUMPS) at these sizes. API:
```python
solve_mixed_poisson_surface(mesh, R, pair, source_fn, exact_fn=None,
                            petsc="mixed_lu") -> (sigma_h, phi_h, errs)
# errs = {"phi_l2_rel": float|None, "flux_residual_max": float}
```
MMS: `mms.surface_source`/`mms.y42_cart` (exact σ = −grad of Y₄² — derive
tangential gradient analytically in `mms.py` as `y42_surface_grad(xyz, R)`;
its correctness test: FD tangential gradient on a lat-lon grid, rtol 1e-3,
pattern of `test_y42_matches_independent_spherical_harmonic`).
**K1 extension:** a second micro-test assembling a pure interior-facet form on
the icosa surface: `assemble_scalar(jump(g)*jump(g)*dS)` for g ∈ DG0
interpolant of z/R — must run and be finite; and an upwind-flux form with
`dot(w, n_f)('+')` for w ∈ RT1 — must assemble. Record any failure verbatim.
**Tests/gates:**
- [S] `test_mixed_assembles_on_manifold` (icosa refine 2; no exception).
- [S] `test_local_conservation`: max cell residual ≤ 1e-10·max|∫f|.
- [S] `test_phi_rate_rt1`: icosa refines {2,3,4,5}, φ slope ≥ 0.9; σ L2 slope
  ≥ 0.9. (RT1×DG0 is first-order.)
- [S] `test_scvt_mixed_runs`: 480+240 km, finite errors, conservation holds.
- [S] `test_ds_forms_on_manifold` (the K1 extension).
**Decision table K1** in report (§5). Fallbacks: PC-2/PC-3 — do NOT attempt
workarounds inside this task; report and stop.

#### SP2-T2 (M) — Quad-split meshes + LANE decision [SPIKE K2] ⏸
**Goal:** build all-quad sphere meshes from both families; run the same mixed
MMS with RTCF1×DG0; measure quality; decide the lane.
**Owns:** `mesh/quadsplit.py`, `tests/test_quadsplit.py`, K2 report.
**Method:**
- `icosa_quad_surface_mesh(refine)`: split each triangle into 3 quads via
  (vertex, two edge midpoints, centroid), project all new points to the
  sphere. Cells: 3·20·4^refine.
- `scvt_quad_surface_mesh(path)`: for each MPAS primal cell (verticesOnCell,
  ordered), quads = (vertex_i midpoint_i, center, midpoint_{i-1})-style fan:
  (cell center, edge-midpoint_i, vertex_{i+1}, edge-midpoint_{i+1}) — one quad
  per cell vertex; pentagon → 5 quads, hexagon → 6. All points projected to
  radius R. **dolfinx quad vertex ordering is tensor-product
  (v00,v10,v01,v11), NOT cyclic** — convert explicitly and test orientation.
- Quality metrics per mesh: min/max interior angle, aspect ratio, skew;
  aggregate + worst-cells-near-pentagons table (pentagon adjacency from
  registry data).
- Repeat SP2-T1's MMS with `pair="RTCF1xDG0"` on both quad families
  (basix "RT" on quadrilateral cell — if the family enum differs, record it).
**Tests/gates:**
- [S] χ=2, all-on-sphere, positive bilinear-cell Jacobians (sample corners),
  conformity (max cells-per-edge = 2).
- [S] `test_quad_mixed_phi_rate`: icosa-quad slope ≥ 0.9; SCVT-quad finite +
  conserving; [B, CALIBRATE] record slopes.
- Quality: no hard gate — REPORT numbers; angles within [30°,150°] expected
  ⚠verify empirically.
**Decision table K2 → LANE choice (PC-1).** Recommendation logic: if
RTCF1×DG0 rates clean AND worst-quality acceptable → LANE-Q (skip SP2-T3;
all later pair references resolve to RTCF1×DG0 and hex 3-D). Else LANE-T.
⏸ David decides with the report in hand.

#### SP2-T3 (M) — BDFM1 custom element [LANE-T ONLY — skipped under LANE-Q]
**Goal:** the Cotter–Shipton triangle pair's velocity space.
**Owns:** `element/bdfm.py`, `tests/test_bdfm.py`.
**Method:** implement BDFM1 (paper naming; = FIAT/Firedrake "BDFM" degree 2 —
**naming trap, document it**) as a basix custom element
(`basix.ufl.custom_element`): space = {v ∈ [P2]²: v·n|_e ∈ P1(e) ∀e}, dim 9.
Port the dual basis from FIAT's `brezzi_douglas_fortin_marini.py` (BSD;
authoritative); map type contravariant Piola, Sobolev H(div), embedded
superdegree 2. The TNT-elements dolfinx demo is the API template.
**Tests/gates [S]:** dim=9; contains [P1]² (interpolate-evaluate exact,
rtol 1e-12); contains curl(λ1λ2λ3); sampled v·n|_e ∈ P1 (fit residual ≤1e-10);
H(div) conformity across a 2-triangle patch; mixed Poisson (3.1) with
BDFM1×DG1 on icosa {2,3,4}: φ slope ≥ 1.9, conservation ≤1e-10.

#### SP2-T4 (M) — Element-pair mode audit on the sphere
**Goal:** verify no spurious branches for the chosen pair; keep an anchor
showing the rejected pair's defect. Env: add `scipy`.
**Owns:** `verify/suites/modes.py`, `tests/test_modes.py`.
**Method:** nonrotating linear SWE eigenproblem on icosa refine 2 (and quad
equivalent): assemble mixed operator L and block mass M; dense
`scipy.linalg.eigh`-compatible reduction (or `eig` on M⁻¹L via factorization —
sizes ≲ 3k DOF). Analytic: gravity multiplets ω² = gH l(l+1)/R², multiplicity
2(2l+1) counting ±ω, plus a zero-frequency (geostrophic/rotational) block of
dimension dim(Vq)−1 (LANE-Q: Q1) ⚠verify count formula empirically at two
refinements before hard-coding.
**Tests/gates:**
- [S] `test_gravity_multiplets_match`: l=1..4 discrete multiplet mean
  frequencies within 5% of analytic (coarse mesh), multiplicities EXACT.
- [S] `test_zero_mode_count`: matches formula.
- [S] `test_rejected_pair_has_spurious_branch` (**xfail-style anchor**):
  RT1×DG0 on triangles exhibits extra/displaced modes between multiplets
  (assert presence — this is the analogue of SP1's lumped anchor). If it
  unexpectedly passes clean, that's Plan feedback, not silent deletion.

#### SP2-T5 (M) — Linear rotating SWE + geostrophic hold
**Goal:** first time-dependent solver; the steady-geostrophic-modes property.
**Owns:** `model/swe_linear.py`, `tests/test_swe_linear.py`.
**Method:** (3.2f), implicit midpoint, monolithic `mixed_lu` per step (sizes
small); Coriolis `f = 2Ω z/R`; balance solve §3.2; ψ = R²Ωs·(z/R) (solid
body) and ψ = Y₂¹-shaped alternative.
**Tests/gates:**
- [S] structural first: perp identities, solid-body vorticity (§3.6).
- [S] `test_balanced_state_holds`: 10 model days, refine {3,4}:
  ‖u(T)−u(0)‖_L2/‖u(0)‖ ≤ 5e-2 at refine 3 AND decreasing with refinement.
- [S] `test_gravity_wave_frequency`: initialize l=2 gravity eigenmode from
  T4 machinery; FFT of ⟨h,Y₂⁰⟩ time series; frequency within 2% (coarse).
- [S] mass drift ≤ 1e-12.

#### SP2-T6 (S) — IO and diagnostics infrastructure
**Owns:** `io/output.py`, `io/latlon.py`, `verify/conservation.py`,
`verify/norms.py` (refactor SP1's error code into it, keeping poisson.py
API), `tests/test_io.py`. Env: add `matplotlib` (Agg).
**Method:** VTX (`dolfinx.io.VTXWriter`) time series; scalar diagnostics CSV;
`latlon.sample(field, nlat, nlon)` via `dolfinx.fem.Function.eval` with
bounding-box trees (serial first; parallel gather later task if needed).
**Gates [S]:** write-read roundtrip; conservation tracker reproduces
assembled integrals to 1e-14; a TC2-style h field plots without exception.

**SP2 HEADLINE ⏸:** machine-exact local mass balance + expected rates for the
mixed operator on SCVT-derived meshes + clean mode audit ⇒ the compatible
machinery is real on MPAS meshes. Discussion before SP3.

---

### SP3 — Nonlinear shallow water + TRiSK head-to-head

#### SP3-T0 (S) — References, mesh registry, core_sw harness
**Owns:** `verify/references.py`, `mesh/registry.py`, `.git/sdd/` provenance
notes; `~/Data/MPAS/fecore_refs/` (data, not committed).
**Method:** (a) mesh registry: named bundles (SP1 SCVT set; quad-split
derivatives; icosa levels) with file checksums. (b) Read
`src/core_sw/Registry.xml` + namelist templates: confirm Williamson
`config_test_case` support (R9) — record exact case list. Build core_sw
(`make gfortran CORE=sw` per BUILD.md; conda env `mpas`, NOT `fecore`).
Run TC2 + TC5 at 480/240/120 km, pin outputs. (c) Fetch/derive TC5 reference
(STSWM spectral reference if obtainable; else 60 km core_sw run labeled
"surrogate reference" — **logged distinction**). (d) Document Galewsky/TC
parameter sources (App. B) against the papers (⚠verify sweep).
**Gates [S]:** registry loads + checksums; core_sw runs complete & are
readable by `io/mpas_compare.py` stub (cell-center h, u fields).

#### SP3-T1 (L) — Nonlinear SWE core
**Owns:** `model/swe.py`, `timestep/ssprk.py`, `tests/test_swe_core.py`.
**Method:** §3.2 complete: q-diagnosis (Vq per lane), F-projection, APVM,
upwind-DG continuity, SSPRK3; per-stage solves from `solver/options.py`
(`mass_v2_cg`: CG+block-Jacobi rtol 1e-10; DG mass: cellwise direct).
MPI-clean (mpirun -n 4 test marker). Δt default `0.3·h_min·R/√(g h0)`.
Diagnostics: mass, energy `∫ ½D|u|²+½g(D+b)² …` (write exact discrete form in
docstring), enstrophy `∫ ½ q² D`, divergence L2.
**Tests/gates:**
- [S] structural suite §3.6 items 1,2,5 for the lane pair.
- [S] `test_tc2_short`: 1 day at refine 3: mass ≤1e-12, energy drift ≤1e-4
  rel [CALIBRATE], solution finite, div noise bounded.
- [S] `test_mpi4_matches_serial`: 6-h TC2, n=1 vs n=4 h-field L2 diff ≤1e-10.

#### SP3-T2 (S) — Williamson TC2 convergence gate
**Owns:** `verify/suites/williamson.py` (TC2 part), `tests/test_tc2.py`.
**Method:** App. B params; day-5 normalized l2(h) per Williamson definition;
icosa (or icosa-quad) refine {3,4,5}.
**Gates:** [S] slope ≥ 1.7; [B CALIBRATE] abs error at finest recorded as
regression band; [S] conservation set.

#### SP3-T3 (M) — TC5 mountain
**Method:** App. B; 15 days; conical mountain in b (L2-projected into V3;
document projection choice). Compare vs SP3-T0 reference on lat-lon grid.
**Gates:** [S] conservation (mass 1e-12; enstrophy monotone post-spinup);
[B] l2(h) at day 15 vs reference within CALIBRATE band; visual wave train
plot archived in report.

#### SP3-T4 (S) — TC6 Rossby–Haurwitz
**Gates:** [S] stable 14 days; [B CALIBRATE] wavenumber-4 pattern speed vs
reference run (loose — document known sensitivity).

#### SP3-T5 (M) — Galewsky barotropic instability
**Method:** App. B (balance integral via dense quadrature, np.trapezoid on
1e-5-rad grid); unperturbed variant FIRST: hold for 5 days (balance gate),
then perturbed 6-day run.
**Gates:** [S] unperturbed hold: vorticity deviation ≤ CALIBRATE, decreasing
with resolution; [B] day-6 vorticity field vs published pattern (archive
plot; quantitative: correlation vs pinned high-res reference ≥ 0.9
CALIBRATE); [S] conservation set.

#### SP3-T6 (L) — Head-to-head vs core_sw + pentagon decomposition ⏸
**Owns:** `io/mpas_compare.py`, `verify/suites/h2h_sw.py`, report with plots.
**Method:** identical `grid.nc` bundles both models (FE on the derived
dual/quad-split of the SAME file — state this framing clearly in the report);
TC2 (analytic truth), TC5 (common reference), Galewsky (high-res reference).
Errors on common lat-lon grid. Pentagon-ring decomposition: BFS rings 0..8
from pentagon cells (graph from `cellsOnCell`); error norms per ring
(SP1/Tier-A.2 methodology). Cost: wall-clock per simulated day both models
(document hardware; Python overhead stated).
**Gates:** [B] FE error ≤ FV error at matched grid for TC2 at 240 km and
finer (expectation from SP1; if it fails, that IS the scientific result —
report, don't hide); [S] FE pentagon-ring profile flat (ring-0:ring-8 error
ratio ≤ 2) while FV shows the known amplification ⚠verify FV behavior against
Tier A.2 diagnostics. **This is the P2-paper centerpiece.** ⏸

#### SP3-T7 (M) — Variable-resolution SCVT
**Method:** a variable-resolution MPAS mesh (registry entry; e.g. X4/X5-type
if available in `~/Data/MPAS`, else generate/request — if none available,
report BLOCKED early, cheap). TC5 + Galewsky on var-res; transition-zone
noise diagnostics (div field, grid-scale energy).
**Gates:** [B CALIBRATE] no transition-zone error spike > 3× surrounding
rings; qualitative comparison vs core_sw noted.

**SP3 HEADLINE ⏸:** TC5 + Galewsky FE-vs-TRiSK on the same grids with
pentagon story quantified.

---

### SP4 — Nonhydrostatic vertical slice (parallel lane; planar quads, no manifold risk)

#### SP4-T0 (M) — Slice infrastructure [SPIKES K3, K5]
**Owns:** `mesh/slice2d.py`, `element/thetaspace.py`, `tests/test_slice_infra.py`.
Env: add `dolfinx_mpc` (K5; drop again if NO-GO).
**Method:**
- `slice_mesh(Lx, Lz, nx, nz, terrain=None)`: structured quads,
  tensor-ordered vertices; terrain map: Gal-Chen
  `z(x,ζ) = h(x) + ζ(Lz−h(x))/Lz` applied to node coords post-build.
- **K3**: Vθ custom element on quad: value size 1, span{1, ẑ} (vertical P1,
  horizontal P0), DOFs = point evals at horizontal-facet midpoints (entity
  dim 1, the two horizontal edges of the reference quad). Verify on the real
  mesh: (a) vertical continuity: interpolate `sin(z)`-like function, jump
  across horizontal facets ≤1e-12, jump across vertical facets ≠ 0 allowed;
  (b) DOF co-location with W2 (RTCF1) horizontal-facet DOFs (coordinates
  match ≤1e-12).
- **K5**: periodic-in-x MPC on V2 and Vθ; assemble+solve a Poisson-type toy.
**Gates [S]:** mesh cell count/orientation; terrain mesh min Jacobian > 0
for Schär params; K3 (a),(b); K5 recorded in decision table (either outcome
acceptable — routes case setup).

#### SP4-T1 (M) — Slice space stack + hydrostatic balance
**Owns:** `init/balanced.py`, `model/slice_euler.py` (state container +
residual assembly only), `tests/test_slice_balance.py`.
**Method:** §3.5 column balance solve; isothermal and N=const θ-profiles;
state container (u∈RTCF1, ρ∈DG0, θ∈Vθ).
**Gates [S]:** `test_at_rest_stays_at_rest`: balanced state, 100 steps
explicit tiny-dt (or 10 SI steps once T2 lands — initial version explicit):
max|u| ≤ 1e-8 m/s; residual of (3.3e)+gravity h-converges (two meshes).

#### SP4-T2 (L) — Semi-implicit quasi-Newton stepper [SPIKE K4]
**Owns:** `timestep/si.py`, `solver/hybrid.py` (optional), `solver/options.py`
entries `si_fieldsplit`, `helmholtz_amg`; `tests/test_si.py`.
**Method:** §3.3 SI algorithm; θ/ρ elimination as specified; Helmholtz via
fieldsplit FIRST (make it correct), then K4 hybridization attempt (correct +
timing comparison at 300×100 mesh). α=0.55 default off-centering ⚠tunable.
**Gates:**
- [S] SK94-linear (small-amplitude): stable at Δt = 12 s where explicit
  acoustic CFL would need Δt ≲ 1 s (i.e., SI is doing its job; exact factor
  CALIBRATE);
- [S] energy drift bounded (report), balanced state still held with SI;
- K4 decision table (correctness: hybridized vs fieldsplit solutions agree
  to 1e-8; timing recorded).

#### SP4-T3 (M) — SK94 / Baldauf–Brdar gravity wave
App. B params. **Gates:** [B] θ′ field vs Baldauf–Brdar semi-analytic
solution: L2 ≤ CALIBRATE and self-convergence slope ≥ 1.5; [S] symmetry
(‖θ′(x)−θ′(−x+2xc)‖ ≤ 1e-10 with symmetric setup); conservation set.

#### SP4-T4 (M) — Straka density current
App. B; explicit viscosity terms `μ ∇²` weak-formed (interior penalty for the
DG-direction of θ — document the form). **Gates:** [B] front location at
t=900 s within 300 m of pinned reference table value at Δx=50 m ⚠verify; θ′_min
within band; [S] self-convergence of front position, no-viscosity run
excluded (it must NOT be run as a gate — under-resolved).

#### SP4-T5 (M) — Schär mountain wave
App. B; terrain-following mesh; Rayleigh sponge (γ(z) ramp, top 8 km ⚠verify
depth) on w implemented as implicit relaxation.
**Gates:** [B] normalized vertical momentum flux within 10% of linear-theory
value across z ∈ [3,8] km CALIBRATE; [S] no 2Δx noise above terrain
(spectral diagnostic: energy in top-octave wavenumbers ≤ CALIBRATE fraction);
comparison note vs MPAS-A's known metric-term sensitivities (qualitative).

#### SP4-T6 (S) — Hydrostatic-limit IGW
Wide-domain (6000 km) variant; gate [B]: matches hydrostatic analytic phase
speeds within 2%; documents that the nonhydrostatic core has the correct
hydrostatic limit.

**SP4 HEADLINE ⏸:** all four canonical slice cases within bands, SI stable at
acoustic-CFL-free Δt.

---

### SP5 — 3-D spherical shell dry core

#### SP5-T1 (M) — Columnar shell meshes + 3-D space stack
**Owns:** `mesh/columns.py`, `tests/test_columns.py`.
**Method (LANE-Q):** extrude quad-split surface meshes radially (reuse
`_shell_utils` layer logic, hexes need NO diagonal splitting — simpler than
SP1's tets): `quad_shell_mesh(surface, n_layers, H, stretch=None,
terrain=None)`; tensor-ordered hex connectivity; optional Gal-Chen terrain.
3-D stack: W2=RT(hex,1), W3=DG0, W1=N1curl(hex,1), Vθ = 3-D analogue
(DOFs on horizontal quad facets; extend `element/thetaspace.py`).
(LANE-T: prisms instead — this task becomes the prism-element spike, sized L,
with the K2 report's guidance; expect plan amendment.)
**Gates [S]:** conformity (max cells-per-facet 2), positive Jacobians
(sampled), volume ratio vs analytic shell ∈ [0.95, 1.01] for quad-split
(bilinear faces recover volume better than flat tris — CALIBRATE), Vθ
vertical-continuity + co-location tests on the shell, at-rest balanced
atmosphere holds (reuse SP4-T1 machinery on one column of cells first, then
globally).

#### SP5-T2 (L) — 3-D dry SI dycore
**Owns:** `model/euler3d.py`, extensions to `timestep/si.py`,
`tests/test_euler3d.py`.
**Method:** §3.3 in 3-D: Coriolis 2Ω×u; vector-invariant with ω ∈ W1;
SI with columnar Δθ elimination; Helmholtz per SP4-T2 outcome; column-ordered
DOF permutation for preconditioning (D10; implement as PETSc ordering or
assembled-block structure — document).
**Gates [S]:** balanced solid-body zonal flow on the shell (small-planet
X=20 config for cheapness) holds 10 days: max|u−u0| ≤ CALIBRATE and
h-decreasing; conservation set; MPI n=4 reproducibility ≤ 1e-10.

#### SP5-T3 (M) — DCMIP small-planet convergence [the 3-D Tier-A gate]
**Method:** DCMIP-2012 case 3-1 (nonhydrostatic gravity wave, X=125) App. B;
coupled refinement (horizontal refine ⇒ n_layers doubles — SP1 Task-8
methodology); self-convergence against finest (Richardson-style).
**Gates:** [S] θ′ L2 self-convergence slope ≥ 1.7; [B] θ′ extrema vs published
DCMIP model range CALIBRATE.

#### SP5-T4 (L) — Jablonowski–Williamson baroclinic wave
**Method:** App. B init (η-coordinate profiles → z via hypsometric iteration;
algorithm given in App. B); unperturbed steady-state hold FIRST (4-day ps
drift gate), then perturbed 10-day run at ~120 km + 30 layers, MPI.
**Gates:** [S] unperturbed: ps drift ≤ CALIBRATE hPa at day 4, h-decreasing;
[B] perturbed: day-9 min ps and wave-train position within the published
multi-model envelope ⚠verify (JW06 tables + MPAS-A reference run from
`~/Data/MPAS/jw_baroclinic_wave` — pin via SP3-T0-style provenance);
conservation set.

#### SP5-T5 (M) — Head-to-head vs MPAS-A (dry JW) ⏸
As SP3-T6 but vs `atmosphere_model` (JW config pinned): error-vs-cost,
KE spectra (diagnostic in `io/latlon.py` + spherical-harmonic transform via
scipy — document truncation), pentagon-ring diagnostics on ps/vorticity.
**Gates:** [B] comparable wave evolution (report); [S] no pentagon signature
in FE fields (ring ratio ≤ 2). ⏸

**SP5 HEADLINE ⏸:** JW on MPAS's grid, FE vs FV, plus a genuine 3-D
convergence rate. This is "a working dry FE dycore".

---

### SP6 — Transport, moisture, supercell

#### SP6-T1 (M) — DG tracer transport + Nair–Lauritzen suite
**Owns:** `model/tracer.py`, `verify/suites/nl_transport.py`, tests.
**Method:** §3.4; DG1 on the sphere (K1's dS result governs; LANE-Q: DG1 on
quads); prescribed velocity fields (solid body; NL deformational cases 1&4 —
App. B); SSPRK3; TVB minmod limiter + clip-and-rescale positivity (log).
**Gates:** [S] Gaussian-hills variant unlimited rate ≥ 1.9; tracer mass
≤1e-12; limited: min ≥ −1e-14, maxima non-amplifying; [B CALIBRATE] cosine
bells standard mixing diagnostics recorded.

#### SP6-T2 (S) — Consistent coupling / free-stream
Couple tracer to SWE then 3-D mass flux; **gate [S]:** m≡1 preserved to 1e-12
over a 5-day TC5 flow and a 1-day 3-D run.

#### SP6-T3 (M) — Moist thermodynamics + Kessler + Bryan–Fritsch bubble
**Owns:** `physics/moist.py`, `physics/kessler.py`, tests.
**Method:** §3.4; Kessler constants ported from the MPAS source file (single
source of truth — record file+hash); θ→θ_ρ (density temperature) coupling in
buoyancy/EOS documented; column split-step after dynamics.
**Gates:** [S] Kessler column unit tests vs hand-computed saturation cases;
water mass (vapor+cloud+rain, incl. surface flux ledger) ≤ 1e-12; [B]
Bryan–Fritsch moist bubble: w_max and cloud-top height at t=1000 s within
published band ⚠verify (BF2002).

#### SP6-T4 (L) — Reduced-radius supercell vs MPAS-A ⏸
**Method:** KSP15 configuration — DO NOT re-derive the sounding: read
`~/Data/MPAS/supercell` namelists/init and initialize fecore by interpolating
the MPAS init fields (u, θ, ρ, qv) into FE spaces (`io/mpas_compare.py`
gains an init-import mode; document interpolation choices). X=120 sphere,
Kessler both models.
**Gates:** [B] storm splits by 90 min; updraft max ∈ published/MPAS-run band
CALIBRATE; rainwater structure compared side-by-side (archived plots);
[S] conservation set (water ledger).

**SP6 HEADLINE ⏸ = Definition-of-Done moist milestone.**

---

### SP7 — Performance, scaling, production decision

- **SP7-T1 (M)** MPI strong/weak scaling: TC5 (SWE), slice GW, 3-D JW at two
  sizes; assembly/solve/other breakdown via PETSc log + timers. Gate:
  report exists with scaling curves to the hardware limit (document).
- **SP7-T2 (M)** Cost per DOF-step vs MPAS-A on JW at matched grid: honest
  table incl. Python overhead; error-vs-cost plot combining SP5-T5.
- **SP7-T3 (S)** Assessment memo: matrix-free/sum-factorization, C++ forms,
  GPU via PETSc, dolfinx C++ driver — with microbenchmarks where cheap.
- **SP7-T4 (M)** Production-path decision memo (D8 options, criteria,
  estimates). ⏸ David decides.

---

### XC — Cross-cutting threads

- **XC-A (M) Isoparametric geometry:** degree-2 coordinate element on
  surface/shell meshes (midpoints snapped to sphere); rerun SP1 matrix + T2-1
  mixed MMS. **Gates [S]:** P2/BDM-class slopes ≥ 2.7 (expect ≈3.0 —
  validates the SP1 cap diagnosis); SP1 tests unchanged at degree-1 (no
  regression).
- **XC-B (L) Electricity payoff [SPIKE K6]:** physical-geometry shell solver:
  `XC-B-T0` (S) K6 preconditioning spike: anisotropic hex columns at
  R=6.371e6, H=6e4?, layers stretched (Δz 100 m→5 km ⚠match the Poisson
  design docs' vertical grid), conductivity σ(z)=σ0 exp(z/h_σ) ⚠take profile
  from `docs/design/poisson/governing-equations.md`; column-ordered
  block-Jacobi/line PC vs hypre baseline; gate: iterations ≤ 200 and
  mesh-stable CALIBRATE. `XC-B-T1` (M): full solve div(σ grad φ)=−S with
  ground Dirichlet/top Dirichlet-or-Neumann per the design doc; MMS at
  physical scale (rescaled Y₄²·R(ζ) — extend `mms.py` with σ(z) variant);
  Tier A.2-equivalent gate ≥1.9 **at physical thin-shell geometry**; then
  Tier-B idealized charge cases from the design doc. Positioning: the
  Phase-2 operator candidate for the Poisson paper thread.
- **XC-C (S, recurring) RTD design narrative:** `docs/design/fecore/`
  mirroring poisson/shdom structure (index, governing-equations,
  discretization, verification-plan pages); update at every SP boundary;
  Zenodo DOI at paper milestones.
- **XC-D (S) Mesh registry hardening:** fold into SP3-T0 if convenient;
  otherwise named-bundle registry + checksums + var-res entries.

---

## 8. Sequencing

```
 SP2-T1 (K1) ──► SP2-T2 (K2) ──⏸ LANE ──► [SP2-T3 if LANE-T] ─► SP2-T4 ─► SP2-T5 ─► SP2-T6 ─⏸
                                   │
 SP4-T0 (K3,K5) ─► SP4-T1 ─► SP4-T2 (K4) ─► SP4-T3..T6 ─⏸        (parallel lane, start anytime)
                                   │
 SP3-T0 ─► SP3-T1 ─► T2..T5 ─► SP3-T6 ⏸ ─► T7 ─⏸
                                   │
 (needs SP3 + SP4 + lane) SP5-T1 ─► SP5-T2 ─► SP5-T3 ─► SP5-T4 ─► SP5-T5 ⏸
                                   │
 SP6-T1..T2 (after SP3; T1 also gated by K1-dS) ─► SP6-T3 ─► SP6-T4 ⏸ ─► SP7 ─⏸
 XC-B-T0 (K6) anytime after SP2-T2;  XC-A after SP2;  XC-C recurring;  XC-D with SP3-T0
```

**Recommended execution order (first ten):** SP2-T1, SP2-T2 ⏸, SP4-T0,
SP2-T4 (or SP2-T3 if LANE-T), SP2-T5, SP2-T6, SP3-T0, SP3-T1, SP4-T1,
XC-B-T0. Two agents can run the SP2/SP3 and SP4 lanes concurrently once the
lane decision is made (disjoint file ownership by design).

Scale estimate: ~38 tasks ≈ 5× SP1 — a multi-month programme; sub-project
boundaries are the re-planning points and this file is amended there
(version bump + change note).

---

## 9. Possible plan changes (anticipated revision register)

The spikes exist because these are the known unknowns. Each entry: trigger →
change → blast radius. Implementers flag triggers in reports; only David
amends the plan.

- **PC-1 — LANE-Q adoption** (trigger: K2 GO on quad-split quality+rates;
  planning bias says likely). Change: SP2-T3 (BDFM1) deleted; all pair
  references → RTCF1×DG0; SP3 runs on quad-split meshes (framing: "same SCVT
  generating set"); SP5 hexes as specced. Blast radius: small — this plan is
  already written LANE-Q-first. Paper framing gains "LFRic-family elements on
  MPAS meshes" angle; keep one triangle-mesh mixed-Poisson test alive as a
  cross-check.
- **PC-2 — Manifold H(div) NO-GO** (trigger: K1 assembly/solve failure).
  Options, in order: (a) dolfinx version bump if upstream fix exists (PC-9
  procedure); (b) hand-rolled RT0/RTCF1 surface assembly in numpy/scipy
  against our own mesh classes (bounded work: lowest-order Piola on flat
  cells is classical; ~2 tasks, spec to be written); (c) Firedrake for
  sphere SPs only (env split; D7 amendment). Slice lane (SP4) proceeds
  regardless — planar is unaffected.
- **PC-3 — Manifold interior-facet (dS) NO-GO** (trigger: K1 extension
  fails). Change: sphere transport switches to CG+SUPG or to a
  vertex/cell FV fallback for D and tracers (SP3-T1, SP6-T1 rescoped);
  slice/3-D unaffected (planar/volume dS fine).
- **PC-4 — Vθ custom element fails** (trigger: K3). Change: θ ∈ Q1 (full CG)
  Lorenz-flavored compromise; SI Δθ elimination becomes approximate (lumped
  vertical mass); note expected extra dispersion in reports; revisit at SP5.
- **PC-5 — Hybridization abandoned** (trigger: K4 slower/incorrect).
  fieldsplit-Schur becomes the permanent SI solver; SP7-T3 revisits; no
  science impact, cost impact measured.
- **PC-6 — SI quasi-Newton inadequate on the shell** (trigger: SP5-T2
  instability/iteration blowup at operational Δt). Change: add SP5-T2b
  HEVI-IMEX (ARK2-class, horizontally-explicit vertically-implicit; column
  tridiagonal solves — the columnar DOF layout from D10 makes this cheap).
- **PC-7 — Lowest-order transport too diffusive** (trigger: SP5-T4 JW wave
  train visibly damped vs MPAS at matched grid, or SP6-T4 storm structure
  lost). Change options: recovered-space transport (Bendall et al. —
  DG0→CG1 recovery, transport in DG1, project back; add task SP6-T1b), or
  order-up the whole hierarchy (RTCF2×DG1 + Vθ degree 2 — bigger, touches
  SP3/SP5 rates for the better).
- **PC-8 — Performance wall pre-SP7** (trigger: any SP5 task can't run in
  wall-clock ≤ a workday at specced resolution). Change: pull SP7-T3
  forward; drop resolutions (small-planet more); consider C++ assembly of
  the hot forms.
- **PC-9 — dolfinx version bump** (trigger: needed bugfix). Procedure, not
  really a change: new env pin + FULL suite re-run (SP1 tests included) + all
  CALIBRATE bands re-validated before any new work lands on the new version.
- **PC-10 — Priority inversion: electricity first** (trigger: David's Q4
  answer). Change: XC-B promoted immediately after SP2-T2; SP3+ shifts right;
  P1 paper fast-tracked on SP1+SP2+XC-B evidence.
- **PC-11 — Both TP-element routes dead AND quad-split bad** (trigger: K2
  fail + SP5-T1 LANE-T fail). Last resort: tets + fully-coupled 3-D solves
  (SP1 shell meshes reused); physics stays idealized-only; SP5 gates relaxed
  to slice-parity; production story weakens — flag loudly.
- **PC-12 — Thin-shell preconditioning fails** (trigger: K6). XC-B proceeds
  at reduced aspect (H exaggerated ×10) with a documented conditioning study
  as the deliverable; real-geometry solve becomes an SP7-adjacent research
  item.
- **PC-13 — Gate-threshold amendments.** Any [S] threshold change is a plan
  amendment recorded here with rationale (e.g., a rate gate lowered from 1.9
  to 1.7 for a case where theory says the limiter caps it). Never silent.

---

## 10. Papers and deliverables (unchanged from v0.1, task-mapped)

| # | Working title | Evidence from | Venue class |
|---|---|---|---|
| P1 | Compatible FEM repairs the SCVT convergence defect | SP1 + SP2 + XC-B | GMD / JAMES |
| P2 | Compatible-FE shallow water on MPAS meshes: head-to-head with TRiSK | SP3 (+T6/T7 Tier C) | QJRMS / MWR / JAMES |
| P3 | A nonhydrostatic compatible-FE prototype on MPAS meshes | SP4 + SP5 + SP6 | JAMES / GMD |

Each milestone: XC-C narrative update, tagged release, Zenodo DOI.

---

## Appendix A — Constants, configs, solver registry

**A.1 `fecore/constants.py`** (create in SP3-T1; values ⚠verify against
`src/framework`/MPAS constants module so FE-vs-FV comparisons share physics):

```python
@dataclass(frozen=True)
class EarthConfig:
    g: float = 9.80616          # m s^-2  (MPAS/Williamson value)
    omega: float = 7.29212e-5   # s^-1
    R: float = 6.37122e6        # m
    X: float = 1.0              # small-planet reduction factor (R/X, omega*X)
R_d = 287.0; c_p = 1004.5      # ⚠verify vs MPAS constants (c_v = c_p - R_d)
p00 = 1.0e5
```
Small-planet rule: radius R/X, rotation Ω·X (DCMIP convention) ⚠verify.

**A.2 Mesh/DOF quick reference** (icosa: cells 20·4^r, edges 30·4^r,
verts 10·4^r+2; quad-split ×3 cells; SCVT 480 km: 2562 cells / 5120 dual
tris). Δx_icosa ≈ R·h with h = 1/2^r (SP1 proxy convention).

**A.3 `solver/options.py` initial registry:**

```python
OPTIONS = {
  "mixed_lu":      {"ksp_type":"preonly","pc_type":"lu","pc_factor_mat_solver_type":"mumps"},
  "mass_v2_cg":    {"ksp_type":"cg","pc_type":"bjacobi","ksp_rtol":1e-10},
  "poisson_hypre": {"ksp_type":"cg","pc_type":"hypre"},          # SP1 heritage
  "si_fieldsplit": {"ksp_type":"gmres","pc_type":"fieldsplit",
                    "pc_fieldsplit_type":"schur","pc_fieldsplit_schur_fact_type":"full",
                    "fieldsplit_0_ksp_type":"preonly","fieldsplit_0_pc_type":"bjacobi",
                    "fieldsplit_1_ksp_type":"cg","fieldsplit_1_pc_type":"hypre"},
  "helmholtz_amg": {"ksp_type":"cg","pc_type":"hypre","ksp_rtol":1e-8},
}
```
(Names are the API; tune values freely within a task, record changes.)

## Appendix B — Case parameter tables (⚠verify each against the cited paper at first use; record in the case module docstring)

- **Williamson TC2:** u0 = 2πR/(12 days); g·h0 = 2.94e4 m² s⁻²;
  h = h0 − (R Ω u0 + u0²/2) sin²φ / g; α=0. Error norm: Williamson normalized
  l2 of h at day 5. [W92 §3.2]
- **TC5:** u0=20 m s⁻¹, h0=5960 m; mountain b = b0(1−r/R_m), b0=2000 m,
  R_m=π/9, center (λc,φc)=(−π/2, π/6), r=min(R_m, great-circle arg); 15 days.
  [W92 §3.5]
- **TC6:** R-H wavenumber 4, ω=K=7.848e-6 s⁻¹, h0=8000 m; 14 days. [W92 §3.6]
- **Galewsky:** u_max=80 m s⁻¹; φ0=π/7, φ1=π/2−φ0, e_n=exp(−4/(φ1−φ0)²),
  u(φ)=u_max/e_n · exp(1/((φ−φ0)(φ−φ1))) on (φ0,φ1) else 0; h by balance
  integral (numeric); mean depth 10 158 m ⚠verify; perturbation
  h′=120 m · cosφ · exp(−(λ/α)²) · exp(−((φ2−φ)/β)²), α=1/3, β=1/15, φ2=π/4;
  6 days. [Galewsky et al. 2004]
- **SK94 gravity wave:** domain 300 km × 10 km; θ0=300 K, N=0.01 s⁻¹
  (θ̄=θ0 e^{N²z/g}); Δθ=1e-2 K, θ′=Δθ sin(πz/H)/(1+((x−xc)/a)²), a=5 km,
  xc=100 km; U=20 m s⁻¹; T=3000 s; periodic-x (or wide+sponge per K5).
  Reference: Baldauf & Brdar (2013) semi-analytic. [SK94]
- **Straka:** domain 51.2 × 6.4 km (or half+symmetry); isentropic θ0=300 K;
  bubble ΔT=−15 K·(cos(πL)+1)/2, L=√(((x−xc)/4 km)²+((z−zc)/2 km)²) ≤ 1,
  zc=3 km; ν=μ_θ=75 m² s⁻¹; t=900 s; reference front ≈15.5 km ⚠verify table.
  [Straka et al. 1993]
- **Schär:** h(x)=h0 e^{−(x/a)²} cos²(πx/λ), h0=250 m, a=5 km, λ=4 km;
  U=10 m s⁻¹, N=0.01 s⁻¹; sponge above ⚠verify domain/sponge from paper;
  steady-state momentum-flux profile vs linear theory. [Schär et al. 2002]
- **DCMIP-2012 3-1:** X=125, N=0.01, u_eq=20 m s⁻¹, θ′ Gaussian d=5 km at
  equator, amplitude 1 K ⚠verify all against the DCMIP document; T=3600 s.
- **JW baroclinic wave:** full JW06 η-profiles; z from hypsometric iteration:
  per column, evaluate T(φ,η),u on η-grid, integrate z(η)=z(1)+∫R_d T dlnp/g,
  invert η(z) by monotone interpolation to mesh levels; perturbation zonal
  wind Gaussian (λ,φ)=(π/9, 2π/9), radius R/10, amp 1 m s⁻¹ ⚠verify. Gates
  vs JW06 tables + pinned MPAS-A run. [JW06]
- **Bryan–Fritsch bubble:** saturated θ_e=320 K background; 2 K perturbation
  ⚠verify all from BF2002.
- **Supercell:** DO NOT re-derive — import init from `~/Data/MPAS/supercell`
  (KSP15 config). X=120.

## Appendix C — dolfinx 0.10 API gotchas (extend as discovered; SP1-verified items marked ✓)

1. ✓ `create_mesh(comm, cells, element, coords)` — element BEFORE coords.
2. ✓ `Function.interpolate` passes x with shape (3, N) — transpose for
   `f(xyz, ...)` conventions.
3. ✓ `LinearProblem(..., petsc_options_prefix=...)` — keyword REQUIRED.
4. ✓ `locate_entities_boundary` marker receives (3, n) array.
5. ✓ `assemble_scalar(fem.form(expr))` — wrap in `form` first.
6. ✓ After manual PETSc solves: `uh.x.scatter_forward()`; destroy PETSc
   objects (SP1 `_assemble_and_solve_lumped` pattern).
7. Quad/hex cells use TENSOR-PRODUCT vertex ordering (v00,v10,v01,v11 /
   +z-layer), not cyclic — converters in mesh builders + orientation tests.
8. Mixed spaces: `basix.ufl.mixed_element([...])`; subfunctions
   `w.sub(i).collapse()`; BCs on subspaces need `(V.sub(i), collapsed_V)`
   pairs.
9. Nullspaces: build `PETSc.NullSpace`, `A.setNullSpace(ns)` (+
   `setNearNullSpace` for AMG); orthogonalize RHS.
10. Manifold forms: use `SpatialCoordinate`-based analytic normals
    (`x/R`); do NOT trust `CellNormal` orientation until tested (K1 records
    behavior). `FacetNormal` on manifolds = in-surface conormal — verify in
    K1's dS probe.
11. `dolfinx.fem.petsc.LinearProblem.solve()` returns the Function; KSP
    convergence must be checked (`prob.solver.getConvergedReason()` — SP1
    checked `ksp.is_converged` on manual paths).
12. Custom elements: `basix.ufl.custom_element(...)` wrapping
    `basix.create_custom_element` — follow the TNT demo; conformity comes
    from entity-attached interpolation DOFs, so co-location tests (K3) are
    behavioral, not assumed.

## Appendix D — Task report template (`.git/sdd/task-SPx-Ty-report.md`)

```markdown
# Task SPx-Ty — <title>
Date / model / env hash (conda list --explicit | sha256sum)
## Goal (one paragraph, from the plan)
## What was done (files created/changed; brief diff narrative)
## API surprises / gotchas (→ also added to Appendix C? y/n)
## Results
   - test counts before/after; gate table: name | threshold | measured | class
   - CALIBRATE values established (if any) — exact numbers
## Verdict: GO / NO-GO (+xfail anchor name) / BLOCKED (+what unblocks)
## Spike decision table (spike tasks only)
## Plan feedback (anything §9 should absorb; threshold amendment requests)
```

## References (implementers: cite-check ⚠verify items against these)

Cotter & Shipton 2012 (JCP, arXiv:1103.2440) — 2:1 ratio, BDFM1.
Staniforth, Melvin & Cotter 2013 (QJRMS) — BDM1×DG1 analysis.
Natale, Shipton & Cotter 2016 — compatible-FE dycore review (W-spaces).
McRae & Cotter 2014 — energy–enstrophy SWE. Shipton, Gibson & Cotter 2018
(JCP, arXiv:1707.00855) — sphere SWE implementation. Gibson et al. 2020
(GMD) — hybridization. Melvin et al. 2019 (QJRMS), 2023 (GMD 16:1265) —
LFRic slice & SWE. Bendall et al. — recovered spaces / Gusto. Rognes, Ham,
Cotter & McRae 2013 (GMD 6:2099) — manifold FEM. Thuburn 2009 / Ringler
et al. 2010 — TRiSK (`core_sw`). Danilov 2010; Peixoto 2016 — triangular
C-grid modes; pentagon accuracy. Williamson et al. 1992; Galewsky et al.
2004. Skamarock & Klemp 1994; Baldauf & Brdar 2013; Straka et al. 1993;
Schär et al. 2002. Jablonowski & Williamson 2006. Ullrich et al. DCMIP
2012/2016. Klemp, Skamarock & Park 2015 — reduced-radius cases (supercell).
Bryan & Fritsch 2002. Dziuk 1988; Bernardi 1989; Demlow 2009 — surface-FEM
geometry error. Baratta et al. 2023 — DOLFINx. Scroggs et al. — basix custom
elements. Dokken — dolfinx_mpc.

---

*Plan of record, implementation edition. Amendments: version bump + §9 note.
Authored by Claude Fable 5 building on SP1 (Claude Opus 4.8 / Sonnet 4.6).
Filename honors the planning model; the gates honor nobody — they just fail
honestly.*
