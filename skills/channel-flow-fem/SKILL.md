---
name: channel-flow-fem
description: "Steady and oscillatory incompressible flow in 2-D microfluidic channels with open-source finite elements (scikit-fem + gmsh), with closed-form benchmarks built in. Use for wall shear stress, cavity vortices, ridged or stepped channel geometry, pulsatile or oscillating flow, Womersley number, TAWSS and OSI, or replacing an ANSYS Fluent / COMSOL laminar-flow model with free software."
---

# Channel flow with finite elements, benchmarked first

Taylor–Hood P2/P1 elements on gmsh meshes. Steady Navier–Stokes by Picard iteration; transient
Navier–Stokes by BDF2 with an implicit-Euler start. Wall shear stress is extracted as the
tangential wall traction. Everything is plain Python; the largest case tested so far
(about 200,000 unknowns) peaks near 2.2 GB of memory and a three-level mesh study with a
four-point flow-rate sweep runs in about four minutes on 16 cores.

## Rule zero: run the benchmarks before believing anything

```bash
pip install numpy scipy matplotlib scikit-fem gmsh
python3 scripts/benchmarks.py all
```

Expected on a correct installation:

| Benchmark | Expected | Graded |
|---|---|---|
| Steady parallel-plate wall shear, 500 µm channel, 6 mL/min | 0.336000 Pa; relative error ~1e-15 | yes, < 1e-8 |
| Womersley flow, 1 Hz (Wo ≈ 0.75), 25/50/100 steps per cycle | errors ≈ 9.1e-5, 1.8e-5, 4.1e-6; ratio ≈ 4–5 per halving | yes |
| Womersley flow, 100 Hz (Wo ≈ 7.5) | errors stall near 3e-3 | no — spatial floor, see below |

If the steady benchmark is not at machine precision, something is wrong with the installation
or a change you made — stop there. The exact Poiseuille profile lies inside the P2 space, so
any visible error is a bug, not under-resolution.

When adapting to a new geometry or fluid, first change the benchmark to the new channel height
and fluid properties and re-run it. A benchmark in someone else's units verifies nothing.

## Scripts

All in `scripts/`, SI units throughout. They import each other, so run them from `scripts/` or
put that directory on `PYTHONPATH`.

| File | Provides |
|---|---|
| `make_mesh.py` | `build(path, n_ridge, lc_bulk, lc_wall, L_IN, L_OUT, H)` — 2-D channel with `n_ridge` rectangular floor ridges (100 µm wide and high, 400 µm gaps, i.e. 500 µm pitch, by default); `n_ridge=0` gives a plain channel. Boundary tags: `inlet`, `outlet`, `bottom`, `top`. Returns a dict of mesh statistics. |
| `solve_ns.py` | `solve_case(mesh_path, Q_mL_min, H)` — steady solve with a developed parabolic inlet; `wall_shear(res, "bottom")` → `(x, z, tau)`, sign = flow direction. Constants `RHO`, `MU`, `W_CHAN`. |
| `solve_unsteady.py` | `solve_transient(mesh_path, inlet_fn, dt, n_steps, H, u0_fn, linear, on_step)`; `womersley(z, t, H, dpdx_amp, omega)` closed-form profile; `womersley_number`, `stokes_layer`; `ShearAccumulator` for cycle-averaged TAWSS and OSI. |
| `benchmarks.py` | The checks above, as a command. Non-zero exit on failure. |

Fluid properties default to values close to cell-culture medium
(ρ = 998 kg m⁻³, μ = 7 × 10⁻⁴ Pa s); change `RHO`, `MU` in `solve_ns.py` for anything else.
Flow rate is converted to mean velocity through the spanwise width `W_CHAN` (5 mm), so a 2-D
result stands for the mid-plane of a wide, shallow channel.

## Resolution requirements found the hard way

- **Mesh convergence on the quantity you report.** Refine `lc_wall` and `lc_bulk` together
  (10, 7 and 5 µm wall spacing is a reasonable first ladder for 500 µm channels) until that quantity — peak position, peak
  shear, recirculation length — stops moving. A field plot that looks converged is not evidence.
- **Oscillatory flow needs the Stokes layer resolved.** Thickness δ = √(2ν/ω): 47 µm at 100 Hz
  in water-like medium. At 8 µm wall spacing the error floor is 3.5e-3; at 3 µm (≈ 15 elements
  per δ) it is 5e-4, after which time stepping limits it again. Check δ before choosing a mesh.
- **Validated range.** Steady solves converged in 11–14 Picard iterations for Re ≈ 40–170
  (hydraulic-diameter based). Higher Re is untested here; expect to need under-relaxation
  (`relax` in `solve_case`) and eventually a turbulence or 3-D question this code does not answer.

## Implementation notes

- scikit-fem mixed assembly needs the same `intorder` on the velocity and pressure bases, or it
  fails with a quadrature mismatch.
- `FacetBasis` normals point out of the fluid, so the traction is negated before taking its
  tangential part; `wall_shear` returns shear as positive when it acts downstream.
- `gmsh` boolean cut fails on an empty tool list; `make_mesh.build` guards the ridge-free case.
- Laminar is not the same as Stokes. A low-Re channel with cavities still needs the convective
  term: Stokes flow produces no cavity vortices at all.
