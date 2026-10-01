# fiducial-lab

Verified computational methods for biomedical-engineering research — for use with Claude, or
as plain Python without it.

A *fiducial* is the reference mark everything else is aligned to. The idea here is the same:
check every computation against something with a known answer before believing it, and write
results up so the document refuses to build when a claim stops holding.

The first release covers steady and oscillating flow in microfluidic channels — the kind of
problem often run in ANSYS Fluent or COMSOL — with free software, plus the general procedure
for reproducing a published simulation so that every step can be checked.

## New to Claude? Two terms

- **Skill** — a written procedure Claude reads before doing a task, plus any scripts it needs.
  You don't call it; Claude picks it up when your request matches what the skill is for.
- **Plugin** — a bundle of skills you install once. This repository is one plugin containing
  two skills.

## Get started

### In the Claude app or Cowork (web or desktop)

1. Open **Customize** in the left sidebar, then the **Plugins** tab.
2. Click **Add**, choose **Add marketplace**, and enter `s27183/fiducial-lab`.
3. Find **fiducial-lab** in the list and click **Add**.

The plugin is saved to your account, so it also appears in Claude Code when you sign in with
the same account.

### In Claude Code (terminal)

```
/plugin marketplace add s27183/fiducial-lab
/plugin install fiducial-lab@fiducial-lab
```

Then run `/reload-plugins`, or start a new session.

### In Claude Science

Claude Science adds skills rather than plugins. Open **Customize**, choose to add skills from a
repository, and enter `s27183/fiducial-lab`. If it asks which folder, point it at
`skills/` — both skills live there, each with its `SKILL.md` and, for the solver, its
`scripts/` folder.

Claude Science runs Python on your own machine, so both skills, solvers included, work there.

### Without Claude

```bash
pip install -r requirements.txt
python3 skills/channel-flow-fem/scripts/benchmarks.py all       # ~10 seconds
```

Tested with Python 3.12; older versions are untested. [Typst](https://github.com/typst/typst) is needed only to build the PDF.

## What to ask first

Once installed, try these in a new conversation:

> Run the fiducial-lab benchmarks and tell me whether they pass.

> Using fiducial-lab, compute the wall shear stress in a 500 µm channel with six 100 µm
> ridges at 5 mL/min, and show me the mesh-convergence check.

> I want to reproduce the simulation in [paper or DOI]. Walk me through recovering its
> parameters before we build anything.

The third is the real use. Claude will follow the `verified-reproduction` procedure: find the
paper's parameters in the publisher's source, flag anything missing or inconsistent, and set up
a check against a known answer before the actual geometry.

## Where it runs

| Skill | Claude app chat | Cowork | Claude Code / Claude Science | Without Claude |
|---|---|---|---|---|
| `verified-reproduction` — the method | yes | yes | yes | read `SKILL.md` |
| `channel-flow-fem` — the solvers | needs Python with pip | needs Python with pip | yes | yes |

The solvers need `pip install` of `scikit-fem` and `gmsh`, so they need an environment where
Claude can run Python and install packages — your own machine through Claude Code or Claude
Science is the dependable choice. **Tested so far:** the scripts, run directly. **Not yet tested:** loading the
plugin through each Claude client. Reports of either working or failing are welcome as issues.

## What is inside

| Skill | What it does |
|---|---|
| `verified-reproduction` | The procedure: recover a paper's parameters from its structured source, verify against a closed-form case, converge on the reported quantity, compare without tuning toward the published answer, and write up with build-time checks. Not specific to any field. |
| `channel-flow-fem` | Steady and oscillating 2-D channel flow with scikit-fem and gmsh — wall shear stress, cavity vortices, Womersley flow, TAWSS and OSI — with closed-form benchmarks that run in seconds. |

A fully worked reproduction of a published study will be added as an example.

## Scope, honestly

The solver skill rests on two closed-form benchmarks, so its verified scope is 2-D flow in
channels, steady and oscillating, up to a Reynolds number of about 170. Nothing here yet covers
3-D flow, flexible walls, non-Newtonian fluids, solid mechanics, imaging or signals. Skills are
added as real problems arrive.

## Contributing

New skills are welcome on one condition: each ships with a check against something whose
answer is known independently, runnable as a single command. `verified-reproduction` explains
what that means in practice.

## Licence

MIT — see `LICENSE`. Copyright © 2026 Son Tran.
