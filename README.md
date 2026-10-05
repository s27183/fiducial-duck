# fiducial-duck

AI that does the computational heavy lifting in research, and checks its own work. Claude skills, plus plain Python that runs without Claude.

A *fiducial* is the reference mark everything else is aligned to. The duck is the programmer's
rubber duck: the one you explain your code to until you spot your own mistake. Together they
are the idea here. Every computation is anchored to something with a known answer, and checked
against it, before anyone builds on it.

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
2. Click **Add**, choose **Add marketplace**, and enter `s27183/fiducial-duck`.
3. Find **fiducial-duck** in the list and click **Add**.

The plugin is saved to your account, so it also appears in Claude Code when you sign in with
the same account.

### In Claude Code (terminal)

```
/plugin marketplace add s27183/fiducial-duck
/plugin install fiducial-duck@fiducial-duck
```

Then run `/reload-plugins`, or start a new session.

### In Claude Science

1. Open **Customize**, then **Skills**, then **Import from GitHub**.
2. Enter `s27183/fiducial-duck` and click **Preview**. Both skills are listed.
3. Click **Import 2 skills**.

Imported skills don't update automatically; import again to pick up a new release. Claude
Science runs Python on your own machine, so both skills, solvers included, work there.

### Without Claude

```bash
pip install -r requirements.txt
python3 skills/channel-flow-fem/scripts/benchmarks.py all       # ~10 seconds
```

Tested with Python 3.12; older versions are untested. [Typst](https://github.com/typst/typst) is needed only to build the PDF.

## What to ask first

Once installed, try these in a new conversation:

> Run the fiducial-duck benchmarks and tell me whether they pass.

> Using fiducial-duck, compute the wall shear stress in a 500 µm channel with six 100 µm
> ridges at 5 mL/min, and show me the mesh-convergence check.

> I want to reproduce the simulation in [paper or DOI]. Walk me through recovering its
> parameters before we build anything.

The third is the real use. Claude will follow the `verified-reproduction` procedure: find the
paper's parameters in the publisher's source, flag anything missing or inconsistent, and set up
a check against a known answer before the actual geometry.

## Where it runs

| Skill | Cowork | Claude Science | Claude Code | Claude app chat | Without Claude |
|---|---|---|---|---|---|
| `verified-reproduction` — the method | yes | yes | yes | yes | read `SKILL.md` |
| `channel-flow-fem` — the solvers | yes, tested | yes | yes | needs Python with pip | yes, tested |

The solvers need `pip install` of `scikit-fem`, `gmsh` and `meshio`, so they need an environment where
Claude can run Python and install packages.

**Tested so far:**

- **Cowork**, installed through the plugin marketplace. A request that never named the plugin
  ("compute the wall shear stress in a 500 µm channel at 6 mL/min and check it against the
  analytical value") loaded `channel-flow-fem` unprompted, installed the dependencies, and
  matched the exact 0.336 Pa to within 10⁻¹³. Its sandbox also needed the X11/OpenGL system
  libraries for `gmsh`, and `meshio`, which v0.1.0 failed to list.
- **Claude Science**, imported from GitHub. Both skills were found and imported intact, and the
  benchmarks pass from the imported copy.
- **The scripts**, run directly from a fresh download of the release.

**Not yet tested:** installing through Claude Code and using the solvers from Claude app chat.
Reports of either working or failing are welcome as issues.

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

MIT, which in plain terms means you can use, change and share this freely, in research,
teaching, papers or anything else, with no permission needed. The only condition is keeping
the copyright notice in copies of the code. Full text in `LICENSE`. Copyright © 2026 Son Tran.
