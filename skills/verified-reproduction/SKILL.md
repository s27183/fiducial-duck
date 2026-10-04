---
name: verified-reproduction
description: "Procedure for reproducing a published simulation or computational result so that each step can be checked: recovering parameters from the publisher's structured source, verifying against a closed-form case, mesh or resolution convergence, comparing with the published result without tuning toward it, and writing it up so the document fails to build when a claim stops holding. Use when replicating a paper's CFD, FEA or numerical model, porting a model from ANSYS or COMSOL to open-source tools, or auditing whether a computed result supports a written claim. Part of the fiducial-duck plugin."
---

# Verified reproduction

A reproduction is worth something only if every link can be checked by someone who was not
there. The order below matters: each step protects the next from a class of error that is
plausible, internally consistent and wrong.

## 1. Recover the parameters from the structured source, and list the gaps

Read geometry, material properties and boundary conditions from the publisher's tagged
full text (JATS XML from PubMed Central or the publisher) rather than a PDF-to-text rendering.
Plain-text extraction silently drops superscripts: `7 × 10<sup>−4</sup> Pa s` becomes
"7 × 10−4" or worse, and a viscosity off by eight orders of magnitude still produces a
plausible-looking field.

For every parameter record the value, where it was found, and whether it is stated once or
several times. Then search for each dimension a second time. Papers contradict themselves —
the Methods and a figure caption giving different values for the same channel height is a
common pattern, and it is not cosmetic: wall shear at a fixed flow rate scales as 1/h², so a
20 % difference in height changes shear by a factor of 1.44, and by more inside cavities.

When a value conflicts, check each candidate against any control measurement the paper reports
(a ridge-free channel at a stated flow rate and shear, say). Report the resolution as an
inference from consistency, not a fact, and carry the other reading as a sensitivity case.

Supplementary material is part of the source. Fetch it before concluding a value is absent.

## 2. Verify the solver against a closed form, in the target's units and regime

Before the real geometry, solve a case with an exact answer at the same channel size, fluid
properties and flow rate: Poiseuille flow for steady channels, Womersley flow for oscillatory
ones, a cantilever or Lamé cylinder for solid mechanics. Report the relative error as a number.

If the discretisation can represent the exact solution, the error should be at machine
precision; anything larger is a bug. Otherwise it should converge at the scheme's design order.

## 3. Converge on the quantity you will report

Refine until the reported quantity stops changing — peak positions, peak shear, a count — and
report the change between the last two levels. Three levels is the minimum that shows a trend.
Converging the field and reporting a derived quantity is not the same thing.

## 4. Compare with the published result, and never tune toward it

Only now compare. When the result disagrees, ask which **unstated choice** could account for
it before changing anything: a peak-detection threshold, an averaging window, which wall was
sampled, which cavity along the channel. Sweep that choice and report the band over which the
published result is recovered. "Matches at an 8 % prominence criterion, and at any criterion
between 5 and 8 %" is a finding; "matches" is a claim.

If the agreement depends on a choice the paper does not state, say so in the write-up. That is
a contribution to the record, not a weakness of the reproduction.

## 5. Make the write-up regenerate its own numbers and refuse to lie

Write the document as a template whose every number is substituted from the results files at
build time. Then add assertions for each load-bearing sentence — not only the headline
sequence, but the specific per-case claims the prose makes. The build must fail when a claim
stops holding. Assertions of this kind catch, among other things, a results table whose
"computed" column was quietly filled from the published values, so it could never disagree.

Check every reference against Crossref or PubMed before it enters the bibliography: author
lists, volume, pages and DOI. Identifiers recalled from memory are wrong often enough to matter.

## 6. State the boundary of the checking

Separate what was verified mechanically (code against closed forms, document against results,
references against registries) from what needs a domain expert (whether the 2-D mid-plane
idealisation, the rigid walls or the inlet assumption suit the biology). Name the class of
error the checking leaves open rather than adding a blanket disclaimer.

## Failure modes observed in practice

- Plain-text extraction losing an exponent in a material property.
- A dimension taken from a figure caption while the control case was validated against the
  Methods value. Every automated check passed, because each was internally consistent; a human
  reader re-checking the source caught it.
- A table column described as computed but populated from the paper.
- Fabricated or misattributed identifiers (DOIs, PMIDs, first authors) in drafts.
- "Laminar due to low Reynolds number" read as licence for Stokes flow, which removes the
  very vortices being studied.
