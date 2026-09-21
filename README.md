# Bounded Adaptive Coherence: a coupling-tensor λ_min criterion as a computational object for aging-versus-cancer failure modes

**Thesis #18.** Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

Can the Principle of Bounded Adaptive Coherence — λ_min(C(t)) against scale-wise entropy production — be stated as a falsifiable computational object such that aging-like global coupling decay and cancer-like selective organism-scale decoupling appear as distinct sectors of the same tensor, without claiming clinical rejuvenation or cure?

The literal minimum eigenvalue of a raw nonnegative matrix is not that object. On C(ρ) = (1−ρ)I + ρ11ᵀ, λ_min falls as ρ rises, while algebraic connectivity of the Laplacian rises. The criterion in the manuscript is λ_min of the grounded Laplacian of the symmetric weights, minus a declared dimensionless proxy. Five closed-form paths run on one 5×5 toy. Global off-diagonal decay is classed aging-like. Decay confined to the edges that touch the organism index is classed cancer-like. An uneven decay, a diagonal-only decay, and a held matrix are classed neither. A negative margin is not by itself a sector label.

The index names are row labels. They are not assays. No entropy production was computed. This deposit does not claim rejuvenation or cure.

A neighbouring deposit asks when load×gain coupling produces a Gompertz-like hazard (Thesis #13). That question is about hazard shape. It is not this spectral margin, and no number here is taken from that deposit.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `sim/bac_toy.py` | Seeded toy spectrum and classifier (seed 20260921) |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Eigenvalue paths, sectors, Fiedler mass, margin |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/bac_toy.py
python3 build_pdf.py
```

NumPy and Matplotlib are required for the toy. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json` and `sim/figures/`.

## Cite

Ogbonna KE. Bounded adaptive coherence: a coupling-tensor λ_min criterion as a computational object for aging-versus-cancer failure modes [Internet]. Thesis #18 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-18-bounded-adaptive-coherence

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.
