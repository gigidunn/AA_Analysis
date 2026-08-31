# AA_Analysis

End-to-end pipelines to fine-tune, run, and evaluate two archetypal analysis methods, ParTI (linear) and MIDAA (deep generative), on single-cell transcriptomic data. MIDAA is fine-tuned via its API with custom hyperparameters and treatment-label supervision; statistical model selection and biological interpretation are original evaluation pipelines built on top of each tool's output, not a one-off comparison.

## Background

Archetypal analysis represents cell populations as mixtures of extreme phenotypes ("archetypes") rather than discrete clusters, often a better fit for continuous biological processes including treatment response (and development of resistance). This repo builds the pipeline to take a lab's transcriptomic data through either method end to end:

- **ParTI**: the established linear method (Hart et al. 2015, Uri Alon Lab), with randomisation-based statistical significance testing.
- **MIDAA**: a variational autoencoder based method ([sottorivalab](https://sottorivalab.github.io/midaa/scMulti_multimodal.html)), supporting supervised or unsupervised, raw or normalised input, with model selection across a range of archetype counts via ELBO.

## Pipeline structure

**1. Data** (`Data/`)
Seurat-based preprocessing of the BacDrop single-cell dataset across three antibiotic treatments (meropenem, ciprofloxacin, gentamicin): QC, filtering, normalisation, clustering, and marker detection. Also handles class balancing across treatments and generates scrambled null data for downstream significance testing, then converts the processed matrix into the input formats each tool needs.

**2. MIDAA** (`MIDAA/`)
Runs MIDAA across a range of archetype counts, supervised on treatment labels, with ELBO-based model selection to identify the archetype count that best balances statistical fit against interpretability. Post-processing covers both statistical evaluation (ELBO, t-ratio between polytope volume and the convex hull of the data) and biological interpretation (mapping archetypes to gene expression and treatment contribution).

**3. ParTI** (`ParTI/`)
Runs ParTI with randomisation controls to test statistical significance of the resulting polytope. Post-processing covers 3D polytope visualisation, t-ratio, density of cells and marker genes relative to each archetype, and comparison against MIDAA's archetype-gene weightings.

## Dependencies

**Python** (MIDAA): `torch`, `numpy`, `scanpy`, `matplotlib`, `scipy`, `scikit-learn`, `midaa`
**R** (Data preprocessing): `Seurat`, `data.table`, `dplyr`, `ggplot2`, `Matrix`, `patchwork`
**MATLAB** (ParTI): ParTI toolbox (Uri Alon Lab), not included, download separately

## Data

BacDrop dataset: https://pubmed.ncbi.nlm.nih.gov/36708705/

## Context

Built during an RA position in the MIT Bioengineering Department (Lauffenburger Lab), evaluating frontier model architectures for interpretable analysis of high-dimensional biological data. This work accompanies a research paper further discussing and evaluating the application of bioinformatic tools for biological interpretation. 

## Citation

This repo builds on the following tools:
- Hart, Y. et al. (2015). ParTI: a toolbox for Pareto Task Inference.
- Milite, S., Caravagna, G., Sottoriva, A. (2025). MIDAA. *Genome Biology*.

## License

MIT

## Contact

Gigi Dunn — https://www.linkedin.com/in/gigi-dunn-/
