Library for processing nanobody (VHH) sequencing data.

This package is intended to be used by the `nbseq-workflow`, a Snakemake workflow for processing raw sequencing from Phage-seq experiments into feature tables. This package can then be used interactively to query, calculate, and visualize the resulting data structures.

## Package organization

- `nbseq`: `Experiment` class that collects and organizes data for a given experiment
- `utils`: utility functions
- `asvs`: process VHH sequences
- `ft`: read and process feature tables (sparse matrices of sample x feature [i.e. VHH, CDR3, etc.])
- `select`: perform calculations relevant to phage display selection (e.g. enrichment, amplification bias)
- `resynth`: choose and resynthesize recombinant VHH genes as fragments
- `norm`: normalize feature table data to remove effect of variable library sizes
- `prep`: utilities to aid in sample preparation
- `ordination`: perform ordination/dimensionality reduction on feature tables
- `predict`: machine learning prediction on feature tables
- `cloning`: simulate cloning recombinant VHHs into destination vectors
- `design`: create design matrices for inference and machine learning
- `pheno`: compare and visualize phenotypes of samples
- `msa`: perform multiple sequence alignment
- `viz`: visualization

## Dependencies

Later versions may work but have not been tested.

Required and recommended dependencies can be installed using `conda` via the included `environment.yml` file

- Required dependencies:
    - anndata=0.9.2
    - biom-format=2.1.15
    - humanize=4.7.0
    - natsort=8.4.0
    - numpy=1.24.4
    - pandas=2.0.3
    - pysam=0.21.0
    - pyyaml=6.0
    - scikit-bio=0.5.9
    - scipy=1.10.0
    - statsmodels=0.14.0
- Optional dependencies:
    - For machine learning:
        - scikit-learn=1.3.0
        - scikit-optimize=0.9.0
        - xgboost=1.5.1
    - For database-accelerated feature queries:
        - connectorx=0.3.1
        - mmseqs2=14.7e284
        - sqlalchemy=2.0.19
        - sqlite=3.42.0
    - For recombinant sequence optimization and cloning:
        - dna_features_viewer=3.1.2
        - dnachisel=3.2.11
        - pydna=3.1.0
        - python-codon-tables=0.1.12
    - For processing Sanger sequencing chromatograms:
        - bioconvert=1.1.1
    - For visualizations:
        - altair==5.1.0.dev0
        - logomaker=0.8
        - matplotlib=3.7.2
        - plotly=5.16.0
        - pygments=2.16.1
        - plotnine=0.12.2
        - seaborn=0.12.2
        - patchworklib=0.6.3
        - pip: mnemonicode=1.4.5
    - For interactive "dashboard" visualizations:
        - altair-transform=0.2.0
        - bokeh=3.2.2
        - ipykernel=6.25.1
        - ipywidgets=8.1.0
        - panel=1.2.1
    - For normalization using `scran` package:
        - r=4.1
        - bioconductor-biomformat=1.22.0
        - bioconductor-scran=1.22.1