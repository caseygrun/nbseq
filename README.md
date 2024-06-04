Library for processing nanobody (VHH) sequencing data.

This package is intended to be used in two ways:

1) by the Snakemake workflows in the [`phage-seq` repository](http://github.com/caseygrun/phage-seq) to batch process raw sequencing data from Phage-seq experiments into feature tables. For this usage, the Snakemake workflows will install `nbseq` automatically as needed at the appropriate steps/ 

2) interactively within Jupyter notebooks to query, calculate, and visualize the resulting data structures.


To explore the code used in our paper, start with the [`phage-seq` repository](http://github.com/caseygrun/phage-seq). That repository also includes several demonstration notebooks and datasets to explore the functionality of this library. Follow the instructions there to create or obtain an example dataset, then return to this repository for instructions on how to install `nbseq` for interactive analysis.

## Installation and usage

1. First, perform preprocessing of raw data using the Snakemake workflow(s) in the [`phage-seq` repository](http://github.com/caseygrun/phage-seq), following instructions there. The relevant steps within the Snakemake workflows will install the `nbseq` package; it is not necessary to manually install the `nbseq` package for this step.

2. Second, for interactive analysis, it is recommended to create a dedicated `conda` environment for use with the `nbseq` package. 

    1.  If you have not already done so, install the [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html#mamba-install) (or [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)) package manager. I recommend using the [`miniforge`](https://github.com/conda-forge/miniforge) distribution. 
    2.  Create and activate a new `conda` environment for `nbseq` and its dependencies. You have two options:

        - Minimal installation: installs only the required core dependencies:

                wget https://github.com/caseygrun/nbseq/raw/main/environment-min.yaml
                conda env create -f environment-min.yaml
                conda activate nbseq-min

        - Full installation of all optional dependencies:

                wget https://github.com/caseygrun/nbseq/raw/main/environment.yaml
                conda env create -f environment.yaml
                conda activate nbseq

        In both cases, you do not need to clone this repository. You only need to download the `.yaml` file(s) using the steps above; the remaining files will be downloaded and installed by `conda`. 

    3.  [Install JupyterLab](https://jupyterlab.readthedocs.io/en/stable/getting_started/installation.html#conda) or the Jupyter Notebook, if you have not already; you also have two choices for this:

        - I recommend creating a separate dedicated `conda` environment for JupyterLab and using the [`nb_conda_kernels`](https://github.com/anaconda/nb_conda_kernels) package; this will allow you to install and update JupyterLab separately from `nbseq` and its many dependencies; the `nb_conda_kernels` package lets you access the `nbseq` environment (and any other conda environments you create) from within JupyterLab:
        
                conda deactivate 
                conda create -n jupyter jupyterlab nb_conda_kernels panel
                conda activate jupyter

        - Alternatively, you can install JupyterLab directly into the same environment as `nbseq`:
        
                conda install jupyterlab

    4.  Launch JupyterLab and follow the instructions below in "Usage:"

            jupyter lab


Note: `nbseq` is tested only on 64-bit Linux.

## Usage

The main entry point for interactive analysis is the `nbseq.Experiment` class, which loads and organizes feature tables, phylogenetic trees, metadata, and databases for a given experiment. `nbseq.Experiment.from_files` can load data from the directory structure created by the [`phage-seq`](http://github.com/caseygrun/phage-seq) Snakemake workflows. Consult the docstring `?nbseq.Experiment.from_files` for a more detailed description of the options.


    >>> import nbseq
    >>> ex = nbseq.Experiment.from_files(
    ...     # skip loading the larger `aa` (e.g. each VHH amino acid sequence is a 
    ...     # distinct column) feature table and # phylogenetic tree; by default, 
    ...     # the function loads the `cdr3` and `aa` feature tables
    ...     ft_aa=None, tree_aa=None, 
    ...     metadata='config/metadata_full.csv') #'intermediate/cdr3/features/all/alpaca/asvs.nwk')
    Loading experiment panning-extended from '/vast/palmer/home.mccleary/cng2/code/phageseq-paper/panning-extended'...
    - Reading metadata from config/metadata_full.csv ...
    - Reading phenotypes from config/phenotypes.csv ...
    - Reading Config from config/config.yaml ...
    - Using SQL database at 'sqlite:////vast/palmer/home.mccleary/cng2/code/phageseq-paper/panning-extended/intermediate/aa/asvs.db'
    - Reading feature data for table 'cdr3' from results/tables/cdr3/asvs.csv (2.6 MB)...
    - Reading aa feature table from results/tables/aa/feature_table.biom (350.4 MB)...
    - Reading cdr3 feature table from results/tables/cdr3/feature_table.biom (8.4 MB)...
    - Warning: phylogeny for space 'aa' at 'intermediate/aa/features/top_asvs/alpaca/asvs.nwk' does not exist!
    - Warning: phylogeny for space 'cdr3' at 'intermediate/cdr3/features/top_asvs/alpaca/asvs.nwk' does not exist!
    - Using mmseqs2 database 'aa' at 'intermediate/aa/features_db/features'
    - Warning: mmseqs2 database for space 'cdr3' at 'intermediate/cdr3/features_db/features' does not exist!
    - Reading enrichment model (conditional ECDF) for space cdr3 from results/tables/cdr3/enrichment/null/ecdf.pickle (307.6 kB)...
    Finished in 20.29 seconds


Displaying the `Experiment` object shows a summary:

    >>> ex
    Experiment('panning-extended') with feature spaces ['aa', 'cdr3']:
    obs: ['plate.x' 'well.x' 'depth' 'expt' 'round' 'sample' 'phage_library'
        'notes' 'r' 'io' 'kind' 'selection' 'replicate' 'name_full' 'name'
        'well_027e' 'sel_plate_027i' 'sel_well_027i' 'selection_027j' 'plate.y'
        'well.y' 'category' 'antigen' 'genotype_pair' 'gene_CS' 'gene_S'
        'genotype_CS' 'background_CS' 'strain_CS' 'loc_CS' 'cond_CS' 'genotype_S'
        'background_S' 'strain_S' 'loc_S' 'cond_S' 'cond_notes' 'bflm' 'swim'
        'twitch' 'swarm' 'PMB-R' 'FEP-R' 'TET-R' 'CIP-R' 'CHL-R' 'GEN-R' 'ERY-R'
        'IPM-R' 'cdiGMP' 'FliC' 'FliCa' 'FliCb' 'FlgEHKL' 'PilQ' 'PilA' 'PilB'
        'LasA' 'LasB' 'Apr' 'XcpQ' 'ToxA' 'EstA' 'LepA' 'PlpD' 'Phz' 'Pcn' 'Pvd'
        'Hcn' 'Rhl' 'T3SS' 'T6SS' 'Pel' 'Psl' 'CdrB' 'SCV' 'Mucoid' 'Alginate'
        'OprM' 'OprJ' 'OprN' 'OprOP' 'OpdH' 'OprD' 'OprL' 'OprF' 'OprG' 'OprH'
        'OprB' 'MexAB' 'MexCD' 'MexEF' 'MexJK' 'MexXY' 'MexGHI' 'PirA' 'Pfu'  'TonB'
        'FptA' 'FpvA' 'PfeA' 'CupB5' 'CupA' 'CupB' 'CupC' 'CupD'  'LPS-LipidA-
        Palmitoyl' 'L-LipidA-Ara4N' 'LPS-CPA' 'LPS-OSA' 'LPS-galU'  'LPS-rough'
        'LPS' 'description']
    - aa      : 439 samples x 5134305 features, database: None
    var: ['reads' 'nsamples']
    - cdr3    : 439 samples x 40292 features, database: None
    var: ['reads' 'nsamples']
    SQL: sqlite:////vast/palmer/home.mccleary/cng2/code/phageseq-paper/panning-extended/intermediate/aa/asvs.db


From there, you can access various visualizations via the experiment visualizer, `ex.viz`, e.g.:

    >>> ex.viz.top_feature_barplot(f"expt == '027j' & FlgEHKL == 1", select_from_round=None, n=100).facet(column='selection')

Or load additional interactive visualizations using the `nbseq.viz` package, e.g.

    >>> import nbseq.viz.dash
    >>> nbseq.viz.dash.selection_group_dashboard(
    ...     ex, starting_phenotype='FlgEHKL', 
    ...     global_query=(
    ...         "expt == '027j' & io == 'i' & kind == '+'")
    ... )

See the [`phage-seq` repository](http://github.com/caseygrun/phage-seq) for additional examples: [`panning-minimal`](http://github.com/caseygrun/phage-seq/panning-minimal/workflow/analysis.ipynb) and [`panning-extended`](http://github.com/caseygrun/phage-seq/panning-extended/workflow/analysis.ipynb).

## Package organization

The `nbseq` package contains the following sub-modules:

- `nbseq`: `Experiment` class that collects and organizes data for one or more Phage-seq experiments. Namely, `Experiment` loads and organizes trees, metadata, and feature tables in multiple feature spaces (e.g. VHH, CDR3, etc.) and facilitates projecting between them. `Experiment` also provides an interface for interactive visualization of the entire experiment or subsets thereof.
- `utils`: utility functions
- `asvs`: process VHH sequences: calculate residue frequencies, consensus sequences, query for similar sequences, project between feature spaces (e.g. CDR3 counts to full length amino acid sequence counts)
- `ft`: read and process feature tables (sparse matrices of _sample_ x _feature_ [i.e. VHH, CDR3, etc.])
- `select`: perform calculations relevant to phage display selection (e.g. enrichment, amplification bias); calculate null models of enrichment probabilities
- `norm`: normalize feature table data to remove effect of variable library sizes
- `ordination`: perform ordination/dimensionality reduction on feature tables
- `design`: create design matrices for inference and machine learning
- `pheno`: compare and visualize phenotypes of samples
- `msa`: perform multiple sequence alignment with `mafft`
- `viz`: generate various visualizations: feature bar plots, rank-abundance curve (Whittaker plots), abundance curves, 2D/3D ordination plots, sequence logos, receiver-operator characteristic curves, etc.
- `predict`: perform machine learning prediction on feature tables
- `resynth`: choose and resynthesize recombinant VHH genes as gene fragments. Includes routines for identifying consensus sequences, trimming and adding adapter sequences, etc. 
- `cloning`: simulate cloning recombinant VHHs into destination vectors
- `prep`: utilities to aid in HTS library preparation

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