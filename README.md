# DISCOVER supplementary notebooks

The [Jupyter](http://jupyter.org) notebooks below contain all the code required to reproduce the figures and results of the paper *A novel independence test for somatic alterations in cancer shows that biology drives mutual exclusivity but chance explains most co-occurrence*.

To work with these files, Jupyter, IPython, and several Python packages should be installed. The easiest way to install these dependencies is by using [Miniconda](http://conda.pydata.org/miniconda.html) or [Anaconda](https://www.continuum.io/why-anaconda). The following command creates a conda environment that contains all required packages to execute the notebooks.

```bash
conda create -n discover-notebooks -c http://ccb.nki.nl/software/discover/repos/conda \
    corclust==0.1 \
    discover==0.9 \
    matplotlib==1.5.1 \
    networkx==1.11 \
    numpy==1.10.4 \
    pandas==0.17.1 \
    pytables==3.2.2 \
    scipy==0.17.0 \
    statsmodels==0.6.1 \
    notebook \
    ipykernel
```

Only for the notebook named *Group test* a few more packages need to be installed using the following command.

```bash
conda install -n discover-notebooks -c http://ccb.nki.nl/software/discover/repos/conda -c r -c msys2 \
    switching==0.1 \
    ccomet-with-timeout==1.0.2 \
    rpy2 \
    ipyparallel
```
	
Next, activate the created environment and start the Jupyter notebook using the following two commands. Make sure `<notebook-dir>` is replaced by the location of the .ipynb files after unzipping the downloaded file.
	
```bash
source activate discover-notebooks
jupyter notebook --notebook-dir=<notebook-dir>
```
	
On Windows, the first command should be replaced by:
	
```bash
activate discover-notebooks
```


## Simulated data analyses

* [Pairwise analyses of simulated data](notebooks/Pairwise%20analyses%20of%20simulated%20data.ipynb)

  Compares the Binomial, Fisher's exact and DISCOVER tests on simulated data.


* [Group test](notebooks/Group%20test.ipynb)

  Compares the DISCOVER group test to six alternative methods (CoMEt, MEGSA, MEMo, muex, mutex, and TiMEx) on simulated data.


## Pan-cancer analyses

* [Download PanCan12 data](notebooks/Download%20PanCan12%20data.ipynb)

  Downloads the mutation and copy number data for the TCGA PANCAN12 studies.


* [Gene selection](notebooks/Gene%20selection.ipynb)

  Selects the genes for use in the pairwise analyses.


* [Pairwise analysis](notebooks/Pairwise%20analysis.ipynb)

  Performs pairwise co-occurrence and mutual exclusivity analyses.


* [Within-chromosome co-occurrence analysis](notebooks/Within-chromosome%20co-occurrence%20analysis.ipynb)

  Tests for co-occurrences between genes located on the same chromosome, in order to assess whether the DISCOVER test will detect these 'positive controls'.


* [STRING enrichment](notebooks/STRING%20enrichment.ipynb)

  Determines the overlap of mutually exclusive gene pairs with the STRING functional interaction network.


* [MSigDb group tests](notebooks/MSigDb%20group%20tests.ipynb)

  Identifies significantly mutually exclusive gene sets based on predefined gene sets extracted from MSigDb.


* [De novo gene set identification](notebooks/De%20novo%20gene%20set%20identification.ipynb)

  Detects de novo mutually exclusive gene sets based on correlation clustering of pairwise mutual exclusivities.
