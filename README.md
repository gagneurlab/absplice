AbSplice: aberrant splicing prediction across human tissues
--------------------------------
This package predicts aberrant splicing across human tissues with AbSplice. 
AbSplice predictions are based on enhanced tissue-specific splice site annotations ([SpliceMaps](https://github.com/gagneurlab/splicemap)).
If purely sequence based information is available, different DNA-based splicing predictions are combined into an integrative model "AbSplice-DNA".
Integration of RNA-seq data from an accessible tissue (e.g. blood or skin) of an individual to predict aberrant splicing in any other tissue from the same individual is supported in "AbSplice-RNA".
Genome-wide AbSplice-DNA scores for all possible SNVs are available [here](https://doi.org/10.5281/zenodo.6408331) for download (currently 10,000 scores for genes uploaded).

## Installation
Clone git repository of splicing_outlier_prediction:
```
git clone https://github.com/gagneurlab/absplice.git
```

cd into repo directory:
```
cd absplice
```

Install conda environment:
```
# Recommended if you have mamba installed
mamba env create -f environment.yaml
# otherwise
conda env create -f environment.yaml
```
Activate conda environment:
```
conda activate absplice
```

## Example usecase
The [example](https://github.com/gagneurlab/splicing-outlier-prediction/tree/master/example) folder contains a snakemake workflow to generate AbSplice predictions, given a vcf file and a fasta file.
The snakemake workflow will download precomputed SpliceMaps from zenodo and run AbSplice based on these annotations.
To generate predictions run:
```
cd example
python -m snakemake -j 1
```
To run this example on your own data, simply change the config:
* vcf: path to your vcf file
* fasta: provide url to download fasta file (e.g. from Gencode or Ensembl website)
* splicemap_dir: directory where precomputed SpliceMaps from Zenodo will be downloaded