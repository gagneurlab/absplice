Splicing Outlier Prediction
--------------------------------

## Installation
Clone git repository of splicing_outlier_prediction:
```
git clone git@gitlab.cmm.in.tum.de:gagneurlab/splicing-outlier-prediction.git
```

cd into repo directory:
```
cd splicing-outlier-prediction
```

Install conda environment:
```
mamba env create -f environment.yaml
```
Install splicing-outlier-prediction into conda environment:
```
conda activate absplice
pip install -e .
```

Clone extra git repository (SpliceMaps):
```
git clone https://github.com/gagneurlab/splicemap.git
```

Install splicemap into conda environment:
```
conda activate absplice
cd splicemap
pip install -e .
```

## Example usecase
The ./example folder contains a snakemake workflow to generate AbSplice predictions, given a vcf file and a fasta file.
The snakemake workflow will download precomputed SpliceMaps from zenodo and run AbSplice based on these annotations.
To generate predictions run:
```
cd example
snakemake -m python -j 1
```