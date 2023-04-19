AbSplice: aberrant splicing prediction across human tissues
--------------------------------
This package predicts aberrant splicing across human tissues with AbSplice. 
AbSplice predictions are based on enhanced tissue-specific splice site annotations ([SpliceMaps](https://github.com/gagneurlab/splicemap)).
If purely sequence based information is available, different DNA-based splicing predictions are combined into an integrative model "AbSplice-DNA".
Integration of RNA-seq data from an accessible tissue (e.g. blood or skin) of an individual to predict aberrant splicing in any other tissue from the same individual is supported in "AbSplice-RNA".
Genome-wide AbSplice-DNA scores for all possible SNVs are available [here](https://doi.org/10.5281/zenodo.6408331) for download.

## Installation
Clone git repo:
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
Install modules from absplice:
```
pip install -e .
```

## Example usecase
The [example](https://github.com/gagneurlab/splicing-outlier-prediction/tree/master/example) folder contains a snakemake workflow to generate AbSplice predictions, given a vcf file and a fasta file (will be downloaded if not provided).
The snakemake workflow will download precomputed SpliceMaps from Zenodo and run AbSplice based on these annotations.
To generate predictions run:
```
cd example
python -m snakemake -j 1 --use-conda
```
To run this example on your own data, specify the genome version that you are going to use (hg19 is compatible with gtex_v7 and hg38 is compatible with gtex_v8) in the field 'genome' of the [config](https://github.com/gagneurlab/splicing-outlier-prediction/tree/master/example/config.yaml) file and store all vcf files for analysis in data/resources/vcf_files/.

In the field 'splicemap_tissues' of the config file you can uncomment the tissues that AbSplice will use to generate predictions.

If you want to run the example on large datasets, you can enable a fast lookup interface spliceai_rocksdb that uses precomputed scores. This precomputed database will be downloaded from Nextcloud (it will take significant time – about 3-4 hours). To enable fast lookup simply change the field 'use_rocksdb' in config file to True.