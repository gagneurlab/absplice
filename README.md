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

## Output

Output of AbSplice is an tabular data which contains following columns:


|     ID     | Column | Description |
|  --------  | ----- | ----------- |
|  `variant` | Variant | ID string of the variant. |
|  `gene_id` | GeneID | Ensembl GeneID of the gene which the variant belongs to. |
|  `tissue`  | Tissue | Name of the tissue that was used from SpliceMap. |
| `AbSplice_DNA` | AbSplice-DNA | The AbSplice score is a probability estimate of how likely aberrant splicing of some sort takes place in a given tissue and reports the splice site with the strongest effect. The model was trained using scores from MMSplice and SpliceAI models as well as annotation features from tissue-specific SpliceMaps. To ease downstream applications we suggest three cutoffs (high: 0.2, medium: 0.05, low: 0.01), which approximately have the same recalls as the high, medium and low cutoffs of SpliceAI. |
| `delta_score` | SpliceAI DeltaScore | Input feature to AbSplice. The main score predicted by SpliceAI, computed as a maximum of acceptor gain, acceptor loss, donor gain, and donor loss delta scores. The score represents probability of the variant being splice-altering. |
| `delta_logit_psi` | MMSplice + SpliceMap score | Input feature to AbSplice. The score is computed by using SpliceMap as an annotation for MMSplice. The score shows the effect of the variant on the inclusion level (PSI – percent spliced in) of the junction. The score is on a logit scale. If the score is positive, it shows that variant leads higher inclusion rate for the junction. If the score is negative, it shows that variant leads higher exclusion rate for the junction. |
| `delta_psi` | MMSplice + SpliceMap + Ψ_ref score | Input feature to AbSplice. The `delta_psi` (∆Ψ)   score is computed by converting `delta_logit_psi` (∆𝑙𝑜𝑔𝑖𝑡(Ψ)) to natural scale with the splicing scaling law and `ref_psi` (Ψ𝑟𝑒𝑓): <br>∆Ψ = σ(∆𝑙𝑜𝑔𝑖𝑡(Ψ) + 𝑙𝑜𝑔𝑖𝑡(Ψ𝑟𝑒𝑓)) − Ψ𝑟𝑒𝑓</br>|
| `splice_site_is_expressed` | Splice site expressed | Input feature to AbSplice. Binary feature indicating if the splice site is expressed in the target tissue using a cutoff of 10 split reads (`median_n` from SpliceMap). |
| `junction` | Intron | Coordinates of the respective intron from SpliceMap. |
| `event_type` | Type of splicing event | The `event_type` takes values either psi3 or psi5 for alternative donor or alternative acceptor site usage (from SpliceMap). |
| `splice_site` | Splice site location | Coordinates of the splice site. |
| `ref_psi` | Ψ reference score | Reference level of site usage from SpliceMap. |
| `median_n` | Median coverage Splice site | The median number of split reads sharing the splice site (from SpliceMap). |
| `acceptor_gain` | SpliceAI Delta score (acceptor gain) | Probability computed by SpliceAI that the variant will lead to an acceptor gain at `acceptor_gain_position` . |
| `acceptor_loss` | SpliceAI Delta score (acceptor loss) | Probability computed by SpliceAI that the variant will lead to an acceptor loss at `acceptor_loss_position`. |
| `donor_gain` | SpliceAI Delta score (donor gain) | Probability computed by SpliceAI that the variant will lead to a donor gain at `donor_gain_position`. |
| `donor_loss` | SpliceAI Delta score (donor loss) | Probability computed by SpliceAI that the variant will lead to a donor loss at `donor_loss_position`. |
| `acceptor_gain_position` | SpliceAI Delta postion (acceptor gain) | Delta position represents the location of respective splicing change  relative to the variant position: positive values are downstream of the variant, negative values are upstream. |
| `acceptor_loss_position` | SpliceAI Delta position (acceptor loss) | See description of `acceptor_gain_position`. |
| `donor_gain_position` | SpliceAI Delta postion (donor gain) | See description of `acceptor_gain_position`. |
| `donor_loss_position` | SpliceAI Delta position (donor loss) | See description of `acceptor_gain_position`. |

## Example usecase
The [example](https://github.com/gagneurlab/splicing-outlier-prediction/tree/master/example) folder contains a snakemake workflow to generate AbSplice predictions, given a vcf file and a fasta file (will be downloaded if not provided).
The snakemake workflow will download precomputed SpliceMaps from Zenodo and run AbSplice based on these annotations.
To generate predictions run:
```
cd example
python -m snakemake -j 1 --use-conda
```
To run this example on your own data do the following:

- Specify the genome version that you are going to use (hg19 is compatible with gtex_v7 and hg38 is compatible with gtex_v8) in the field `genome` of the [config](https://github.com/gagneurlab/splicing-outlier-prediction/tree/master/example/config.yaml) file.

- Store all vcf files for analysis to `data/resources/vcf_files/`.

Optionally:

- In the field `splicemap_tissues` of the [config](https://github.com/gagneurlab/splicing-outlier-prediction/tree/master/example/config.yaml) file you can uncomment the tissues that AbSplice will use to generate predictions (by default only fibroblasts).

- If you want to run the example on large datasets, you can enable a fast lookup interface spliceai_rocksdb that uses precomputed scores. This precomputed database will be downloaded from Nextcloud (it will take significant time – about 3-4 hours). To enable fast lookup simply change the field `use_rocksdb` in [config](https://github.com/gagneurlab/splicing-outlier-prediction/tree/master/example/config.yaml) file to `True`.