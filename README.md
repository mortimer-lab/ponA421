# ponA421
Genomic analysis of the distribution of *ponA* 421 alleles in *Neisseria gonorrhoeae*

## data/

**pcn_res_alleles_metadata.tsv**

Tab-separated table containing metadata, minimum inhibitory concentrations (MICs), genomic quality metrics, and resistance alleles for publicly available genomic data with available penicillin MICs

### data/trees/

**isolates_with_pcn_susceptibility_subsample.final_tree.tre**

Newick file with final phylogeny of 6,082 representative gonococcal genomes after recombination detection with Gubbins

**parsnp.tree**

Newick file with ParSNP generated phylogeny of 10,282 gonococcal genomes

**parsnp.tree_trimmed_list_RTL_0.999**

List of representative sample from ParSNP tree using Treemmer with a relative tree length of 0.999

**parsnp.tree_trimmed_tree_RTL_0.999**

Newick file of ParSNP tree trimmed by Treemmer using a relative tree length of 0.999

#### data/trees/itol_annotations

ITOL annotation files of penicillin MICs and ponA alleles


## scripts

**filter_combine_data.R**

R script that combined original data frames and filtered genomic data based on quality metrics
