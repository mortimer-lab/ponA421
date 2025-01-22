library(tidyverse)
library(glue)

ng_metadata <- read_tsv("../../data/Ng-Combined-Metadata.txt")
ng_metadata <- ng_metadata %>% filter(!reference %in% c("unpublished", "Ryan2018"))

# filter for isolates with PCN MICs

ng_pcn <- ng_metadata %>% filter(!is.na(penicillin))

ng_pcn <- ng_pcn %>% filter(!penicillin %in% c("I", "NoMIC", "R", "S"))

ng_pcn <- ng_pcn %>% mutate(pcn_numeric = parse_number(penicillin))

# filter by sequencing quality metrics

# assembly length
ng_pcn <- ng_pcn %>% filter(assembly_length > 2100000*0.9 & assembly_length < 2100000*1.1)

# contigs
ng_pcn <- ng_pcn %>% filter(contigs < 250)

# coverage
ng_pcn <- ng_pcn %>% filter(reference_coverage > 40)

# percentage mapped
ng_pcn <- ng_pcn %>% filter(reference_percentage_mapped > 80)

# percent ambiguous calls in reference genome
ng_pcn <- ng_pcn %>% filter(percent_missing < 12)

## get resistance alleles

snps <- read_tsv("/n/grad_lab2/Lab/gonococcus/analyses/resistance_alleles/2021-10-07_gc_resistance_alleles.tsv")
blaTEM <- read_tsv("/n/grad_lab2/Lab/gonococcus/analyses/resistance_alleles/2021-10-07_resistance_gene_presence_absence.tsv")
indels <- read_tsv("/n/grad_lab2/Lab/gonococcus/analyses/resistance_alleles/2021-10-07_resistance_inframe_indels.tsv")
penA <- read_tsv("/n/grad_lab2/Lab/gonococcus/analyses/resistance_alleles/2021-10-07_gc_penA_withDescriptions.txt")
porB1a <- read_tsv("/n/grad_lab2/Lab/gonococcus/analyses/resistance_alleles/2021-10-07_porB1a_SNPs.txt")
porB1b <- read_tsv("/n/grad_lab2/Lab/gonococcus/analyses/resistance_alleles/2021-10-07_porB1b_SNPs.txt")
mtr <- read_tsv("/n/grad_lab2/Lab/gonococcus/analyses/resistance_alleles/2021-10-07_mtr_alleles.tsv")


resistance_all <- snps %>% left_join(porB1a) %>% 
        left_join(porB1b) %>%
        left_join(indels) %>%
        left_join(blaTEM) %>%
        left_join(penA) %>%
        left_join(mtr)
 
# references to add: Golparian 2022 (Thailand), Liao 2023 (China), Reimche 2023 (US), Salmeron 2021 (Barcelona)

compile_resistance_alleles <- function(dataset_path){
    snp_file_name <- list.files(path = dataset_path, pattern = "gc_resistance_alleles.tsv$")[[1]]
    snps <- read_tsv(glue("{dataset_path}/{snp_file_name}"), col_types = cols(.default = "c"))
    genes <- read_tsv(glue("{dataset_path}/resistance_gene_presence_absence.tsv"))
    indels <- read_tsv(glue("{dataset_path}/resistance_inframe_indels.tsv"))
    penA <- read_tsv(glue("{dataset_path}/gc_penA.txt"))
    porB1a <- read_tsv(glue("{dataset_path}/porB1a_SNPs.txt"))
    porB1b <- read_tsv(glue("{dataset_path}/porB1b_SNPs.txt"))
    rRNA <- read_tsv(glue("{dataset_path}/rRNA_allele_summary.txt"))
    mtr <- read_tsv(glue("{dataset_path}/mtr_alleles.tsv"))
    resistance <- snps %>% left_join(porB1a) %>% 
        left_join(porB1b) %>%
        left_join(indels) %>%
        left_join(genes) %>%
        left_join(penA) %>%
        left_join(mtr) %>%
        left_join(rRNA)
    resistance
}

res_alleles_liao <- compile_resistance_alleles("/n/grad_lab2/Lab/gonococcus/datasets/liao_2023_china_penA60/resistance") 
res_alleles_golparian <- compile_resistance_alleles("/n/grad_lab2/Lab/gonococcus/datasets/golparian_2022_thailand/resistance")
res_alleles_reimche <- compile_resistance_alleles("/n/grad_lab2/Lab/gonococcus/datasets/reimche_2023_US/resistance")
res_alleles_salmeron <- compile_resistance_alleles("/n/grad_lab2/Lab/gonococcus/datasets/salmeron_2021_barcelona/resistance")

res_alleles_new <- res_alleles_liao %>% rows_insert(res_alleles_golparian) %>%
    rows_insert(res_alleles_reimche) %>%
    rows_insert(res_alleles_salmeron)

res_alleles_new <- res_alleles_new %>% select(-GyrB_467, -starts_with("23S"), -starts_with("16S"))
resistance_all <- resistance_all %>% rows_insert(res_alleles_new)

# combine with metadata

ng_pcn <- ng_pcn %>% left_join(resistance_all)

# remove isolates with uncertain resistance alleles at relevant loci

ng_pcn %>% select(starts_with("mtr"), starts_with("porB"), starts_with("PBP"), starts_with("penA")) %>% summarise_all(list(~sum(is.na(.)))) %>% pivot_longer(everything(), names_to = "locus", values_to = "num_missing") %>% print(n=21)

ng_pcn <- ng_pcn %>% filter(!is.na(PBP1_421))
ng_pcn <- ng_pcn %>% filter(!(mtrR == "GC_allele" & (is.na(MtrR_45) | is.na(MtrR_39))))
ng_pcn <- ng_pcn %>% filter(penA != "unknown")

ng_pcn %>% write_tsv("../data/pcn_res_alleles_metadata.tsv")

# combined all metadata will resistance alleles to include isolates without penicillin MICs
ponA_421_dates <- ng_metadata %>% left_join(resistance_all) %>% filter(!is.na(PBP1_421)) %>% separate(date, c("Year", "Month", "Day", "-")) %>% select(wgs_id, Year, PBP1_421) %>% filter(!Year %in% c("01", "02", "03", "04", "05", "06"))
ponA_421_dates %>% count(Year, PBP1_421) %>% filter(!is.na(Year)) %>% write_tsv("../data/ponA_alleles_by_year.tsv")
