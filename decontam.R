suppressMessages(suppressWarnings(library('tidyverse')))
suppressMessages(suppressWarnings(library('decontam')))
suppressMessages(suppressWarnings(library('qiime2R')))
suppressMessages(suppressWarnings(library('phyloseq')))

args <- commandArgs(trailingOnly = TRUE)

metadata_file <- args[1]
Decontam_Column <- args[2]

phylo <- qza_to_phyloseq(features = 'feature-table-predecontam.qza', metadata = metadata_file)

contamdf.prev <- decontam::isContaminant(phylo, method = "prevalence", neg = Decontam_Column)

noncontamOTUs <- contamdf.prev %>%
    filter(contaminant == FALSE) %>%
    row.names()

contamdf.prev %>%
    filter(contaminant == TRUE) %>%
    row.names() %>%
    write.csv('contaminant_otus.csv')

phylo_decontam <- prune_taxa(noncontamOTUs, phylo)

phylo_decontam %>%
    otu_table() %>%
    as.matrix() %>%
    biomformat::make_biom() %>%
    biomformat::write_biom('feature-table-decontam.biom')
