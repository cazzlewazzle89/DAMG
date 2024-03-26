library(rstatix)
library(phyloseq)
library(tidyverse)
library(qiime2R)
library(mixOmics)
library(PLSDAbatch)

# import data and build phyloseq object
featuretable <- read_qza('feature-table.qza')$data

taxonomy <- read_qza('taxonomy.qza')$data %>%
  select(-Confidence) %>%
  separate(Taxon,
           into = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
           sep = '; ', remove = F) %>%
  column_to_rownames('Feature.ID')

metadata <- bind_rows(read_tsv('metadata_run1.tsv'),
                      read_tsv('metadata_run2.tsv')) %>%
  column_to_rownames('Sample')
saveRDS(metadata, 'metadata.RDS')

tree <- ape::read.tree('tree.nwk')

phylo <- phyloseq(otu_table(featuretable, taxa_are_rows = T),
                  tax_table(as.matrix(taxonomy)),
                  tree,
                  sample_data(metadata))

saveRDS(phylo, 'phylo_raw.Rds')

rm(featuretable, taxonomy, tree)

metadata <- metadata %>%
  rownames_to_column('Sample')

# filtering featuretable
keep <- otu_table(phylo) %>%
  data.frame() %>%
  rownames_to_column('OTU') %>%
  pivot_longer(cols = -OTU, names_to = 'Sample', values_to = 'Count') %>%
  group_by(Sample) %>%
  mutate(RA = 100*(Count/sum(Count))) %>%
  ungroup() %>%
  left_join(metadata) %>%
  group_by(OTU, Source) %>%
  summarise(MeanRA = mean(RA)) %>%
  filter(MeanRA >= 0.1) %>%
  pull(OTU) %>%
  unique()

phylo_filt <- prune_taxa(keep, phylo)

saveRDS(phylo_filt, 'phylo_filt.RDS')

# looking at features per sample
phylo <- readRDS('phylo_filt.RDS')

p <- sample_sums(phylo) %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  rename(FeaturesPerSample = 2) %>%
  left_join(metadata) %>%
  ggplot(aes(x = Genotype, y = FeaturesPerSample)) +
  facet_wrap(~Experiment, nrow = 1) +
  geom_boxplot(aes(colour = Source), fill = NA, outlier.shape = NA, show.legend = F) +
  geom_point(aes(fill = Source, shape = Genotype), 
             colour = 'black', size = 3, position = position_jitterdodge(jitter.width = 0.2)) +
  theme_bw() +
  scale_y_log10() +
  scale_shape_manual(values = c(21, 23, 24)) +
  guides(fill = guide_legend(override.aes = list(shape = 22))) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = '')

pdf('featurespersample.pdf', width = 6, height = 3)
print(p)
dev.off()

featurespersource <- sample_sums(phylo) %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  rename(FeaturesPerSample = 2) %>%
  left_join(metadata) %>%
  group_by(Source) %>%
  summarise(Mean = mean(FeaturesPerSample),
            Median = median(FeaturesPerSample),
            SD = sd(FeaturesPerSample))

featurespergenotypesource <- sample_sums(phylo) %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  rename(FeaturesPerSample = 2) %>%
  left_join(metadata) %>%
  group_by(Source, Genotype) %>%
  summarise(Mean = mean(FeaturesPerSample),
            Median = median(FeaturesPerSample),
            SD = sd(FeaturesPerSample))

list(Source = featurespersource,
     SourceGenotype = featurespergenotypesource) %>%
  writexl::write_xlsx('featurespersample_summary.xlsx')

# rarefying to lowest sample depth for alpha diversity analysis
# doing this overall and per-sampletype
phylo_ear <- prune_samples(sample_data(phylo)$Source == 'Ear', phylo)
phylo_faecal <- prune_samples(sample_data(phylo)$Source == 'Faecal', phylo)
phylo_lung <- prune_samples(sample_data(phylo)$Source == 'Lung', phylo)

raredepth_all <- phylo %>% sample_sums() %>% min()
raredepth_ear <- phylo_ear %>% sample_sums() %>% min()
raredepth_faecal <- phylo_faecal %>% sample_sums() %>% min()
raredepth_lung <- phylo_lung %>% sample_sums() %>% min()

phylo_rare_all <- rarefy_even_depth(phylo, raredepth_all, rngseed = 42)
phylo_ear_rare <- rarefy_even_depth(phylo_ear, raredepth_ear, rngseed = 42)
phylo_faecal_rare <- rarefy_even_depth(phylo_faecal, raredepth_faecal, rngseed = 42)
phylo_lung_rare <- rarefy_even_depth(phylo_lung, raredepth_lung, rngseed = 42)

saveRDS(phylo_rare_all, 'phylo_rare_all.Rds')
saveRDS(phylo_ear_rare, 'phylo_rare_ear.Rds')
saveRDS(phylo_faecal_rare, 'phylo_rare_faecal.Rds')
saveRDS(phylo_lung_rare, 'phylo_rare_lung.Rds')

# calculating alpha diversity
richness_all <- estimate_richness(phylo_rare_all)
pd_all <- btools::estimate_pd(phylo_rare_all)
alpha_all <- bind_cols(richness_all, pd_all) %>% 
  dplyr::select(-SR) %>% 
  rownames_to_column('Sample') %>% 
  left_join(metadata)

richness_ear <- estimate_richness(phylo_ear_rare)
pd_ear <- btools::estimate_pd(phylo_ear_rare)
alpha_ear <- bind_cols(richness_ear, pd_ear) %>% 
  dplyr::select(-SR) %>% 
  rownames_to_column('Sample') %>% 
  left_join(metadata)

richness_faecal <- estimate_richness(phylo_faecal_rare)
pd_faecal <- btools::estimate_pd(phylo_faecal_rare)
alpha_faecal <- bind_cols(richness_faecal, pd_faecal) %>%
  dplyr::select(-SR) %>%
  rownames_to_column('Sample') %>% 
  left_join(metadata)

richness_lung <- estimate_richness(phylo_lung_rare)
pd_lung <- btools::estimate_pd(phylo_lung_rare)
alpha_lung <- bind_cols(richness_lung, pd_lung) %>% 
  dplyr::select(-SR) %>% 
  rownames_to_column('Sample') %>% 
  left_join(metadata)

list(All = alpha_all,
     Ear = alpha_ear,
     Faecal = alpha_faecal,
     Lung = alpha_lung) %>%
  writexl::write_xlsx('alpha.xlsx')

# alpha diversity comparisons
alpha_all <- readxl::read_xlsx('alpha.xlsx', sheet = 'All')
alpha_ear <- readxl::read_xlsx('alpha.xlsx', sheet = 'Ear')
alpha_faecal <- readxl::read_xlsx('alpha.xlsx', sheet = 'Faecal')
alpha_lung <- readxl::read_xlsx('alpha.xlsx', sheet = 'Lung')

alpha_all %>%
  filter(Source == 'Ear') %>%
  anova_test(Shannon ~ Experiment + Genotype)
alpha_all %>%
  filter(Source == 'Ear') %>%
  anova_test(Shannon ~ Genotype)

alpha_all %>%
  filter(Source == 'Faecal') %>%
  anova_test(Shannon ~ Experiment + Genotype)
alpha_all %>%
  filter(Source == 'Faecal') %>%
  anova_test(Shannon ~ Genotype)

alpha_all %>%
  filter(Source == 'Lung') %>%
  anova_test(Shannon ~ Experiment + Genotype)
alpha_all %>%
  filter(Source == 'Lung') %>%
  anova_test(Shannon ~ Genotype)

list('Shannon_Ear' = alpha_ear %>% anova_test(Shannon ~ Experiment * Genotype) %>% data.frame(check.names = F),
     'Shannon_Faecal' = alpha_faecal %>% anova_test(Shannon ~ Experiment * Genotype) %>% data.frame(check.names = F),
     'Shannon_Lung' = alpha_lung %>% anova_test(Shannon ~ Experiment * Genotype) %>% data.frame(check.names = F),
     'PD_Ear' = alpha_ear %>% anova_test(PD ~ Experiment * Genotype) %>% data.frame(check.names = F),
     'PD_Faecal' = alpha_faecal %>% anova_test(PD ~ Experiment * Genotype) %>% data.frame(check.names = F),
     'PD_Lung' = alpha_lung %>% anova_test(PD ~ Experiment * Genotype) %>% data.frame(check.names = F)) %>%
  writexl::write_xlsx('alpha_anova.xlsx')

p <- bind_rows(alpha_ear %>% dplyr::select(Sample, Genotype, Source, Experiment, Shannon, PD),
          alpha_faecal %>% dplyr::select(Sample, Genotype, Source, Experiment,  Shannon, PD),
          alpha_lung %>% dplyr::select(Sample, Genotype, Source, Experiment,  Shannon, PD)) %>%
  pivot_longer(cols = Shannon:PD, names_to = 'Measure', values_to = 'Diversity') %>%
  mutate(Measure = str_replace(Measure, 'PD', 'Phylogenetic Diversity')) %>%
  mutate(Measure = str_replace(Measure, 'Shannon', 'Shannon Entropy')) %>%
  ggplot(aes(x = 1, y = Diversity)) +
  geom_boxplot(aes(colour = Genotype), fill = NA, outlier.shape = NA, show.legend = F) +
  geom_point(aes(fill = Genotype, shape = Experiment), colour = 'black', size = 3,
             position = position_jitterdodge(jitter.width = 0.2)) +
  facet_grid(Measure~Source, scales = 'free_y') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = '') +
  scale_shape_manual(values = c(21, 23)) +
  guides(fill = guide_legend(override.aes = list(shape = 22)))

pdf('alpha.pdf', width = 8, height = 4)
print(p)
dev.off()

# batch correcting data
# trying both PLSDA-batch and Limma
library(mixOmics)
library(Biobase)
library(PLSDAbatch)
featuretable_clr <- logratio.transfo(X = t(otu_table(phylo)), logratio = 'CLR', offset = 1)
sourcegenotype <- paste(sample_data(phylo)$Source, sample_data(phylo)$Genotype, sep = '')

pca_before <- pca(featuretable_clr, ncomp = 3, scale = TRUE)
Scatter_Density(object = pca_before, batch = sample_data(phylo)$Experiment, trt = sample_data(phylo)$Source, 
                batch.legend.title = 'Experiment', trt.legend.title = 'Source')

treatment_tune <- plsda(X = featuretable_clr, Y = sourcegenotype, ncomp = 5)
treatment_tune$prop_expl_var

batch_tune <- PLSDA_batch(X = featuretable_clr, 
                          Y.trt = sourcegenotype, Y.bat = sample_data(phylo)$Experiment,
                          ncomp.trt = 3, ncomp.bat = 5)
batch_tune$explained_variance.bat 

batch_res <- PLSDA_batch(X = featuretable_clr, 
                         Y.trt = sourcegenotype, Y.bat = sample_data(phylo)$Experiment,
                         ncomp.trt = 3, ncomp.bat = 1)
batch <- batch_res$X.nobatch

pca_after <- pca(batch, ncomp = 3, scale = TRUE)
Scatter_Density(object = pca_after, batch = sample_data(phylo)$Experiment, trt = sample_data(phylo)$Source, 
                batch.legend.title = 'Cohort', trt.legend.title = 'Treatment')

library(limma)
clr_nobatch <- removeBatchEffect(t(featuretable_clr), batch = sample_data(phylo)$Experiment, model = model.matrix(~ Source + Genotype, data = sample_data(phylo)))
clr_nobatch <- t(clr_nobatch)

pca_before <- pca(featuretable_clr, ncomp = 3, scale = TRUE)
pca_after <- pca(clr_nobatch, ncomp = 3, scale = TRUE)

Scatter_Density(object = pca_before, batch = sample_data(phylo)$Experiment, trt = sample_data(phylo)$Source, 
                batch.legend.title = 'Cohort', trt.legend.title = 'Treatment')
Scatter_Density(object = pca_after, batch = sample_data(phylo)$Experiment, trt = sample_data(phylo)$Source, 
                batch.legend.title = 'Cohort', trt.legend.title = 'Treatment')

saveRDS(clr_nobatch, 'featuretable_clr_limmacorrected.RDS')

# limma seemed to work the best so going to use that for the single bodysite objects
clr_nobatch_ear <- removeBatchEffect(t(featuretable_clr_ear), batch = sample_data(phylo_ear)$Experiment, model = model.matrix(~ Source + Genotype, data = sample_data(phylo_ear)))
clr_nobatch_ear <- t(clr_nobatch_ear)
pca_before <- pca(featuretable_clr_ear, ncomp = 3, scale = TRUE)
pca_after <- pca(clr_nobatch_ear, ncomp = 3, scale = TRUE)
Scatter_Density(object = pca_before, batch = sample_data(phylo_ear)$Experiment, trt = sample_data(phylo_ear)$Source, 
                batch.legend.title = 'Cohort', trt.legend.title = 'Treatment')
Scatter_Density(object = pca_after, batch = sample_data(phylo_ear)$Experiment, trt = sample_data(phylo_ear)$Source, 
                batch.legend.title = 'Cohort', trt.legend.title = 'Treatment')
saveRDS(clr_nobatch_ear, 'featuretable_clr_limmacorrected_ear.RDS')

clr_nobatch_faecal <- removeBatchEffect(t(featuretable_clr_faecal), batch = sample_data(phylo_faecal)$Experiment, model = model.matrix(~ Source + Genotype, data = sample_data(phylo_faecal)))
clr_nobatch_faecal <- t(clr_nobatch_faecal)
pca_before <- pca(featuretable_clr_faecal, ncomp = 3, scale = TRUE)
pca_after <- pca(clr_nobatch_faecal, ncomp = 3, scale = TRUE)
Scatter_Density(object = pca_before, batch = sample_data(phylo_faecal)$Experiment, trt = sample_data(phylo_faecal)$Source, 
                batch.legend.title = 'Cohort', trt.legend.title = 'Treatment')
Scatter_Density(object = pca_after, batch = sample_data(phylo_faecal)$Experiment, trt = sample_data(phylo_faecal)$Source, 
                batch.legend.title = 'Cohort', trt.legend.title = 'Treatment')
saveRDS(clr_nobatch_faecal, 'featuretable_clr_limmacorrected_faecal.RDS')

clr_nobatch_lung <- removeBatchEffect(t(featuretable_clr_lung), batch = sample_data(phylo_lung)$Experiment, model = model.matrix(~ Source + Genotype, data = sample_data(phylo_lung)))
clr_nobatch_lung <- t(clr_nobatch_lung)
pca_before <- pca(featuretable_clr_lung, ncomp = 3, scale = TRUE)
pca_after <- pca(clr_nobatch_lung, ncomp = 3, scale = TRUE)
Scatter_Density(object = pca_before, batch = sample_data(phylo_lung)$Experiment, trt = sample_data(phylo_lung)$Source, 
                batch.legend.title = 'Cohort', trt.legend.title = 'Treatment')
Scatter_Density(object = pca_after, batch = sample_data(phylo_lung)$Experiment, trt = sample_data(phylo_lung)$Source, 
                batch.legend.title = 'Cohort', trt.legend.title = 'Treatment')
saveRDS(clr_nobatch_lung, 'featuretable_clr_limmacorrected_lung.RDS')

# calculating beta diversity using aitchison distance
library(vegan)

dist_clr <- vegdist(clr_nobatch, method = 'euclidean')
nmds_clr <- metaMDS(dist_clr, k = 2)
points_clr <- nmds_clr$points %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(metadata)

p <- points_clr %>%
  ggplot(aes(x = MDS1, y = MDS2)) +
  geom_point(aes(fill = Source, shape = Genotype), colour = 'black', size = 3) +
  theme_bw() +
  scale_shape_manual(values = c(21, 23, 24)) +
  guides(fill = guide_legend(override.aes = list(shape = 22))) 

pdf('beta_all.pdf', width = 6, height = 4)
p
dev.off()

dist_clr_ear <- vegdist(clr_nobatch_ear, method = 'euclidean')
nmds_clr_ear <- metaMDS(dist_clr_ear, k = 2)
points_clr_ear <- nmds_clr_ear$points %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(metadata)

dist_clr_faecal <- vegdist(clr_nobatch_faecal, method = 'euclidean')
nmds_clr_faecal <- metaMDS(dist_clr_faecal, k = 2)
points_clr_faecal <- nmds_clr_faecal$points %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(metadata)

dist_clr_lung <- vegdist(clr_nobatch_lung, method = 'euclidean')
nmds_clr_lung <- metaMDS(dist_clr_lung, k = 2)
points_clr_lung <- nmds_clr_lung$points %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(metadata)

points_clr_combined <- bind_rows(points_clr_ear,
                                 points_clr_faecal,
                                 points_clr_lung)

p <- points_clr_combined %>%
  ggplot(aes(x = MDS1, y = MDS2)) +
  facet_wrap(~Source, nrow = 1, scales = 'free') +
  geom_point(aes(fill = Genotype, shape = Experiment), colour = 'black', size = 3) +
  stat_ellipse(aes(colour = Genotype), show.legend = F) +
  theme_bw() +
  scale_shape_manual(values = c(21, 23)) +
  guides(fill = guide_legend(override.aes = list(shape = 22)))

pdf('beta_faceted.pdf', width = 10, height = 3)
print(p)
dev.off()

# permanova on rclr distances
dist_clr <- usedist::dist_subset(dist_clr, sample_names(phylo))
permout_all <- vegan::adonis2(dist_clr ~ Source * Genotype, data = sample_data(phylo) %>%
                                data.frame() %>% 
                                rownames_to_column('Sample'))

dist_clr_ear <- usedist::dist_subset(dist_clr_ear, sample_names(phylo_ear))
permout_ear <- vegan::adonis2(dist_clr_ear ~ Genotype, data = sample_data(phylo_ear) %>%
                                data.frame() %>% 
                                rownames_to_column('Sample'))

dist_clr_faecal <- usedist::dist_subset(dist_clr_faecal, sample_names(phylo_faecal))
permout_faecal <- vegan::adonis2(dist_clr_faecal ~ Genotype, data = sample_data(phylo_faecal) %>%
                                   data.frame() %>% 
                                   rownames_to_column('Sample'))

dist_clr_lung <- usedist::dist_subset(dist_clr_lung, sample_names(phylo_lung))
permout_lung <- vegan::adonis2(dist_clr_lung ~ Genotype, data = sample_data(phylo_lung) %>%
                                 data.frame() %>% 
                                 rownames_to_column('Sample'))

list(All = permout_all %>% data.frame(check.names = F) %>% rownames_to_column('Variable'),
     Ear = permout_ear %>% data.frame(check.names = F) %>% rownames_to_column('Variable'),
     Faecal = permout_faecal %>% data.frame(check.names = F) %>% rownames_to_column('Variable'),
     Lung = permout_lung %>% data.frame(check.names = F) %>% rownames_to_column('Variable')) %>%
  writexl::write_xlsx('permout.xlsx')

# comparing betadisper between groups
betadisper_ear <- betadisper(dist_clr_ear, sample_data(phylo_ear)$Genotype, type = 'centroid')
betadisper_faecal <- betadisper(dist_clr_faecal, sample_data(phylo_faecal)$Genotype, type = 'centroid')
betadisper_lung <- betadisper(dist_clr_lung, sample_data(phylo_lung)$Genotype, type = 'centroid')

list(anova(betadisper_ear) %>% data.frame(check.names = F) %>% rownames_to_column('Variable'),
     anova(betadisper_faecal) %>% data.frame(check.names = F) %>% rownames_to_column('Variable'),
     anova(betadisper_lung) %>% data.frame(check.names = F) %>% rownames_to_column('Variable')) %>%
  writexl::write_xlsx('betadisper.xlsx')

# PLSDA
library(mixOmics)

# EAR SAMPLES
in_mat <- phylo_ear %>%
  otu_table() %>%
  data.frame() %>%
  t()

in_class <- phylo_ear %>%
  sample_data() %>%
  pull(Genotype)

plsda_ear <- plsda(in_mat, in_class, ncomp = 10)

background_ear <- background.predict(plsda_ear,
                                     comp.predicted = 2,
                                     dist = 'max.dist')

# FAECAL SAMPLES
in_mat <- phylo_faecal %>%
  otu_table() %>%
  data.frame() %>%
  t()

in_class <- phylo_faecal %>%
  sample_data() %>%
  pull(Genotype)

plsda_faecal <- plsda(in_mat, in_class, ncomp = 10)

background_faecal <- background.predict(plsda_faecal,
                                        comp.predicted = 2,
                                        dist = 'max.dist')

# LUNG SAMPLES
in_mat <- phylo_lung %>%
  otu_table() %>%
  data.frame() %>%
  t()

in_class <- phylo_lung %>%
  sample_data() %>%
  pull(Genotype)

plsda_lung <- plsda(in_mat, in_class, ncomp = 10)

background_lung <- background.predict(plsda_lung,
                                      comp.predicted = 2,
                                      dist = 'max.dist')


# COMBINING PLOTS
pdf('plsda_ear.pdf', width = 5, height = 4)
plotIndiv(plsda_ear, comp = c(1,2),
          group = in_class, ind.names = FALSE,
          ellipse = TRUE, legend = TRUE,
          title = 'Ear',
          background = background_ear)
dev.off()

pdf('plsda_faecal.pdf', width = 5, height = 4)
plotIndiv(plsda_faecal, comp = c(1,2),
          group = in_class, ind.names = FALSE,
          ellipse = TRUE, legend = TRUE,
          title = 'Faecal',
          background = background_faecal)
dev.off()

pdf('plsda_lung.pdf', width = 5, height = 4)
plotIndiv(plsda_lung, comp = c(1,2),
          group = in_class, ind.names = FALSE,
          ellipse = TRUE, legend = TRUE,
          title = 'Lung',
          background = background_lung)
dev.off()

# plotting boxplots of X and Y axes of batch-corrected PLS-DA of lung samples
df <- data.frame(X = plsda_lung$variates$X[,1],
                 Y = plsda_lung$variates$X[,2]) %>%
  rownames_to_column('Sample') %>%
  left_join(metadata) %>%
  rename(Xvar = X,
         Yvar = Y)

df %>%
  anova_test(Xvar ~ Genotype)

df %>%
  pairwise_t_test(Xvar ~ Genotype)

df %>%
  anova_test(Yvar ~ Genotype)

df %>%
  pairwise_t_test(Yvar ~ Genotype)

genotype_colours <- c('CKO' = 'blue', 'KO' = 'orange', 'WT' = 'grey')

p_scatter <- df %>%
  ggplot(aes(x = Xvar, y = Yvar)) +
  geom_point(aes(fill = Genotype), shape = 21, colour = 'black', size = 3, show.legend = F) +
  stat_ellipse(aes(colour = Genotype), show.legend = F) +
  theme_bw() +
  labs(x = paste0('X-variate 1: ', round(100 * plsda_lung$prop_expl_var$X[1], 0), '% expl. var'),
       y = paste0('X-variate 2: ', round(100 * plsda_lung$prop_expl_var$X[2], 0), '% expl. var')) +
  scale_fill_manual(values = genotype_colours) +
  scale_colour_manual(values = genotype_colours)

p_box_x <- df %>%
  ggplot() +
  geom_boxplot(aes(x = Xvar, y = 1, fill = Genotype), colour = 'black', outlier.shape = NA, show.legend = F) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank()) +
  labs(x = '',
       y = '') +
  scale_fill_manual(values = genotype_colours)

p_box_y <- df %>%
  ggplot() +
  geom_boxplot(aes(y = Yvar, x = 1, fill = Genotype), colour = 'black', outlier.shape = NA, show.legend = F) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank()) +
  labs(x = '',
       y = '') +
  scale_fill_manual(values = genotype_colours)

pdf('plsda_lung_withboxplots.pdf', width = 5, height = 4)
gridExtra::grid.arrange(p_box_x, NULL, p_scatter, p_box_y,
                        ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
dev.off()

# plotting boxplots of X and Y axes of batch-corrected PLS-DA of faecal samples
df <- data.frame(X = plsda_faecal$variates$X[,1],
                 Y = plsda_faecal$variates$X[,2]) %>%
  rownames_to_column('Sample') %>%
  left_join(metadata) %>%
  rename(Xvar = X,
         Yvar = Y)

df %>%
  anova_test(Xvar ~ Genotype)

df %>%
  pairwise_t_test(Xvar ~ Genotype)

df %>%
  anova_test(Yvar ~ Genotype)

df %>%
  pairwise_t_test(Yvar ~ Genotype)

genotype_colours <- c('CKO' = 'blue', 'KO' = 'orange', 'WT' = 'grey')

p_scatter <- df %>%
  ggplot(aes(x = Xvar, y = Yvar)) +
  geom_point(aes(fill = Genotype), shape = 21, colour = 'black', size = 3, show.legend = F) +
  stat_ellipse(aes(colour = Genotype), show.legend = F) +
  theme_bw() +
  labs(x = paste0('X-variate 1: ', round(100 * plsda_faecal$prop_expl_var$X[1], 0), '% expl. var'),
       y = paste0('X-variate 2: ', round(100 * plsda_faecal$prop_expl_var$X[2], 0), '% expl. var')) +
  scale_fill_manual(values = genotype_colours) +
  scale_colour_manual(values = genotype_colours)

p_box_x <- df %>%
  ggplot() +
  geom_boxplot(aes(x = Xvar, y = 1, fill = Genotype), colour = 'black', outlier.shape = NA, show.legend = F) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank()) +
  labs(x = '',
       y = '') +
  scale_fill_manual(values = genotype_colours)

p_box_y <- df %>%
  ggplot() +
  geom_boxplot(aes(y = Yvar, x = 1, fill = Genotype), colour = 'black', outlier.shape = NA, show.legend = F) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank()) +
  labs(x = '',
       y = '') +
  scale_fill_manual(values = genotype_colours)

pdf('plsda_faecal_withboxplots.pdf', width = 5, height = 4)
gridExtra::grid.arrange(p_box_x, NULL, p_scatter, p_box_y,
                        ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
dev.off()

# plotting boxplots of X and Y axes of batch-corrected PLS-DA of ear samples
df <- data.frame(X = plsda_ear$variates$X[,1],
                 Y = plsda_ear$variates$X[,2]) %>%
  rownames_to_column('Sample') %>%
  left_join(metadata) %>%
  rename(Xvar = X,
         Yvar = Y)

df %>%
  anova_test(Xvar ~ Genotype)

df %>%
  pairwise_t_test(Xvar ~ Genotype)

df %>%
  anova_test(Yvar ~ Genotype)

df %>%
  pairwise_t_test(Yvar ~ Genotype)


genotype_colours <- c('CKO' = 'blue', 'KO' = 'orange', 'WT' = 'grey')

p_scatter <- df %>%
  ggplot(aes(x = Xvar, y = Yvar)) +
  geom_point(aes(fill = Genotype), shape = 21, colour = 'black', size = 3, show.legend = F) +
  stat_ellipse(aes(colour = Genotype), show.legend = F) +
  theme_bw() +
  labs(x = paste0('X-variate 1: ', round(100 * plsda_ear$prop_expl_var$X[1], 0), '% expl. var'),
       y = paste0('X-variate 2: ', round(100 * plsda_ear$prop_expl_var$X[2], 0), '% expl. var')) +
  scale_fill_manual(values = genotype_colours) +
  scale_colour_manual(values = genotype_colours)

p_box_x <- df %>%
  ggplot() +
  geom_boxplot(aes(x = Xvar, y = 1, fill = Genotype), colour = 'black', outlier.shape = NA, show.legend = F) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank()) +
  labs(x = '',
       y = '') +
  scale_fill_manual(values = genotype_colours)

p_box_y <- df %>%
  ggplot() +
  geom_boxplot(aes(y = Yvar, x = 1, fill = Genotype), colour = 'black', outlier.shape = NA, show.legend = F) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank()) +
  labs(x = '',
       y = '') +
  scale_fill_manual(values = genotype_colours)

pdf('plsda_ear_withboxplots.pdf', width = 5, height = 4)
gridExtra::grid.arrange(p_box_x, NULL, p_scatter, p_box_y,
                        ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
dev.off()

# same plots as above but with sample labels
df <- data.frame(X = plsda_lung$variates$X[,1],
                 Y = plsda_lung$variates$X[,2]) %>%
  rownames_to_column('Sample') %>%
  left_join(metadata) %>%
  rename(Xvar = X,
         Yvar = Y)

genotype_colours <- c('CKO' = 'blue', 'KO' = 'orange', 'WT' = 'grey')

p_scatter <- df %>%
  ggplot(aes(x = Xvar, y = Yvar)) +
  geom_point(aes(fill = Genotype), shape = 21, colour = 'black', size = 3, show.legend = F) +
  geom_text(aes(label = MouseNo), nudge_x = 0.25, nudge_y = 0.25) +
  stat_ellipse(aes(colour = Genotype), show.legend = F) +
  theme_bw() +
  labs(x = paste0('X-variate 1: ', round(100 * plsda_lung$prop_expl_var$X[1], 0), '% expl. var'),
       y = paste0('X-variate 2: ', round(100 * plsda_lung$prop_expl_var$X[2], 0), '% expl. var')) +
  scale_fill_manual(values = genotype_colours) +
  scale_colour_manual(values = genotype_colours)

p_box_x <- df %>%
  ggplot() +
  geom_boxplot(aes(x = Xvar, y = 1, fill = Genotype), colour = 'black', outlier.shape = NA, show.legend = F) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank()) +
  labs(x = '',
       y = '') +
  scale_fill_manual(values = genotype_colours)

p_box_y <- df %>%
  ggplot() +
  geom_boxplot(aes(y = Yvar, x = 1, fill = Genotype), colour = 'black', outlier.shape = NA, show.legend = F) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank()) +
  labs(x = '',
       y = '') +
  scale_fill_manual(values = genotype_colours)

pdf('plsda_lung_withboxplots_labelled.pdf', width = 5, height = 4)
gridExtra::grid.arrange(p_box_x, NULL, p_scatter, p_box_y,
                        ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
dev.off()

# plotting boxplots of X and Y axes of batch-corrected PLS-DA of faecal samples
df <- data.frame(X = plsda_faecal$variates$X[,1],
                 Y = plsda_faecal$variates$X[,2]) %>%
  rownames_to_column('Sample') %>%
  left_join(metadata) %>%
  rename(Xvar = X,
         Yvar = Y)

genotype_colours <- c('CKO' = 'blue', 'KO' = 'orange', 'WT' = 'grey')

p_scatter <- df %>%
  ggplot(aes(x = Xvar, y = Yvar)) +
  geom_point(aes(fill = Genotype), shape = 21, colour = 'black', size = 3, show.legend = F) +
  geom_text(aes(label = MouseNo), nudge_x = 0.25, nudge_y = 0.25) +
  stat_ellipse(aes(colour = Genotype), show.legend = F) +
  theme_bw() +
  labs(x = paste0('X-variate 1: ', round(100 * plsda_faecal$prop_expl_var$X[1], 0), '% expl. var'),
       y = paste0('X-variate 2: ', round(100 * plsda_faecal$prop_expl_var$X[2], 0), '% expl. var')) +
  scale_fill_manual(values = genotype_colours) +
  scale_colour_manual(values = genotype_colours)

p_box_x <- df %>%
  ggplot() +
  geom_boxplot(aes(x = Xvar, y = 1, fill = Genotype), colour = 'black', outlier.shape = NA, show.legend = F) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank()) +
  labs(x = '',
       y = '') +
  scale_fill_manual(values = genotype_colours)

p_box_y <- df %>%
  ggplot() +
  geom_boxplot(aes(y = Yvar, x = 1, fill = Genotype), colour = 'black', outlier.shape = NA, show.legend = F) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank()) +
  labs(x = '',
       y = '') +
  scale_fill_manual(values = genotype_colours)

pdf('plsda_faecal_withboxplots_labelled.pdf', width = 5, height = 4)
gridExtra::grid.arrange(p_box_x, NULL, p_scatter, p_box_y,
                        ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
dev.off()

# plotting boxplots of X and Y axes of batch-corrected PLS-DA of ear samples
df <- data.frame(X = plsda_ear$variates$X[,1],
                 Y = plsda_ear$variates$X[,2]) %>%
  rownames_to_column('Sample') %>%
  left_join(metadata) %>%
  rename(Xvar = X,
         Yvar = Y)

genotype_colours <- c('CKO' = 'blue', 'KO' = 'orange', 'WT' = 'grey')

p_scatter <- df %>%
  ggplot(aes(x = Xvar, y = Yvar)) +
  geom_point(aes(fill = Genotype), shape = 21, colour = 'black', size = 3, show.legend = F) +
  geom_text(aes(label = MouseNo), nudge_x = 0.25, nudge_y = 0.25) +
  stat_ellipse(aes(colour = Genotype), show.legend = F) +
  theme_bw() +
  labs(x = paste0('X-variate 1: ', round(100 * plsda_ear$prop_expl_var$X[1], 0), '% expl. var'),
       y = paste0('X-variate 2: ', round(100 * plsda_ear$prop_expl_var$X[2], 0), '% expl. var')) +
  scale_fill_manual(values = genotype_colours) +
  scale_colour_manual(values = genotype_colours)

p_box_x <- df %>%
  ggplot() +
  geom_boxplot(aes(x = Xvar, y = 1, fill = Genotype), colour = 'black', outlier.shape = NA, show.legend = F) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank()) +
  labs(x = '',
       y = '') +
  scale_fill_manual(values = genotype_colours)

p_box_y <- df %>%
  ggplot() +
  geom_boxplot(aes(y = Yvar, x = 1, fill = Genotype), colour = 'black', outlier.shape = NA, show.legend = F) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank()) +
  labs(x = '',
       y = '') +
  scale_fill_manual(values = genotype_colours)

pdf('plsda_ear_withboxplots_labelled.pdf', width = 5, height = 4)
gridExtra::grid.arrange(p_box_x, NULL, p_scatter, p_box_y,
                        ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
dev.off()

# comparing within-mouse distances between different sample sources
dist_clr <- usedist::dist_subset(d = dist_clr, idx = metadata$Sample)

df <- dist_clr %>%
  as.matrix() %>%
  data.frame() %>%
  rownames_to_column('Sample1') %>%
  pivot_longer(cols = -Sample1, names_to = 'Sample2', values_to = 'Dist') %>%
  filter(Sample1 < Sample2) %>%
  left_join(metadata, by = c('Sample1' = 'Sample')) %>%
  left_join(metadata, by = c('Sample2' = 'Sample')) %>%
  filter(MouseNo.x == MouseNo.y) %>%
  mutate(Comparison = paste0(Source.x, '_', Source.y))

df %>%
  anova_test(Dist ~ Comparison * Genotype.x)
df %>%
  pairwise_t_test(Dist ~ Comparison)

df %>% filter(Comparison == 'Ear_Faecal') %>%  anova_test(Dist ~ Genotype.x)
df %>% filter(Comparison == 'Ear_Lung') %>% anova_test(Dist ~ Genotype.x)
df %>% filter(Comparison == 'Faecal_Lung') %>%anova_test(Dist ~ Genotype.x)  

p <- df %>%
  ggplot(aes(x = Comparison, y = Dist)) +
  geom_boxplot(aes(colour = Genotype.x), fill = NA, outlier.shape = NA, show.legend = F) +
  geom_point(aes(fill = Genotype.x), 
             colour = 'black', size = 3, shape = 21, position = position_jitterdodge(jitter.width = 0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = '')
pdf('pairwisedist.pdf', width = 6, height = 4)
p
dev.off()

df %>% 
  dplyr::select(Sample1, Sample2, Dist, Comparison, Genotype.x, Mouse.x, MouseNo.x, Experiment.x) %>%
  rename(Genotype = Genotype.x,
         Mouse = Mouse.x,
         MouseNo = MouseNo.x,
         Experiment = Experiment.x) %>%
  arrange(Comparison, Genotype, Dist) %>%
  write_tsv('pairwisedist.tsv')

# features driving separation of samples in PLS-DA 
taxonomy <- read_qza('taxonomy.qza')$data %>%
  dplyr::select(-Confidence)

plsda_lung$loadings$X %>%
  data.frame() %>%
  dplyr::select(comp1, comp2) %>%
  rownames_to_column('OTU') %>%
  left_join(taxonomy, by = c('OTU' = 'Feature.ID')) %>%
  arrange(-abs(comp1)) %>%
  write_tsv('plsda_lung_loadings.tsv')

plsda_faecal$loadings$X %>%
  data.frame() %>%
  dplyr::select(comp1, comp2) %>%
  rownames_to_column('OTU') %>%
  left_join(taxonomy, by = c('OTU' = 'Feature.ID')) %>%
  arrange(-abs(comp1)) %>%
  write_tsv('plsda_faecal_loadings.tsv')

plsda_ear$loadings$X %>%
  data.frame() %>%
  dplyr::select(comp1, comp2) %>%
  rownames_to_column('OTU') %>%
  left_join(taxonomy, by = c('OTU' = 'Feature.ID')) %>%
  arrange(-abs(comp1)) %>%
  write_tsv('plsda_ear_loadings.tsv')

df <- read_tsv('plsda_lung_loadings.tsv')
p <- df %>%
   pivot_longer(cols = starts_with('comp'), names_to = 'Component', values_to = 'Loading') %>%
   group_by(Component) %>%
   slice_max(order_by = abs(Loading), n = 10) %>%
   ggplot(aes(x = Taxon, y = Loading)) +
   ggforce::facet_col(~Component, scales = 'free', space = 'free') +
   geom_bar(colour = 'black', fill = 'grey', stat = 'identity') +
   theme_bw() +
   coord_flip()
pdf('plsda_lung_loadings.pdf', width = 10, height = 5)
p
dev.off()

df <- read_tsv('plsda_faecal_loadings.tsv')
p <- df %>%
  pivot_longer(cols = starts_with('comp'), names_to = 'Component', values_to = 'Loading') %>%
  group_by(Component) %>%
  slice_max(order_by = abs(Loading), n = 10) %>%
  ggplot(aes(x = Taxon, y = Loading)) +
  ggforce::facet_col(~Component, scales = 'free', space = 'free') +
  geom_bar(colour = 'black', fill = 'grey', stat = 'identity') +
  theme_bw() +
  coord_flip()
pdf('plsda_faecal_loadings.pdf', width = 10, height = 4)
p
dev.off()

df <- read_tsv('plsda_ear_loadings.tsv')
p <- df %>%
  pivot_longer(cols = starts_with('comp'), names_to = 'Component', values_to = 'Loading') %>%
  group_by(Component) %>%
  slice_max(order_by = abs(Loading), n = 10) %>%
  ggplot(aes(x = Taxon, y = Loading)) +
  ggforce::facet_col(~Component, scales = 'free', space = 'free') +
  geom_bar(colour = 'black', fill = 'grey', stat = 'identity') +
  theme_bw() +
  coord_flip()
pdf('plsda_ear_loadings.pdf', width = 10, height = 5)
p
dev.off()
   
# extracting/plotting RA of OTUs with strongest X loadings in each body site
phylo <- readRDS('phylo_raw.RDS')

otutable_ra <- otu_table(phylo) %>%
  data.frame() %>%
  rownames_to_column('OTU') %>%
  pivot_longer(cols = -OTU, names_to = 'Sample', values_to = 'Count') %>%
  group_by(Sample) %>%
  mutate(RA = 100*(Count/sum(Count))) %>%
  ungroup() %>%
  left_join(metadata, by = 'Sample') %>% 
  arrange(Sample, OTU)

loadings_ear <- read_tsv('plsda_ear_loadings.tsv') %>%
  pivot_longer(cols = starts_with('comp'), names_to = 'Component', values_to = 'Loading') %>%
  filter(Component == 'comp1') %>%
  slice_max(order_by = abs(Loading), n = 10)
loadings_faecal <- read_tsv('plsda_faecal_loadings.tsv') %>%
  pivot_longer(cols = starts_with('comp'), names_to = 'Component', values_to = 'Loading') %>%
  filter(Component == 'comp1') %>%
  slice_max(order_by = abs(Loading), n = 10)
loadings_lung <- read_tsv('plsda_lung_loadings.tsv') %>%
  pivot_longer(cols = starts_with('comp'), names_to = 'Component', values_to = 'Loading') %>%
  filter(Component == 'comp1') %>%
  slice_max(order_by = abs(Loading), n = 10)

ra_ear <- otutable_ra %>% 
  filter(Source == 'Ear' & OTU %in% loadings_ear$OTU)
ra_faecal <- otutable_ra %>% 
  filter(Source == 'Faecal' & OTU %in% loadings_faecal$OTU)
ra_lung <- otutable_ra %>% 
  filter(Source == 'Lung' & OTU %in% loadings_lung$OTU)

taxonomy <- read_qza('taxonomy.qza')$data %>%
  dplyr::select(-Confidence)

list(Ear = ra_ear %>%
       pivot_wider(id_cols = OTU, names_from = Sample, values_from = RA, values_fill = 0) %>%
       left_join(taxonomy, by = c('OTU' = 'Feature.ID')) %>%
                   relocate(Taxon, .after = OTU),
     Faecal = ra_faecal %>%
       pivot_wider(id_cols = OTU, names_from = Sample, values_from = RA, values_fill = 0) %>%
       left_join(taxonomy, by = c('OTU' = 'Feature.ID')) %>%
       relocate(Taxon, .after = OTU),
     Lung = ra_lung %>%
       pivot_wider(id_cols = OTU, names_from = Sample, values_from = RA, values_fill = 0) %>%
       left_join(taxonomy, by = c('OTU' = 'Feature.ID')) %>%
       relocate(Taxon, .after = OTU)) %>%
  writexl::write_xlsx('ra_toploadings.xlsx')

# trying Taxon Set Enrichment Analysis (applying the GSEA approach to OTUs and taxonomic labels)
# looking for taxa (genera, families, etc.) that are enriched in either end of the position along the PLSDA X axis
library(clusterProfiler)
loadings_ear <- read_tsv('plsda_ear_loadings.tsv') 
loadings_faecal <- read_tsv('plsda_faecal_loadings.tsv') 
loadings_lung <- read_tsv('plsda_lung_loadings.tsv') 

taxonomy <- read_qza('taxonomy.qza')$data %>%
  dplyr::select(-Confidence)

genus2otu <- taxonomy %>% dplyr::select(2,1)%>% mutate(Taxon = str_remove(Taxon, '; s__.*'))
family2otu <- taxonomy %>% dplyr::select(2,1)%>% mutate(Taxon = str_remove(Taxon, '; g__.*'))
order2otu <- taxonomy %>% dplyr::select(2,1)%>% mutate(Taxon = str_remove(Taxon, '; f__.*'))
class2otu <- taxonomy %>% dplyr::select(2,1) %>% mutate(Taxon = str_remove(Taxon, '; o__.*'))
phylum2otu <-  taxonomy %>% dplyr::select(2,1) %>% mutate(Taxon = str_remove(Taxon, '; c__.*'))

otu_background <- loadings_ear$OTU 
loadings_ear <- loadings_ear %>%
  arrange(-comp1)
otu_list <- loadings_ear$comp1
names(otu_list) <- loadings_ear$OTU

tsea_genus <- GSEA(geneList = otu_list,
                   TERM2GENE = genus2otu,
                   TERM2NAME = NA,
                   minGSSize = 5,
                   maxGSSize = length(otu_background) / 2,
                   pvalueCutoff = 1,
                   pAdjustMethod = "fdr") %>%
  data.frame(check.names = F) 

tsea_family <- GSEA(geneList = otu_list,
                   TERM2GENE = family2otu,
                   TERM2NAME = NA,
                   minGSSize = 5,
                   maxGSSize = length(otu_background) / 2,
                   pvalueCutoff = 1,
                   pAdjustMethod = "fdr") %>%
  data.frame(check.names = F) 

tsea_order <- GSEA(geneList = otu_list,
                   TERM2GENE = order2otu,
                   TERM2NAME = NA,
                   minGSSize = 5,
                   maxGSSize = length(otu_background) / 2,
                   pvalueCutoff = 1,
                   pAdjustMethod = "fdr") %>%
  data.frame(check.names = F) 

tsea_class <- GSEA(geneList = otu_list,
                   TERM2GENE = class2otu,
                   TERM2NAME = NA,
                   minGSSize = 5,
                   maxGSSize = length(otu_background) / 2,
                   pvalueCutoff = 1,
                   pAdjustMethod = "fdr") %>%
  data.frame(check.names = F) 

tsea_phylum <- GSEA(geneList = otu_list,
                   TERM2GENE = phylum2otu,
                   TERM2NAME = NA,
                   minGSSize = 5,
                   maxGSSize = length(otu_background) / 2,
                   pvalueCutoff = 1,
                   pAdjustMethod = "fdr") %>%
  data.frame(check.names = F)

list(Genus = tsea_genus,
     Family = tsea_family,
     Order = tsea_order,
     Class = tsea_class,
     Phylum = tsea_phylum) %>%
  writexl::write_xlsx('tsea_ear.xlsx')

otu_background <- loadings_faecal$OTU 
loadings_faecal <- loadings_faecal %>%
  arrange(-comp1)
otu_list <- loadings_faecal$comp1
names(otu_list) <- loadings_faecal$OTU

tsea_genus <- GSEA(geneList = otu_list,
                   TERM2GENE = genus2otu,
                   TERM2NAME = NA,
                   minGSSize = 5,
                   maxGSSize = length(otu_background) / 2,
                   pvalueCutoff = 1,
                   pAdjustMethod = "fdr") %>%
  data.frame(check.names = F) 

tsea_family <- GSEA(geneList = otu_list,
                    TERM2GENE = family2otu,
                    TERM2NAME = NA,
                    minGSSize = 5,
                    maxGSSize = length(otu_background) / 2,
                    pvalueCutoff = 1,
                    pAdjustMethod = "fdr") %>%
  data.frame(check.names = F) 

tsea_order <- GSEA(geneList = otu_list,
                   TERM2GENE = order2otu,
                   TERM2NAME = NA,
                   minGSSize = 5,
                   maxGSSize = length(otu_background) / 2,
                   pvalueCutoff = 1,
                   pAdjustMethod = "fdr") %>%
  data.frame(check.names = F)

tsea_class <- GSEA(geneList = otu_list,
                   TERM2GENE = class2otu,
                   TERM2NAME = NA,
                   minGSSize = 5,
                   maxGSSize = length(otu_background) / 2,
                   pvalueCutoff = 1,
                   pAdjustMethod = "fdr") %>%
  data.frame(check.names = F) 

tsea_phylum <- GSEA(geneList = otu_list,
                    TERM2GENE = phylum2otu,
                    TERM2NAME = NA,
                    minGSSize = 5,
                    maxGSSize = length(otu_background) / 2,
                    pvalueCutoff = 1,
                    pAdjustMethod = "fdr") %>%
  data.frame(check.names = F)

list(Genus = tsea_genus,
     Family = tsea_family,
     Phylum = tsea_phylum) %>%
  writexl::write_xlsx('tsea_faecal.xlsx')

otu_background <- loadings_lung$OTU 
loadings_lung <- loadings_lung %>%
  arrange(-comp1)
otu_list <- loadings_lung$comp1
names(otu_list) <- loadings_lung$OTU

tsea_genus <- GSEA(geneList = otu_list,
                   TERM2GENE = genus2otu,
                   TERM2NAME = NA,
                   minGSSize = 5,
                   maxGSSize = length(otu_background) / 2,
                   pvalueCutoff = 1,
                   pAdjustMethod = "fdr") %>%
  data.frame(check.names = F) 

tsea_family <- GSEA(geneList = otu_list,
                    TERM2GENE = family2otu,
                    TERM2NAME = NA,
                    minGSSize = 5,
                    maxGSSize = length(otu_background) / 2,
                    pvalueCutoff = 1,
                    pAdjustMethod = "fdr") %>%
  data.frame(check.names = F) 

tsea_order <- GSEA(geneList = otu_list,
                   TERM2GENE = order2otu,
                   TERM2NAME = NA,
                   minGSSize = 5,
                   maxGSSize = length(otu_background) / 2,
                   pvalueCutoff = 1,
                   pAdjustMethod = "fdr") %>%
  data.frame(check.names = F) 

tsea_class <- GSEA(geneList = otu_list,
                   TERM2GENE = class2otu,
                   TERM2NAME = NA,
                   minGSSize = 5,
                   maxGSSize = length(otu_background) / 2,
                   pvalueCutoff = 1,
                   pAdjustMethod = "fdr") %>%
  data.frame(check.names = F)

tsea_phylum <- GSEA(geneList = otu_list,
                    TERM2GENE = phylum2otu,
                    TERM2NAME = NA,
                    minGSSize = 5,
                    maxGSSize = length(otu_background) / 2,
                    pvalueCutoff = 1,
                    pAdjustMethod = "fdr") %>%
  data.frame(check.names = F)

list(Genus = tsea_genus,
     Family = tsea_family,
     Order = tsea_order,
     Class = tsea_class,
     Phylum = tsea_phylum) %>%
  writexl::write_xlsx('tsea_lung.xlsx')

# plotting loading distribution
library(tidyverse)
library(ggridges)

tsea <- readxl::read_xlsx('tsea_ear.xlsx', sheet = 'Phylum')
p <- loadings_ear %>%
  mutate(Phylum = str_remove(Taxon, '; c__.*')) %>%
  filter(Phylum %in% tsea$Description) %>%
  left_join(tsea, by = c('Phylum' = 'Description')) %>%
  mutate(Sig = if_else(qvalue <= 0.05 & NES > 0, 'Up', 
                       if_else(qvalue <= 0.05 & NES < 0, 'Down', 'No'))) %>%
  ggplot(aes(x = comp1, y = Phylum, group = Phylum)) +
  geom_density_ridges(aes(fill = Sig)) +
  theme_bw() +
  scale_fill_manual(values = c('No' = 'grey',
                               'Up' = 'forestgreen',
                               'Down' = 'firebrick')) +
  labs(x = 'Loading', y = 'Phylum', title = 'Loading Distribution: Component 1')
pdf('tsea_ear_loadings_phylum.pdf', width = 10, height = 4)
p
dev.off()

tsea <- readxl::read_xlsx('tsea_ear.xlsx', sheet = 'Class')
p <- loadings_ear %>%
  mutate(Class = str_remove(Taxon, '; o__.*')) %>%
  filter(Class %in% tsea$Description) %>%
  left_join(tsea, by = c('Class' = 'Description')) %>%
  mutate(Sig = if_else(qvalue <= 0.05 & NES > 0, 'Up', 
                       if_else(qvalue <= 0.05 & NES < 0, 'Down', 'No'))) %>%
  ggplot(aes(x = comp1, y = Class, group = Class)) +
  geom_density_ridges(aes(fill = Sig)) +
  theme_bw() +
  scale_fill_manual(values = c('No' = 'grey',
                               'Up' = 'forestgreen',
                               'Down' = 'firebrick')) +
  labs(x = 'Loading', y = 'Class', title = 'Loading Distribution: Component 1')
pdf('tsea_ear_loadings_class.pdf', width = 10, height = 4)
p
dev.off()

tsea <- readxl::read_xlsx('tsea_ear.xlsx', sheet = 'Order')
p <- loadings_ear %>%
  mutate(Order = str_remove(Taxon, '; f__.*')) %>%
  filter(Order %in% tsea$Description) %>%
  left_join(tsea, by = c('Order' = 'Description')) %>%
  mutate(Sig = if_else(qvalue <= 0.05 & NES > 0, 'Up', 
                       if_else(qvalue <= 0.05 & NES < 0, 'Down', 'No'))) %>%
  ggplot(aes(x = comp1, y = Order, group = Order)) +
  geom_density_ridges(aes(fill = Sig)) +
  theme_bw() +
  scale_fill_manual(values = c('No' = 'grey',
                               'Up' = 'forestgreen',
                               'Down' = 'firebrick')) +
  labs(x = 'Loading', y = 'Order', title = 'Loading Distribution: Component 1')
pdf('tsea_ear_loadings_order.pdf', width = 10, height = 4)
p
dev.off()

tsea <- readxl::read_xlsx('tsea_ear.xlsx', sheet = 'Family')
p <- loadings_ear %>%
  mutate(Family = str_remove(Taxon, '; g__.*')) %>%
  filter(Family %in% tsea$Description) %>%
  left_join(tsea, by = c('Family' = 'Description')) %>%
  mutate(Sig = if_else(qvalue <= 0.05 & NES > 0, 'Up', 
                       if_else(qvalue <= 0.05 & NES < 0, 'Down', 'No'))) %>%
  ggplot(aes(x = comp1, y = Family, group = Family)) +
  geom_density_ridges(aes(fill = Sig)) +
  theme_bw() +
  scale_fill_manual(values = c('No' = 'grey',
                               'Up' = 'forestgreen',
                               'Down' = 'firebrick')) +
  labs(x = 'Loading', y = 'Family', title = 'Loading Distribution: Component 1')
pdf('tsea_ear_loadings_family.pdf', width = 10, height = 4)
p
dev.off()

tsea <- readxl::read_xlsx('tsea_ear.xlsx', sheet = 'Genus')
p <- loadings_ear %>%
  mutate(Genus = str_remove(Taxon, '; s__.*')) %>%
  filter(Genus %in% tsea$Description) %>%
  left_join(tsea, by = c('Genus' = 'Description')) %>%
  mutate(Sig = if_else(qvalue <= 0.05 & NES > 0, 'Up', 
                       if_else(qvalue <= 0.05 & NES < 0, 'Down', 'No'))) %>%
  ggplot(aes(x = comp1, y = Genus, group = Genus)) +
  geom_density_ridges(aes(fill = Sig)) +
  theme_bw() +
  scale_fill_manual(values = c('No' = 'grey',
                               'Up' = 'forestgreen',
                               'Down' = 'firebrick')) +
  labs(x = 'Loading', y = 'Genus', title = 'Loading Distribution: Component 1')
pdf('tsea_ear_loadings_genus.pdf', width = 10, height = 4)
p
dev.off()

tsea <- readxl::read_xlsx('tsea_faecal.xlsx', sheet = 'Phylum')
p <- loadings_faecal %>%
  mutate(Phylum = str_remove(Taxon, '; c__.*')) %>%
  filter(Phylum %in% tsea$Description) %>%
  left_join(tsea, by = c('Phylum' = 'Description')) %>%
  mutate(Sig = if_else(qvalue <= 0.05 & NES > 0, 'Up', 
                       if_else(qvalue <= 0.05 & NES < 0, 'Down', 'No'))) %>%
  ggplot(aes(x = comp1, y = Phylum, group = Phylum)) +
  geom_density_ridges(aes(fill = Sig)) +
  theme_bw() +
  scale_fill_manual(values = c('No' = 'grey',
                               'Up' = 'forestgreen',
                               'Down' = 'firebrick')) +
  labs(x = 'Loading', y = 'Phylum', title = 'Loading Distribution: Component 1')
pdf('tsea_faecal_loadings_phylum.pdf', width = 10, height = 4)
p
dev.off()

tsea <- readxl::read_xlsx('tsea_faecal.xlsx', sheet = 'Class')
p <- loadings_faecal %>%
  mutate(Class = str_remove(Taxon, '; o__.*')) %>%
  filter(Class %in% tsea$Description) %>%
  left_join(tsea, by = c('Class' = 'Description')) %>%
  mutate(Sig = if_else(qvalue <= 0.05 & NES > 0, 'Up', 
                       if_else(qvalue <= 0.05 & NES < 0, 'Down', 'No'))) %>%
  ggplot(aes(x = comp1, y = Class, group = Class)) +
  geom_density_ridges(aes(fill = Sig)) +
  theme_bw() +
  scale_fill_manual(values = c('No' = 'grey',
                               'Up' = 'forestgreen',
                               'Down' = 'firebrick')) +
  labs(x = 'Loading', y = 'Class', title = 'Loading Distribution: Component 1')
pdf('tsea_faecal_loadings_class.pdf', width = 10, height = 4)
p
dev.off()

tsea <- readxl::read_xlsx('tsea_faecal.xlsx', sheet = 'Order')
p <- loadings_faecal %>%
  mutate(Order = str_remove(Taxon, '; f__.*')) %>%
  filter(Order %in% tsea$Description) %>%
  left_join(tsea, by = c('Order' = 'Description')) %>%
  mutate(Sig = if_else(qvalue <= 0.05 & NES > 0, 'Up', 
                       if_else(qvalue <= 0.05 & NES < 0, 'Down', 'No'))) %>%
  ggplot(aes(x = comp1, y = Order, group = Order)) +
  geom_density_ridges(aes(fill = Sig)) +
  theme_bw() +
  scale_fill_manual(values = c('No' = 'grey',
                               'Up' = 'forestgreen',
                               'Down' = 'firebrick')) +
  labs(x = 'Loading', y = 'Order', title = 'Loading Distribution: Component 1')
pdf('tsea_faecal_loadings_order.pdf', width = 10, height = 4)
p
dev.off()

tsea <- readxl::read_xlsx('tsea_faecal.xlsx', sheet = 'Family')
p <- loadings_faecal %>%
  mutate(Family = str_remove(Taxon, '; g__.*')) %>%
  filter(Family %in% tsea$Description) %>%
  left_join(tsea, by = c('Family' = 'Description')) %>%
  mutate(Sig = if_else(qvalue <= 0.05 & NES > 0, 'Up', 
                       if_else(qvalue <= 0.05 & NES < 0, 'Down', 'No'))) %>%
  ggplot(aes(x = comp1, y = Family, group = Family)) +
  geom_density_ridges(aes(fill = Sig)) +
  theme_bw() +
  scale_fill_manual(values = c('No' = 'grey',
                               'Up' = 'forestgreen',
                               'Down' = 'firebrick')) +
  labs(x = 'Loading', y = 'Family', title = 'Loading Distribution: Component 1')
pdf('tsea_faecal_loadings_family.pdf', width = 10, height = 4)
p
dev.off()

tsea <- readxl::read_xlsx('tsea_faecal.xlsx', sheet = 'Genus')
p <- loadings_faecal %>%
  mutate(Genus = str_remove(Taxon, '; s__.*')) %>%
  filter(Genus %in% tsea$Description) %>%
  left_join(tsea, by = c('Genus' = 'Description')) %>%
  mutate(Sig = if_else(qvalue <= 0.05 & NES > 0, 'Up', 
                       if_else(qvalue <= 0.05 & NES < 0, 'Down', 'No'))) %>%
  ggplot(aes(x = comp1, y = Genus, group = Genus)) +
  geom_density_ridges(aes(fill = Sig)) +
  theme_bw() +
  scale_fill_manual(values = c('No' = 'grey',
                               'Up' = 'forestgreen',
                               'Down' = 'firebrick')) +
  labs(x = 'Loading', y = 'Genus', title = 'Loading Distribution: Component 1')
pdf('tsea_faecal_loadings_genus.pdf', width = 10, height = 4)
p
dev.off()

tsea <- readxl::read_xlsx('tsea_lung.xlsx', sheet = 'Phylum')
p <- loadings_lung %>%
  mutate(Phylum = str_remove(Taxon, '; c__.*')) %>%
  filter(Phylum %in% tsea$Description) %>%
  left_join(tsea, by = c('Phylum' = 'Description')) %>%
  mutate(Sig = if_else(qvalue <= 0.05 & NES > 0, 'Up', 
                       if_else(qvalue <= 0.05 & NES < 0, 'Down', 'No'))) %>%
  ggplot(aes(x = comp1, y = Phylum, group = Phylum)) +
  geom_density_ridges(aes(fill = Sig)) +
  theme_bw() +
  scale_fill_manual(values = c('No' = 'grey',
                               'Up' = 'forestgreen',
                               'Down' = 'firebrick')) +
  labs(x = 'Loading', y = 'Phylum', title = 'Loading Distribution: Component 1')
pdf('tsea_lung_loadings_phylum.pdf', width = 10, height = 4)
p
dev.off()

tsea <- readxl::read_xlsx('tsea_lung.xlsx', sheet = 'Class')
p <- loadings_lung %>%
  mutate(Class = str_remove(Taxon, '; o__.*')) %>%
  filter(Class %in% tsea$Description) %>%
  left_join(tsea, by = c('Class' = 'Description')) %>%
  mutate(Sig = if_else(qvalue <= 0.05 & NES > 0, 'Up', 
                       if_else(qvalue <= 0.05 & NES < 0, 'Down', 'No'))) %>%
  ggplot(aes(x = comp1, y = Class, group = Class)) +
  geom_density_ridges(aes(fill = Sig)) +
  theme_bw() +
  scale_fill_manual(values = c('No' = 'grey',
                               'Up' = 'forestgreen',
                               'Down' = 'firebrick')) +
  labs(x = 'Loading', y = 'Class', title = 'Loading Distribution: Component 1')
pdf('tsea_lung_loadings_class.pdf', width = 10, height = 4)
p
dev.off()

tsea <- readxl::read_xlsx('tsea_lung.xlsx', sheet = 'Order')
p <- loadings_lung %>%
  mutate(Order = str_remove(Taxon, '; f__.*')) %>%
  filter(Order %in% tsea$Description) %>%
  left_join(tsea, by = c('Order' = 'Description')) %>%
  mutate(Sig = if_else(qvalue <= 0.05 & NES > 0, 'Up', 
                       if_else(qvalue <= 0.05 & NES < 0, 'Down', 'No'))) %>%
  ggplot(aes(x = comp1, y = Order, group = Order)) +
  geom_density_ridges(aes(fill = Sig)) +
  theme_bw() +
  scale_fill_manual(values = c('No' = 'grey',
                               'Up' = 'forestgreen',
                               'Down' = 'firebrick')) +
  labs(x = 'Loading', y = 'Order', title = 'Loading Distribution: Component 1')
pdf('tsea_lung_loadings_order.pdf', width = 10, height = 4)
p
dev.off()

tsea <- readxl::read_xlsx('tsea_lung.xlsx', sheet = 'Family')
p <- loadings_lung %>%
  mutate(Family = str_remove(Taxon, '; g__.*')) %>%
  filter(Family %in% tsea$Description) %>%
  left_join(tsea, by = c('Family' = 'Description')) %>%
  mutate(Sig = if_else(qvalue <= 0.05 & NES > 0, 'Up', 
                       if_else(qvalue <= 0.05 & NES < 0, 'Down', 'No'))) %>%
  ggplot(aes(x = comp1, y = Family, group = Family)) +
  geom_density_ridges(aes(fill = Sig)) +
  theme_bw() +
  scale_fill_manual(values = c('No' = 'grey',
                               'Up' = 'forestgreen',
                               'Down' = 'firebrick')) +
  labs(x = 'Loading', y = 'Family', title = 'Loading Distribution: Component 1')
pdf('tsea_lung_loadings_family.pdf', width = 10, height = 4)
p
dev.off()

tsea <- readxl::read_xlsx('tsea_lung.xlsx', sheet = 'Genus')
p <- loadings_lung %>%
  mutate(Genus = str_remove(Taxon, '; s__.*')) %>%
  filter(Genus %in% tsea$Description) %>%
  left_join(tsea, by = c('Genus' = 'Description')) %>%
  mutate(Sig = if_else(qvalue <= 0.05 & NES > 0, 'Up', 
                       if_else(qvalue <= 0.05 & NES < 0, 'Down', 'No'))) %>%
  ggplot(aes(x = comp1, y = Genus, group = Genus)) +
  geom_density_ridges(aes(fill = Sig)) +
  theme_bw() +
  scale_fill_manual(values = c('No' = 'grey',
                               'Up' = 'forestgreen',
                               'Down' = 'firebrick')) +
  labs(x = 'Loading', y = 'Genus', title = 'Loading Distribution: Component 1')
pdf('tsea_lung_loadings_genus.pdf', width = 10, height = 4)
p
dev.off()


# outputting genus table
phylo <- readRDS('phylo_filt.RDS')

taxonomy <- read_qza('taxonomy.qza')$data %>%
  dplyr::select(-Confidence) %>%
  mutate(Genus = str_remove(Taxon, '; s__.*')) %>%
  rename(OTU = 1) %>%
  select(OTU, Genus)

genus <- phylo %>%
  otu_table() %>%
  data.frame(check.names = F) %>%
  rownames_to_column('OTU') %>%
  left_join(taxonomy) %>%
  pivot_longer(cols = starts_with('DMG'), names_to = 'Sample', values_to = 'Count') %>%
  group_by(Sample, Genus) %>%
  summarise(Count = sum(Count)) %>%
  ungroup() %>%
  group_by(Sample) %>%
  mutate(RA = 100*(Count/sum(Count))) %>%
  ungroup()

genus %>%
  pivot_wider(id_cols = Genus, names_from = Sample, values_from = Count, values_fill = 0) %>%
  write_tsv('genus_count.tsv')
genus %>%
  pivot_wider(id_cols = Genus, names_from = Sample, values_from = RA, values_fill = 0) %>%
  write_tsv('genus_ra.tsv')
  
  
