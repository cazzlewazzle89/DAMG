library('tidyverse')
library('rstatix')
library('vegan')

# importing metadata
metadata <- readxl::read_xlsx('metadata.xlsx') 

# plotting read count
reads <- readxl::read_xlsx('readcount.xlsx') %>%
  left_join(metadata)

reads %>%
  writexl::write_xlsx('readcounts.xlsx')

# importing Emu profiles
df <- read_tsv('emu-combined-species.tsv')

df <- df %>%
  select(species, starts_with('DMG')) %>%
  rename(Species = species)
df[is.na(df)] <- 0

df_long <- df %>%
  pivot_longer(cols = -Species, names_to = 'Sample', values_to = 'RA')

meanRA <- df_long %>%
  group_by(Species) %>%
  summarise(MeanRA = mean(RA)) %>%
  arrange(-MeanRA)

speciestokeep <- meanRA %>%
  filter(MeanRA > 0.001)

df_wide <- df_long %>%
  filter(Species %in% speciestokeep$Species) %>%
  group_by(Species) %>%
  mutate(MeanRA = mean(RA)) %>%
  arrange(-MeanRA) %>%
  pivot_wider(id_cols = Sample, names_from = Species, values_from = RA, values_fill = 0)

df_wide %>% 
  writexl::write_xlsx('speciestable.xlsx')

df_wide <- df_wide %>%
  column_to_rownames('Sample')

# calculating diversity metrics
shannon <- diversity(df_wide, index = 'shannon')
bray <- vegdist(df_wide, method = 'bray')
pcoa_bray <- cmdscale(bray, k = 2, eig = T)
points_bray <- pcoa_bray$points %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(metadata)
eigs <- 100*(pcoa_bray$eig/sum(pcoa_bray$eig))

# defining colour palettes
colours <- hues::iwanthue(4)
treatmentColours <- colours[1:2]
names(treatmentColours) <- unique(metadata$treatment)
sexColours <- colours[3:4]
names(sexColours) <- unique(metadata$sex)

# all data together
p_shannon <- shannon %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  rename('Shannon' = 2) %>%
  left_join(metadata) %>%
  ggplot(aes(x = treatment, y = Shannon)) +
  geom_boxplot(aes(colour = treatment), fill = NA, outlier.shape = NA, show.legend = F) +
  geom_point(aes(colour = treatment, shape = sex), size = 3, 
             position = position_jitterdodge(jitter.width = 0.5)) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(shape = 23))) +
  scale_shape_manual(values = list(21, 22)) +
  scale_color_manual(values = treatmentColours)

p_shannon %>% saveRDS('shannon_all.RDS')
pdf('shannon_all.pdf', width = 5, height = 4)
p_shannon
dev.off()

shannon %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  rename('Shannon' = 2) %>%
  left_join(metadata) %>%
  t_test(Shannon ~ treatment)

shannon %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  rename('Shannon' = 2) %>%
  left_join(metadata) %>%
  anova_test(Shannon ~ treatment + sex)

p_bray <- points_bray %>%
  ggplot(aes(x = X1, y = X2)) +
  geom_point(aes(fill = treatment, shape = sex), colour = 'black', size = 3) +
  theme_bw() +
  scale_shape_manual(values = list(21, 22)) +
  scale_fill_manual(values = treatmentColours) +
  guides(fill = guide_legend(override.aes = list(shape = 23))) +
  labs(x = paste0('PC1 [', round(eigs[1], 2), '% Var. Explained]'),
       y = paste0('PC2 [', round(eigs[2], 2), '% Var. Explained]'))

p_bray %>% saveRDS('bray_all.RDS')
pdf('bray_all.pdf', width = 6, height = 4)
p_bray
dev.off()

perm_meta <- metadata %>% filter(Sample %in% colnames(df))
perm_bray <- usedist::dist_subset(bray, perm_meta$Sample)
adonis2(perm_bray ~ perm_meta$treatment)
adonis2(perm_bray ~ perm_meta$treatment * perm_meta$sex)

# male mice only
p_shannon <- shannon %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  rename('Shannon' = 2) %>%
  left_join(metadata) %>%
  filter(sex == 'm') %>%
  ggplot(aes(x = treatment, y = Shannon)) +
  geom_boxplot(aes(colour = treatment), fill = NA, outlier.shape = NA, show.legend = F) +
  geom_point(aes(colour = treatment, shape = sex), size = 3, 
             position = position_jitterdodge(jitter.width = 0.5)) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(shape = 23))) +
  scale_shape_manual(values = list(21, 22)) +
  scale_color_manual(values = treatmentColours)

p_shannon %>% saveRDS('shannon_male.RDS')
pdf('shannon_male.pdf', width = 5, height = 4)
p_shannon
dev.off()

shannon %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  rename('Shannon' = 2) %>%
  left_join(metadata) %>%
  filter(sex == 'm') %>%
  t_test(Shannon ~ treatment)

p_bray <- points_bray %>%
  filter(sex == 'm') %>%
  ggplot(aes(x = X1, y = X2)) +
  geom_point(aes(fill = treatment), shape = 21, colour = 'black', size = 3) +
  theme_bw() +
  scale_fill_manual(values = treatmentColours) +
  labs(x = paste0('PC1 [', round(eigs[1], 2), '% Var. Explained]'),
       y = paste0('PC2 [', round(eigs[2], 2), '% Var. Explained]'))

p_bray %>% saveRDS('bray_male.RDS')
pdf('bray_male.pdf', width = 6, height = 4)
p_bray
dev.off()

perm_meta <- metadata %>% filter(Sample %in% colnames(df) & sex == 'm')
perm_bray <- usedist::dist_subset(bray, perm_meta$Sample)
adonis2(perm_bray ~ perm_meta$treatment)

# female mice only
p_shannon <- shannon %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  rename('Shannon' = 2) %>%
  left_join(metadata) %>%
  filter(sex == 'f') %>%
  ggplot(aes(x = treatment, y = Shannon)) +
  geom_boxplot(aes(colour = treatment), fill = NA, outlier.shape = NA, show.legend = F) +
  geom_point(aes(colour = treatment, shape = sex), size = 3, 
             position = position_jitterdodge(jitter.width = 0.5)) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(shape = 23))) +
  scale_shape_manual(values = list(21, 22)) +
  scale_color_manual(values = treatmentColours)

p_shannon %>% saveRDS('shannon_female.RDS')
pdf('shannon_female.pdf', width = 5, height = 4)
p_shannon
dev.off()

shannon %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  rename('Shannon' = 2) %>%
  left_join(metadata) %>%
  filter(sex == 'f') %>%
  t_test(Shannon ~ treatment)

p_bray <- points_bray %>%
  filter(sex == 'f') %>%
  ggplot(aes(x = X1, y = X2)) +
  geom_point(aes(fill = treatment), shape = 21, colour = 'black', size = 3) +
  theme_bw() +
  scale_fill_manual(values = treatmentColours) +
  labs(x = paste0('PC1 [', round(eigs[1], 2), '% Var. Explained]'),
       y = paste0('PC2 [', round(eigs[2], 2), '% Var. Explained]'))

p_bray %>% saveRDS('bray_female.RDS')
pdf('bray_female.pdf', width = 6, height = 4)
p_bray
dev.off()

perm_meta <- metadata %>% filter(Sample %in% colnames(df) & sex == 'f')
perm_bray <- usedist::dist_subset(bray, perm_meta$Sample)
adonis2(perm_bray ~ perm_meta$treatment)

# comparing by sex within treatment group
# Saline mice
p_shannon <- shannon %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  rename('Shannon' = 2) %>%
  left_join(metadata) %>%
  filter(treatment == 'Saline') %>%
  ggplot(aes(x = sex, y = Shannon)) +
  geom_boxplot(aes(colour = sex), fill = NA, outlier.shape = NA, show.legend = F) +
  geom_point(aes(colour = sex, shape = sex), size = 3, 
             position = position_jitterdodge(jitter.width = 0.5)) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(shape = 23))) +
  scale_shape_manual(values = list(21, 22)) +
  scale_color_manual(values = sexColours)

p_shannon %>% saveRDS('shannon_saline.RDS')
pdf('shannon_saline.pdf', width = 5, height = 4)
p_shannon
dev.off()

shannon %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  rename('Shannon' = 2) %>%
  left_join(metadata) %>%
  filter(treatment == 'Saline') %>%
  t_test(Shannon ~ sex)

p_bray <- points_bray %>%
  filter(treatment == 'Saline') %>%
  ggplot(aes(x = X1, y = X2)) +
  geom_point(aes(fill = sex), shape = 21, colour = 'black', size = 3) +
  theme_bw() +
  scale_fill_manual(values = sexColours) +
  labs(x = paste0('PC1 [', round(eigs[1], 2), '% Var. Explained]'),
       y = paste0('PC2 [', round(eigs[2], 2), '% Var. Explained]'))

p_bray %>% saveRDS('bray_saline.RDS')
pdf('bray_saline.pdf', width = 6, height = 4)
p_bray
dev.off()

perm_meta <- metadata %>% filter(Sample %in% colnames(df) & treatment == 'Saline')
perm_bray <- usedist::dist_subset(bray, perm_meta$Sample)
adonis2(perm_bray ~ perm_meta$sex)

# Poly I:C mice
p_shannon <- shannon %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  rename('Shannon' = 2) %>%
  left_join(metadata) %>%
  filter(treatment == 'Poly I:C') %>%
  ggplot(aes(x = sex, y = Shannon)) +
  geom_boxplot(aes(colour = sex), fill = NA, outlier.shape = NA, show.legend = F) +
  geom_point(aes(colour = sex, shape = sex), size = 3, 
             position = position_jitterdodge(jitter.width = 0.5)) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(shape = 23))) +
  scale_shape_manual(values = list(21, 22)) +
  scale_color_manual(values = sexColours)

p_shannon %>% saveRDS('shannon_polyIC.RDS')
pdf('shannon_polyIC.pdf', width = 5, height = 4)
p_shannon
dev.off()

shannon %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  rename('Shannon' = 2) %>%
  left_join(metadata) %>%
  filter(treatment == 'Poly I:C') %>%
  t_test(Shannon ~ sex)

p_bray <- points_bray %>%
  filter(treatment == 'Poly I:C') %>%
  ggplot(aes(x = X1, y = X2)) +
  geom_point(aes(fill = sex), shape = 21, colour = 'black', size = 3) +
  theme_bw() +
  scale_fill_manual(values = sexColours) +
  labs(x = paste0('PC1 [', round(eigs[1], 2), '% Var. Explained]'),
       y = paste0('PC2 [', round(eigs[2], 2), '% Var. Explained]'))

p_bray %>% saveRDS('bray_polyIC.RDS')
pdf('bray_polyIC.pdf', width = 6, height = 4)
p_bray
dev.off()

perm_meta <- metadata %>% filter(Sample %in% colnames(df) & treatment == 'Poly I:C')
perm_bray <- usedist::dist_subset(bray, perm_meta$Sample)
adonis2(perm_bray ~ perm_meta$sex)

# adding sample number to beta diversity of Poly I:C samples
p_bray <- points_bray %>%
  filter(treatment == 'Poly I:C') %>%
  ggplot(aes(x = X1, y = X2)) +
  geom_label(aes(colour = sex, label = description), fill = NA) +
  theme_bw() +
  scale_colour_manual(values = sexColours) +
  labs(x = paste0('PC1 [', round(eigs[1], 2), '% Var. Explained]'),
       y = paste0('PC2 [', round(eigs[2], 2), '% Var. Explained]'))

p_bray %>% saveRDS('bray_polyIC_sampleID.RDS')
pdf('bray_polyIC_sampleID.pdf', width = 6, height = 4)
p_bray
dev.off()

# plotting top 20 species abundances
topspecies <- meanRA %>%
  slice_max(order_by = MeanRA, n = 20)

p <- df_long %>%
  mutate(Species = if_else(Species %in% topspecies$Species, Species, 'Other')) %>%
  group_by(Sample, Species) %>%
  summarise(RA = sum(RA)) %>%
  mutate(RA = 100*RA) %>%
  left_join(metadata, by = 'Sample') %>%
  ggplot(aes(x = Sample, y = RA)) +
  geom_bar(aes(fill = Species), colour = 'black', stat = 'identity') +
  ggforce::facet_row(~treatment+sex, scales = 'free_x', space = 'free') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = 'bottom') +
  scale_fill_manual(values = hues::iwanthue(21)) +
  labs(x = 'Relative Abundance (%)')

pdf('barplots.pdf', width = 15, height = 4)
p
dev.off()

p %>% saveRDS('barplots.RDS')

# plotting BC distance between samples
df <- bray %>%
  as.matrix() %>%
  data.frame() %>%
  rownames_to_column('Sample1') %>%
  pivot_longer(-Sample1, names_to = 'Sample2', values_to = 'Dist') %>%
  left_join(metadata, by = c('Sample1' = 'Sample')) %>%
  left_join(metadata, by = c('Sample2' = 'Sample')) %>%
  filter(Sample1 != Sample2)

pA <- df %>%
  filter(treatment.x == 'Saline' & treatment.y == 'Poly I:C') %>%
  filter(sex.x == sex.y) %>%
  ggplot(aes(x = sex.x, y = Dist)) +
  geom_boxplot(aes(colour = sex.x), fill = NA, outlier.shape = NA, show.legend = F) +
  geom_point(aes(colour = sex.x), size = 3, 
             position = position_jitterdodge(jitter.width = 0.5)) +
  theme_bw() +
  scale_colour_manual(values = sexColours)

pB <- df %>%
  filter(sex.x == 'm' & sex.y == 'f') %>%
  filter(treatment.x == treatment.y) %>%
  ggplot(aes(x = treatment.x, y = Dist)) +
  geom_boxplot(aes(colour = treatment.x), fill = NA, outlier.shape = NA, show.legend = F) +
  geom_point(aes(colour = treatment.x), size = 3, 
             position = position_jitterdodge(jitter.width = 0.5)) +
  theme_bw() +
  scale_colour_manual(values = treatmentColours)

p <- cowplot::plot_grid(plotlist = list(pA, pB), nrow = 1, labels = LETTERS[1:2])
p %>% saveRDS('braydist_combined.RDS')
pdf('braydist_combined.pdf', width = 8, height = 4)
p
dev.off()

df %>%
  filter(treatment.x == 'Saline' & treatment.y == 'Poly I:C') %>%
  filter(sex.x == sex.y) %>%
  t_test(Dist ~ sex.x)

df %>%
  filter(sex.x == 'm' & sex.y == 'f') %>%
  filter(treatment.x == treatment.y) %>%
  t_test(Dist ~ treatment.x)




