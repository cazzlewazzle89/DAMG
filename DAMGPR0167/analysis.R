library('tidyverse')
library('rstatix')
library('vegan')
library('lme4')

# importing metadata
metadata <- readxl::read_xlsx('metadata.xlsx') 

# plotting read count
reads <- readxl::read_xlsx('numreads.xlsx') %>%
  left_join(metadata)

reads %>%
  writexl::write_xlsx('readcounts.xlsx')

# importing Emu profiles
df <- read_tsv('emu-combined-species.tsv')

df <- df %>%
  select(species, starts_with('DMG')) %>%
  rename(Species = species)
df[is.na(df)] <- 0

# removing a sample with low read counts
df <- df %>%
  select(-DMG2303277)
metadata <- metadata %>%
  filter(Sample != 'DMG2303277')

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
richness <- specnumber(df_wide)
bray <- vegdist(df_wide, method = 'bray')
bray <- usedist::dist_subset(bray, metadata$Sample)
pcoa_bray <- cmdscale(bray, k = 2, eig = T)
points_bray <- pcoa_bray$points %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(metadata)
eigs <- 100*(pcoa_bray$eig/sum(pcoa_bray$eig))

# defining colour palettes
colours <- hues::iwanthue(4)
treatmentColours <- colours[1:2]
names(treatmentColours) <- unique(metadata$Treatment)
sexColours <- colours[3:4]
names(sexColours) <- unique(metadata$Sex)

shannon_df <- shannon %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  rename('Shannon' = 2) %>%
  left_join(metadata)

# all data together
p_shannon <- shannon_df %>%
  ggplot(aes(x = Treatment, y = Shannon)) +
  geom_boxplot(aes(colour = Treatment), fill = NA, outlier.shape = NA, show.legend = F) +
  geom_point(aes(colour = Treatment, shape = Sex), size = 3, 
             position = position_dodge(width = 0.25)) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(shape = 23))) +
  scale_shape_manual(values = list(21, 22)) +
  scale_color_manual(values = treatmentColours)

p_shannon %>% saveRDS('shannon_all.RDS')
pdf('shannon_all.pdf', width = 5, height = 4)
p_shannon
dev.off()

p_shannon <- shannon_df %>%
  ggplot(aes(x = Treatment, y = Shannon)) +
  geom_boxplot(aes(colour = Cohort), fill = NA, outlier.shape = NA, show.legend = F) +
  geom_point(aes(colour = Cohort, shape = Sex), size = 3, 
             position = position_dodge(width = 0.25)) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(shape = 23))) +
  scale_shape_manual(values = list(21, 22))

p_shannon %>% saveRDS('shannon_all_cohort.RDS')
pdf('shannon_all_cohort.pdf', width = 5, height = 4)
p_shannon
dev.off()

shannon_df %>%
  t_test(Shannon ~ Treatment)

shannon_df %>%
  anova_test(Shannon ~ Treatment + Sex)
shannon_df %>%
  anova_test(Shannon ~ Treatment * Sex)

shannon_df %>%
  anova_test(Shannon ~ Treatment + Cohort)
shannon_df %>%
  anova_test(Shannon ~ Treatment * Cohort)

shannon_df %>%
  anova_test(Shannon ~ Treatment + Sex + Cohort)
shannon_df %>%
  anova_test(Shannon ~ Treatment * Sex * Cohort)

shannon_lm <- lm(Shannon ~ Treatment + Cohort, data = shannon_df)
summary(shannon_lm)

richness_df <- richness %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  rename('Richness' = 2) %>%
  left_join(metadata)

p_richness <- richness_df %>%
  ggplot(aes(x = Treatment, y = Richness)) +
  geom_boxplot(aes(colour = Treatment), fill = NA, outlier.shape = NA, show.legend = F) +
  geom_point(aes(colour = Treatment, shape = Sex), size = 3, 
             position = position_dodge(width = 0.25)) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(shape = 23))) +
  scale_shape_manual(values = list(21, 22)) +
  scale_color_manual(values = treatmentColours)

p_richness %>% saveRDS('richness_all.RDS')
pdf('richness_all.pdf', width = 5, height = 4)
p_richness
dev.off()

p_richness <- richness_df %>%
  ggplot(aes(x = Treatment, y = Richness)) +
  geom_boxplot(aes(colour = Cohort), fill = NA, outlier.shape = NA, show.legend = F) +
  geom_point(aes(colour = Cohort, shape = Sex), size = 3, 
             position = position_dodge(width = 0.25)) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(shape = 23))) +
  scale_shape_manual(values = list(21, 22))

p_richness %>% saveRDS('richness_all_cohort.RDS')
pdf('richness_all_cohort.pdf', width = 5, height = 4)
p_richness
dev.off()

richness_df %>%
  t_test(Richness ~ Treatment)

richness_df %>%
  anova_test(Richness ~ Treatment + Sex)
richness_df %>%
  anova_test(Richness ~ Treatment * Sex)

richness_df %>%
  anova_test(Richness ~ Treatment + Cohort)
richness_df %>%
  anova_test(Richness ~ Treatment * Cohort)

richness_df %>%
  anova_test(Richness ~ Treatment + Sex + Cohort)
richness_df %>%
  anova_test(Richness ~ Treatment * Sex * Cohort)

richness_lm <- lm(Richness ~ Treatment + Cohort, data = richness_df)
summary(richness_lm)

# looking at only Cohort1 data
shannon_cohort1_df <- shannon_df %>%
  filter(Cohort == 'Cohort1')
richness_cohort1_df <- richness_df %>%
  filter(Cohort == 'Cohort1')

p_shannon <- shannon_cohort1_df %>%
  ggplot(aes(x = Treatment, y = Shannon)) +
  geom_boxplot(aes(colour = Treatment), fill = NA, outlier.shape = NA, show.legend = F) +
  geom_point(aes(colour = Treatment, shape = Sex), size = 3, 
             position = position_dodge(width = 0.25)) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(shape = 23))) +
  scale_shape_manual(values = list(21, 22)) +
  scale_color_manual(values = treatmentColours)

p_shannon %>% saveRDS('shannon_cohort1.RDS')
pdf('shannon_cohort1.pdf', width = 5, height = 4)
p_shannon
dev.off()

p_richness <- richness_cohort1_df %>%
  ggplot(aes(x = Treatment, y = Richness)) +
  geom_boxplot(aes(colour = Treatment), fill = NA, outlier.shape = NA, show.legend = F) +
  geom_point(aes(colour = Treatment, shape = Sex), size = 3, 
             position = position_dodge(width = 0.25)) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(shape = 23))) +
  scale_shape_manual(values = list(21, 22)) +
  scale_color_manual(values = treatmentColours)

p_richness %>% saveRDS('richness_cohort1.RDS')
pdf('richness_cohort1.pdf', width = 5, height = 4)
p_richness
dev.off()

shannon_cohort1_df %>%
  anova_test(Shannon ~ Treatment*Sex)

richness_cohort1_df %>%
  anova_test(Richness ~ Treatment*Sex)

# looking at only Cohort3 data
shannon_cohort3_df <- shannon_df %>%
  filter(Cohort == 'Cohort3')
richness_cohort3_df <- richness_df %>%
  filter(Cohort == 'Cohort3')

p_shannon <- shannon_cohort3_df %>%
  ggplot(aes(x = Treatment, y = Shannon)) +
  geom_boxplot(aes(colour = Treatment), fill = NA, outlier.shape = NA, show.legend = F) +
  geom_point(aes(colour = Treatment, shape = Sex), size = 3, 
             position = position_dodge(width = 0.25)) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(shape = 23))) +
  scale_shape_manual(values = list(21, 22)) +
  scale_color_manual(values = treatmentColours)

p_shannon %>% saveRDS('shannon_cohort3.RDS')
pdf('shannon_cohort3.pdf', width = 5, height = 4)
p_shannon
dev.off()

p_richness <- richness_cohort3_df %>%
  ggplot(aes(x = Treatment, y = Richness)) +
  geom_boxplot(aes(colour = Treatment), fill = NA, outlier.shape = NA, show.legend = F) +
  geom_point(aes(colour = Treatment, shape = Sex), size = 3, 
             position = position_dodge(width = 0.25)) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(shape = 23))) +
  scale_shape_manual(values = list(21, 22)) +
  scale_color_manual(values = treatmentColours)

p_richness %>% saveRDS('richness_cohort3.RDS')
pdf('richness_cohort3.pdf', width = 5, height = 4)
p_richness
dev.off()

shannon_cohort3_df %>%
  anova_test(Shannon ~ Treatment*Sex)

richness_cohort3_df %>%
  anova_test(Richness ~ Treatment*Sex)

# combined alpha diversity plot of Cohort3 data
alpha_cohort3 <- bind_rows(shannon_cohort3_df %>%
                             mutate(Measure = 'Shannon') %>%
                             rename(Diversity = Shannon),
                           richness_cohort3_df %>%
                             mutate(Measure = 'Richness') %>%
                             rename(Diversity = Richness))
p_alpha <- alpha_cohort3 %>%
  ggplot(aes(x = Treatment, y = Diversity)) +
  ggforce::facet_row(~Measure, scales = 'free') +
  geom_boxplot(aes(colour = Treatment), fill = NA, outlier.shape = NA, show.legend = F) +
  geom_point(aes(colour = Treatment, shape = Sex), size = 3, 
             position = position_dodge(width = 0.25)) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(shape = 23))) +
  scale_shape_manual(values = list(21, 22)) +
  scale_color_manual(values = treatmentColours)

p_alpha %>% saveRDS('alpha_cohort3.RDS')
pdf('alpha_cohort3.pdf', width = 8, height = 4)
p_alpha
dev.off()

# beta diversity using all samples
adonis2(bray ~ Treatment, data = metadata)
adonis2(bray ~ Treatment + Sex, data = metadata)
adonis2(bray ~ Treatment + Sex + Cohort, data = metadata)
adonis2(bray ~ Treatment * Sex * Cohort, data = metadata)

# beta diversity separated by Cohort
metadata_cohort1 <- metadata %>% filter(Cohort == 'Cohort1')
metadata_cohort3 <- metadata %>% filter(Cohort == 'Cohort3')
bray_cohort1 <- usedist::dist_subset(bray, metadata_cohort1$Sample)
bray_cohort3 <- usedist::dist_subset(bray, metadata_cohort3$Sample)

adonis2(bray_cohort1 ~ Treatment, data = metadata_cohort1)
adonis2(bray_cohort1 ~ Treatment * Sex, data = metadata_cohort1)
adonis2(bray_cohort3 ~ Treatment, data = metadata_cohort3)
adonis2(bray_cohort3 ~ Treatment * Sex, data = metadata_cohort3)

pcoa_bray_cohort1 <- cmdscale(bray_cohort1, k = 2, eig = T)
points_bray_cohort1 <- pcoa_bray_cohort1$points %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(metadata)
eigs_cohort1 <- 100*(pcoa_bray_cohort1$eig/sum(pcoa_bray_cohort1$eig))

pcoa_bray_cohort3 <- cmdscale(bray_cohort3, k = 2, eig = T)
points_bray_cohort3 <- pcoa_bray_cohort3$points %>%
  data.frame() %>%
  rownames_to_column('Sample') %>%
  left_join(metadata)
eigs_cohort3 <- 100*(pcoa_bray_cohort3$eig/sum(pcoa_bray_cohort3$eig))

p <- points_bray_cohort1 %>%
  ggplot(aes(x = X1, y = X2)) +
  geom_point(aes(fill = Treatment, shape = Sex), colour = 'black', size = 3) +
  theme_bw() +
  scale_shape_manual(values = c(21, 23)) +
  guides(fill = guide_legend(override.aes = list(shape = 22))) +
  labs(x = paste0('PC1 [', round(eigs_cohort1[1], 2), '%]'),
       y = paste0('PC2 [', round(eigs_cohort1[2], 2), '%]'))
pdf('bray_cohort1.pdf', width = 5, height = 3)
p
dev.off()

p <- points_bray_cohort3 %>%
  ggplot(aes(x = X1, y = X2)) +
  geom_point(aes(fill = Treatment, shape = Sex), colour = 'black', size = 3) +
  theme_bw() +
  scale_shape_manual(values = c(21, 23)) +
  guides(fill = guide_legend(override.aes = list(shape = 22))) +
  labs(x = paste0('PC1 [', round(eigs_cohort3[1], 2), '%]'),
       y = paste0('PC2 [', round(eigs_cohort3[2], 2), '%]'))
pdf('bray_cohort3.pdf', width = 5, height = 3)
p
dev.off()

saveRDS(p, 'bray_cohort3.RDS')

# barplots of the top 10 most abundant species
topspecies <- meanRA %>%
  slice_max(order_by = MeanRA, n = 10)
p <- df_long %>%
  mutate(Species = if_else(Species %in% topspecies$Species, Species, 'Other')) %>%
  group_by(Sample, Species) %>%
  summarise(RA= sum(RA)) %>%
  left_join(metadata, by = 'Sample') %>%
  ggplot(aes(x = Sample, y = RA)) +
  geom_bar(aes(fill = Species), colour = 'black', stat = 'identity') +
  ggforce::facet_row(~Cohort+Treatment, scales = 'free_x') +
  theme_bw() +
  scale_fill_manual(values = hues::iwanthue(11)) +
  labs(x = '',
       y = 'Relative Abundance') +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

pdf('barplots.pdf', width = 12, height = 4)
p
dev.off()

saveRDS(p, 'barplots.RDS')

# cohort 3 looks more homogeneous than cohort1
# possible effect of being in the freezer?
# will check with pairwise bray distances within cohorts
bray_pairwise <- bray %>%
  as.matrix() %>%
  data.frame() %>%
  rownames_to_column('Sample1') %>%
  pivot_longer(cols = -Sample1, names_to = 'Sample2', values_to = 'Dist') %>%
  left_join(metadata, by = c('Sample1' = 'Sample')) %>%
  left_join(metadata, by = c('Sample2' = 'Sample')) %>%
  filter(Cohort.x == Cohort.y & Treatment.x == Treatment.y & Sample1 < Sample2) 

bray_pairwise %>%
  anova_test(Dist ~ Cohort.x*Treatment.x)

p <- bray_pairwise %>%
  ggplot(aes(x = Treatment.x, y = Dist)) +
  geom_boxplot(aes(colour = Cohort.x), fill = NA, outlier.shape = NA) +
  geom_point(aes(fill = Cohort.x), shape = 21, colour = 'black', size = 3, position = position_jitter(width = 0.2)) +
  theme_bw() +
  labs(x = '', y = 'Pairwise Bray-Curtis Distance')

pdf('pairwise_dist.pdf', width = 6, height = 4)
p
dev.off()

saveRDS(p, 'pairwise_dist.RDS')




  
  
  
  
  
