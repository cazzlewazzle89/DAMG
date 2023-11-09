
# load libraries
{
  library('ALDEx2')
  library('qiime2R')
  library('rstatix')
  library('vegan')
  library('nlme')
  library('emmeans')
  library('tidyverse')
}

# import metadata and define group colours
{
  metadata <- read_delim('metadata.txt') %>%
    unite(col = Group_Time, c('group', 'timepoint'), sep = '_', remove = F)
  
  groupColours <- c('Chow_T1' = '#d30000',
                    'Chow_T2' = '#fa8072',
                    'IA_T1' = '#3944bc',
                    'IA_T2' = '#63c5da',
                    'Ad lib_T1' = '#028a0f',
                    'Ad lib_T2' = '#98bf64')
}

# import and combine alpha diversity metrics
{
  alpha_measures <- c('evenness', 'faith_pd', 'observed_features', 'shannon')
  alpha_df <- data.frame()
  
  for (i in alpha_measures){
    df <- read_qza(paste0('CoreMetricsPhylogenetic/', i, '_vector.qza'))$data %>%
      rename('Diversity' = 1) %>%
      mutate(Measure = i) %>%
      rownames_to_column('sample_name')
  
    alpha_df <- rbind(alpha_df, df)
  }
  
  alpha_df <- alpha_df %>%
    left_join(metadata)

  alpha_df %>%
    pivot_wider(id_cols = c('sample_name', 'description', 'group', 'timepoint', 'Group_Time', 'mouse'),
                names_from = Measure, values_from = Diversity, values_fill = NA) %>%
    write.table('alpha_diversity.tsv', quote = F, row.names = F, sep = '\t')
}

# test for differences in shannon diversity between groups and timepoints 
{ 
  test_in <- alpha_df %>%
    filter(Measure == 'shannon') %>%
    arrange(sample_name)
  
  shapiro_res <- test_in %>%
    group_by(group, timepoint) %>%
    shapiro_test(Diversity)
  
  outliers_res <- test_in %>%
    group_by(group, timepoint) %>%
    identify_outliers(Diversity)
  
  levene_res <- test_in %>%
    group_by(timepoint) %>%
    levene_test(Diversity ~ group)
  
  box_res <- box_m(test_in[, "Diversity", drop = FALSE], test_in$group)
  
  anova_res <- anova_test(data = test_in,
                          dv = Diversity,
                          wid = mouse,
                          between = group,
                          within = timepoint)
  
  get_anova_table(anova_res) %>%
    data.frame() %>%
    dplyr::select(-ges) %>%
    rename('p.sig' = 'p..05') %>%
    mutate(p.sig = if_else(p <= 0.05, '*', 'ns')) %>%
    write_delim('alpha_anova_res.tsv', delim = '\t')
}

# test for differences in Faith's Phylogenetic Diversity between groups and timepoints 
{ 
  test_in <- alpha_df %>%
    filter(Measure == 'faith_pd') %>%
    arrange(sample_name)
  
  anova_res <- anova_test(data = test_in,
                          dv = Diversity,
                          wid = mouse,
                          between = group,
                          within = timepoint)
  
  get_anova_table(anova_res) %>%
    data.frame() %>%
    dplyr::select(-ges) %>%
    rename('p.sig' = 'p..05') %>%
    mutate(p.sig = if_else(p <= 0.05, '*', 'ns')) 
  
  t_test(data = test_in %>% filter(group == 'Chow'),
         formula = Diversity ~ timepoint)
  t_test(data = test_in %>% filter(group == 'Ad lib'),
         formula = Diversity ~ timepoint)
  t_test(data = test_in %>% filter(group == 'IA'),
         formula = Diversity ~ timepoint)
  
  test_in %>%
    group_by(group, timepoint) %>%
    summarise(Mean = mean(Diversity),
              Median = median(Diversity))
}

# import and combine beta diversity metrics
{
  beta_measures <- c('bray_curtis', 'jaccard', 'weighted_unifrac', 'unweighted_unifrac')
  beta_df <- data.frame()
  for (i in beta_measures){
    df <- read_qza(paste0('CoreMetricsPhylogenetic/', i, '_pcoa_results.qza'))$data$Vectors %>%
      rename('sample_name' = 'SampleID') %>%
      mutate(Measure = i)
    
    beta_df <- rbind(beta_df, df)
  }
  
  beta_df <- beta_df %>%
    left_join(metadata)
  
  beta_df %>%
    write_delim('beta_coordinates.tsv', delim = '\t')
}

# test for differences in Bray-Curtis diversity between groups and timepoints
{
  bray_dist <- read_qza('CoreMetricsPhylogenetic/bray_curtis_distance_matrix.qza')$data
  
  bray_meta <- bray_dist %>%
    as.matrix() %>%
    data.frame() %>%
    rownames_to_column('sample_name') %>%
    dplyr::select(sample_name) %>%
    left_join(metadata)
  
  permout <- adonis(bray_dist ~ group * timepoint, data = bray_meta)
  
  data.frame(Variable = c('Group', 'Timepoint', 'Group:Timepoint'),
             R2 = permout$aov.tab$R2[1:3],
             p = permout$aov.tab$`Pr(>F)`[1:3]) %>%
    mutate(p.sig = if_else(p <= 0.05, '*', 'ns')) %>%
    write_delim('bray_permout.tsv', delim = '\t')
}

# test for differences in Unweighted Unifrac diversity between groups and timepoints
{
  unifrac_dist <- read_qza('CoreMetricsPhylogenetic/unweighted_unifrac_distance_matrix.qza')$data
  
  unifrac_meta <- unifrac_dist %>%
    as.matrix() %>%
    data.frame() %>%
    rownames_to_column('sample_name') %>%
    dplyr::select(sample_name) %>%
    left_join(metadata)
  
  permout <- adonis(unifrac_dist ~ group * timepoint, data = unifrac_meta)
  
  data.frame(Variable = c('Group', 'Timepoint', 'Group:Timepoint'),
             R2 = permout$aov.tab$R2[1:3],
             p = permout$aov.tab$`Pr(>F)`[1:3]) %>%
    mutate(p.sig = if_else(p <= 0.05, '*', 'ns')) %>%
    write_delim('unifrac_permout.tsv', delim = '\t')
}

# plot Unifrac PCoA
{
  a <- read_qza('CoreMetricsPhylogenetic/unweighted_unifrac_pcoa_results.qza')$data$ProportionExplained
  beta_df <- read_delim('beta_coordinates.tsv', delim = '\t')
  
  p <- beta_df %>%
    filter(Measure == 'unweighted_unifrac') %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point(aes(fill = Group_Time), shape = 21, colour = 'black', size = 3) +
    stat_ellipse(aes(colour = Group_Time), size = 1) +
    theme_bw() +
    labs(x = paste0('PC1 [', round(100*(a[1]), digits = 2), '%]'),
         y = paste0('PC2 [', round(100*(a[2]), digits = 2), '%]')) +
    scale_fill_manual(values = groupColours) +
    scale_colour_manual(values = groupColours)
  
  print(p)
}

# posthoc testing for differences in Bray-Curtis diversity between groups at each timepoint
{
  bray_dist_T1 <- usedist::dist_subset(bray_dist, bray_meta$timepoint == 'T1')
  bray_dist_T2 <- usedist::dist_subset(bray_dist, bray_meta$timepoint == 'T2')
  
  bray_meta_T1 <- bray_meta %>% filter(timepoint == 'T1')
  bray_meta_T2 <- bray_meta %>% filter(timepoint == 'T2')
  
  permout_T1 <- RVAideMemoire::pairwise.perm.manova(bray_dist_T1, bray_meta_T1$group, R2 = T)
  permout_T2 <- RVAideMemoire::pairwise.perm.manova(bray_dist_T2, bray_meta_T2$group, R2 = T)
  
  temp_T1p <- permout_T1$p.value %>%
    data.frame() %>%
    rownames_to_column('Group1') %>%
    pivot_longer(cols = -Group1, names_to = 'Group2', values_to = 'p', values_drop_na = T)
  temp_T1r <- permout_T1$R2.value %>%
    data.frame() %>%
    rownames_to_column('Group1') %>%
    pivot_longer(cols = -Group1, names_to = 'Group2', values_to = 'R2', values_drop_na = T)
  
  temp_T1 <- temp_T1r %>%
    left_join(temp_T1p,
              by = c('Group1', 'Group2')) %>%
    mutate(Timepoint = 'T1', .before = Group1)
  
  temp_T2p <- permout_T2$p.value %>%
    data.frame() %>%
    rownames_to_column('Group1') %>%
    pivot_longer(cols = -Group1, names_to = 'Group2', values_to = 'p', values_drop_na = T)
  temp_T2r <- permout_T2$R2.value %>%
    data.frame() %>%
    rownames_to_column('Group1') %>%
    pivot_longer(cols = -Group1, names_to = 'Group2', values_to = 'R2', values_drop_na = T)
  
  temp_T2 <- temp_T2r %>%
    left_join(temp_T2p,
              by = c('Group1', 'Group2')) %>%
    mutate(Timepoint = 'T2', .before = Group1)
  
  bind_rows(temp_T1, temp_T2) %>%
    mutate(p.sig = if_else(p <= 0.05, '*', 'ns')) %>%
    write_delim('bray_permout_groups_pairwise.tsv', delim = '\t')
}

# posthoc testing for differences in Unweighted Unifrac diversity between groups at each timepoint
{
  unifrac_dist_T1 <- usedist::dist_subset(unifrac_dist, unifrac_meta$timepoint == 'T1')
  unifrac_dist_T2 <- usedist::dist_subset(unifrac_dist, unifrac_meta$timepoint == 'T2')
  
  unifrac_meta_T1 <- unifrac_meta %>% filter(timepoint == 'T1')
  unifrac_meta_T2 <- unifrac_meta %>% filter(timepoint == 'T2')
  
  permout_T1 <- RVAideMemoire::pairwise.perm.manova(unifrac_dist_T1, unifrac_meta_T1$group, R2 = T)
  permout_T2 <- RVAideMemoire::pairwise.perm.manova(unifrac_dist_T2, unifrac_meta_T2$group, R2 = T)
  
  temp_T1p <- permout_T1$p.value %>%
    data.frame() %>%
    rownames_to_column('Group1') %>%
    pivot_longer(cols = -Group1, names_to = 'Group2', values_to = 'p', values_drop_na = T)
  temp_T1r <- permout_T1$R2.value %>%
    data.frame() %>%
    rownames_to_column('Group1') %>%
    pivot_longer(cols = -Group1, names_to = 'Group2', values_to = 'R2', values_drop_na = T)
  
  temp_T1 <- temp_T1r %>%
    left_join(temp_T1p,
              by = c('Group1', 'Group2')) %>%
    mutate(Timepoint = 'T1', .before = Group1)
  
  temp_T2p <- permout_T2$p.value %>%
    data.frame() %>%
    rownames_to_column('Group1') %>%
    pivot_longer(cols = -Group1, names_to = 'Group2', values_to = 'p', values_drop_na = T)
  temp_T2r <- permout_T2$R2.value %>%
    data.frame() %>%
    rownames_to_column('Group1') %>%
    pivot_longer(cols = -Group1, names_to = 'Group2', values_to = 'R2', values_drop_na = T)
  
  temp_T2 <- temp_T2r %>%
    left_join(temp_T2p,
              by = c('Group1', 'Group2')) %>%
    mutate(Timepoint = 'T2', .before = Group1)
  
  bind_rows(temp_T1, temp_T2) %>%
    mutate(p.sig = if_else(p <= 0.05, '*', 'ns')) %>%
    write_delim('unifrac_permout_groups_pairwise.tsv', delim = '\t')
}

# posthoc testing for differences in Bray-Curtis diversity between timepoints in each group
{
  bray_dist_A <- usedist::dist_subset(bray_dist, bray_meta$group == 'Chow')
  bray_dist_B <- usedist::dist_subset(bray_dist, bray_meta$group == 'IA')
  bray_dist_C <- usedist::dist_subset(bray_dist, bray_meta$group == 'Ad lib')
  
  bray_meta_A <- bray_meta %>% filter(group == 'Chow')
  bray_meta_B <- bray_meta %>% filter(group == 'IA')
  bray_meta_C <- bray_meta %>% filter(group == 'Ad lib')
  
  permout_A <- adonis(bray_dist_A ~ timepoint, data = bray_meta_A)
  permout_B <- adonis(bray_dist_B ~ timepoint, data = bray_meta_B)
  permout_C <- adonis(bray_dist_C ~ timepoint, data = bray_meta_C)
  
  temp_A <- data.frame(Group = 'Chow',
                       R2 = permout_A$aov.tab$R2[1],
                       p = permout_A$aov.tab$`Pr(>F)`[1])
  
  temp_B <- data.frame(Group = 'IA',
                       R2 = permout_B$aov.tab$R2[1],
                       p = permout_B$aov.tab$`Pr(>F)`[1])
  
  temp_C <- data.frame(Group = 'Ad lib',
                       R2 = permout_C$aov.tab$R2[1],
                       p = permout_C$aov.tab$`Pr(>F)`[1])
  
  bind_rows(temp_A, temp_B, temp_C) %>%
    mutate(p.sig = if_else(p <= 0.05, '*', 'ns')) %>%
    write_delim('bray_permout_time_pairwise.tsv', delim = '\t')
}

# posthoc testing for differences in Unweighted Unifrac diversity between timepoints in each group
{
  unifrac_dist_A <- usedist::dist_subset(unifrac_dist, unifrac_meta$group == 'Chow')
  unifrac_dist_B <- usedist::dist_subset(unifrac_dist, unifrac_meta$group == 'IA')
  unifrac_dist_C <- usedist::dist_subset(unifrac_dist, unifrac_meta$group == 'Ad lib')
  
  unifrac_meta_A <- unifrac_meta %>% filter(group == 'Chow')
  unifrac_meta_B <- unifrac_meta %>% filter(group == 'IA')
  unifrac_meta_C <- unifrac_meta %>% filter(group == 'Ad lib')
  
  permout_A <- adonis(unifrac_dist_A ~ timepoint, data = unifrac_meta_A)
  permout_B <- adonis(unifrac_dist_B ~ timepoint, data = unifrac_meta_B)
  permout_C <- adonis(unifrac_dist_C ~ timepoint, data = unifrac_meta_C)
  
  temp_A <- data.frame(Group = 'Chow',
                       R2 = permout_A$aov.tab$R2[1],
                       p = permout_A$aov.tab$`Pr(>F)`[1])
  
  temp_B <- data.frame(Group = 'IA',
                       R2 = permout_B$aov.tab$R2[1],
                       p = permout_B$aov.tab$`Pr(>F)`[1])
  
  temp_C <- data.frame(Group = 'Ad lib',
                       R2 = permout_C$aov.tab$R2[1],
                       p = permout_C$aov.tab$`Pr(>F)`[1])
  
  bind_rows(temp_A, temp_B, temp_C) %>%
    mutate(p.sig = if_else(p <= 0.05, '*', 'ns')) %>%
    write_delim('unifrac_permout_time_pairwise.tsv', delim = '\t')
}

# import feature table 
# filter to only include sOTUs that are > 0.1% in at least one group-timepoint combo
{
  feature_table <- phyloseq::import_biom('feature-table_json.biom')
  taxonomy <- read_delim('taxonomy.tsv', delim = '\t') %>%
    dplyr::select(-Confidence)
  
  otus_to_keep <- feature_table %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    dplyr::select(OTU, Sample, RA) %>%
    filter(RA >= 0.1) %>%
    pull(OTU) %>%
    unique()
  
  feature_table_filt <- feature_table %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    filter(OTU %in% otus_to_keep) %>%
    column_to_rownames('OTU')
}

# identify OTUs that are differentially abundant between timepoints within each group
{
  samples_A <- metadata %>%
    filter(group == 'Chow')
  samples_B <- metadata %>%
    filter(group == 'IA')
  samples_C <- metadata %>%
    filter(group == 'Ad lib')
  
  aldex_in_A <- feature_table_filt %>%
    dplyr::select(samples_A$sample_name)
  aldex_in_B <- feature_table_filt %>%
    dplyr::select(samples_B$sample_name)
  aldex_in_C <- feature_table_filt %>%
    dplyr::select(samples_C$sample_name)
  
  aldex_out_A <- aldex(aldex_in_A,
                       samples_A$timepoint,
                       mc.samples = 128,
                       test = "t",
                       paired.test = TRUE,
                       effect = TRUE,
                       include.sample.summary = FALSE, 
                       denom = "all",
                       verbose = FALSE) %>%
    rownames_to_column('OTU') %>%
    left_join(taxonomy,
              by = c('OTU' = 'Feature ID'))
  
  aldex_out_B <- aldex(aldex_in_B,
                       samples_B$timepoint,
                       mc.samples = 128,
                       test = "t",
                       paired.test = TRUE,
                       effect = TRUE,
                       include.sample.summary = FALSE, 
                       denom = "all",
                       verbose = FALSE) %>%
    rownames_to_column('OTU') %>%
    left_join(taxonomy,
              by = c('OTU' = 'Feature ID'))
  
  aldex_out_C <- aldex(aldex_in_C,
                       samples_C$timepoint,
                       mc.samples = 128,
                       test = "t",
                       paired.test = TRUE,
                       effect = TRUE,
                       include.sample.summary = FALSE, 
                       denom = "all",
                       verbose = FALSE) %>%
    rownames_to_column('OTU') %>%
    left_join(taxonomy,
              by = c('OTU' = 'Feature ID'))
  
  list(GroupA = aldex_out_A,
       GroupB = aldex_out_B,
       GroupC = aldex_out_C) %>%
    writexl::write_xlsx('aldex2_timepoint_all.xlsx')
  
  list(GroupA = aldex_out_A %>%
         filter(wi.eBH <= 0.05),
       GroupB = aldex_out_B %>%
         filter(wi.eBH <= 0.05),
       GroupC = aldex_out_C %>%
         filter(wi.eBH <= 0.05)) %>%
    writexl::write_xlsx('aldex2_timepoint_sig.xlsx')
  
  sig_otus_A <- aldex_out_A %>%
    filter(wi.eBH <= 0.05) %>%
    pull(OTU)
  sig_otus_B <- aldex_out_B %>%
    filter(wi.eBH <= 0.05) %>%
    pull(OTU)
  sig_otus_C <- aldex_out_C %>%
    filter(wi.eBH <= 0.05) %>%
    pull(OTU)
  
  ra_summary_A <- feature_table_filt %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(OTU %in% sig_otus_A) %>%
    dplyr::select(OTU, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(group == 'Chow') %>%
    group_by(OTU, timepoint) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'OTU', names_from = 'timepoint', names_prefix = 'RA_', values_from = 'MeanRA') %>%
    left_join(taxonomy,
              by = c('OTU' = 'Feature ID'))
  
  ra_summary_B <- feature_table_filt %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(OTU %in% sig_otus_B) %>%
    dplyr::select(OTU, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(group == 'IA') %>%
    group_by(OTU, timepoint) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'OTU', names_from = 'timepoint', names_prefix = 'RA_', values_from = 'MeanRA') %>%
    left_join(taxonomy,
              by = c('OTU' = 'Feature ID'))
  
  ra_summary_C <- feature_table_filt %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(OTU %in% sig_otus_C) %>%
    dplyr::select(OTU, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(group == 'Ad lib') %>%
    group_by(OTU, timepoint) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'OTU', names_from = 'timepoint', names_prefix = 'RA_', values_from = 'MeanRA') %>%
    left_join(taxonomy,
              by = c('OTU' = 'Feature ID'))
  
  list(GroupA = ra_summary_A,
       GroupB = ra_summary_B,
       GroupC = ra_summary_C) %>%
    writexl::write_xlsx('aldex2_timepoint_sig_rasummary.xlsx')
}

# identify OTUs that are differentially abundant between groups within each timepoint
{
  samples_T1AB <- metadata %>%
    filter(group %in% c('Chow', 'IA') & timepoint == 'T1')
  samples_T1AC <- metadata %>%
    filter(group %in% c('Chow', 'Ad lib') & timepoint == 'T1')
  samples_T1BC <- metadata %>%
    filter(group %in% c('IA', 'Ad lib') & timepoint == 'T1')
  samples_T2AB <- metadata %>%
    filter(group %in% c('Chow', 'IA') & timepoint == 'T2')
  samples_T2AC <- metadata %>%
    filter(group %in% c('Chow', 'Ad lib') & timepoint == 'T2')
  samples_T2BC <- metadata %>%
    filter(group %in% c('IA', 'Ad lib') & timepoint == 'T2')
  
  aldex_in_T1AB <- feature_table_filt %>%
    dplyr::select(samples_T1AB$sample_name)
  aldex_in_T1AC <- feature_table_filt %>%
    dplyr::select(samples_T1AC$sample_name)
  aldex_in_T1BC <- feature_table_filt %>%
    dplyr::select(samples_T1BC$sample_name)
  aldex_in_T2AB <- feature_table_filt %>%
    dplyr::select(samples_T1AB$sample_name)
  aldex_in_T2AC <- feature_table_filt %>%
    dplyr::select(samples_T2AC$sample_name)
  aldex_in_T2BC <- feature_table_filt %>%
    dplyr::select(samples_T2BC$sample_name)
  
  aldex_out_T1AB <- aldex(aldex_in_T1AB, samples_T1AB$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>%
    rownames_to_column('OTU') %>% left_join(taxonomy, by = c('OTU' = 'Feature ID'))
  aldex_out_T1AC <- aldex(aldex_in_T1AC, samples_T1AC$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>%
    rownames_to_column('OTU') %>% left_join(taxonomy, by = c('OTU' = 'Feature ID'))
  aldex_out_T1BC <- aldex(aldex_in_T1BC, samples_T1BC$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>%
    rownames_to_column('OTU') %>% left_join(taxonomy, by = c('OTU' = 'Feature ID'))
  aldex_out_T2AB <- aldex(aldex_in_T2AB, samples_T2AB$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>%
    rownames_to_column('OTU') %>% left_join(taxonomy, by = c('OTU' = 'Feature ID'))
  aldex_out_T2AC <- aldex(aldex_in_T2AC, samples_T2AC$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>%
    rownames_to_column('OTU') %>% left_join(taxonomy, by = c('OTU' = 'Feature ID'))
  aldex_out_T2BC <- aldex(aldex_in_T2BC, samples_T2BC$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>%
    rownames_to_column('OTU') %>% left_join(taxonomy, by = c('OTU' = 'Feature ID'))
  
  list(T1AB = aldex_out_T1AB,
       T1AC = aldex_out_T1AC,
       T1BC = aldex_out_T1BC,
       T2AB = aldex_out_T2AB,
       T2AC = aldex_out_T2AC,
       T2BC = aldex_out_T2BC) %>%
    writexl::write_xlsx('aldex2_group_all.xlsx')
  
  list(T1AB = aldex_out_T1AB %>%
         filter(wi.eBH <= 0.05),
       T1AC = aldex_out_T1AC %>%
         filter(wi.eBH <= 0.05),
       T1BC = aldex_out_T1BC %>%
         filter(wi.eBH <= 0.05),
       T2AB = aldex_out_T2AB %>%
         filter(wi.eBH <= 0.05),
       T2AC = aldex_out_T2AC %>%
         filter(wi.eBH <= 0.05),
       T2BC = aldex_out_T2BC %>%
         filter(wi.eBH <= 0.05)) %>%
    writexl::write_xlsx('aldex2_group_sig.xlsx')
  
  sig_otus_T1AB <- aldex_out_T1AB %>%
    filter(wi.eBH <= 0.05) %>%
    pull(OTU)
  sig_otus_T1AC <- aldex_out_T1AC %>%
    filter(wi.eBH <= 0.05) %>%
    pull(OTU)
  sig_otus_T1BC <- aldex_out_T1BC %>%
    filter(wi.eBH <= 0.05) %>%
    pull(OTU)
  sig_otus_T2AB <- aldex_out_T2AB %>%
    filter(wi.eBH <= 0.05) %>%
    pull(OTU)
  sig_otus_T2AC <- aldex_out_T2AC %>%
    filter(wi.eBH <= 0.05) %>%
    pull(OTU)
  sig_otus_T2BC <- aldex_out_T2BC %>%
    filter(wi.eBH <= 0.05) %>%
    pull(OTU)
  
  ra_summary_T1AB <- feature_table_filt %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(OTU %in% sig_otus_T1AB) %>%
    dplyr::select(OTU, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T1' & group %in% c('Chow', 'IA')) %>%
    group_by(OTU, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'OTU', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') %>%
    left_join(taxonomy,
              by = c('OTU' = 'Feature ID'))
  
  ra_summary_T1AC <- feature_table_filt %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(OTU %in% sig_otus_T1AC) %>%
    dplyr::select(OTU, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T1' & group %in% c('Chow', 'Ad lib')) %>%
    group_by(OTU, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'OTU', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') %>%
    left_join(taxonomy,
              by = c('OTU' = 'Feature ID'))
  
  ra_summary_T1BC <- feature_table_filt %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(OTU %in% sig_otus_T1BC) %>%
    dplyr::select(OTU, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T1' & group %in% c('IA', 'Ad lib')) %>%
    group_by(OTU, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'OTU', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') %>%
    left_join(taxonomy,
              by = c('OTU' = 'Feature ID'))
  
  ra_summary_T2AB <- feature_table_filt %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(OTU %in% sig_otus_T2AB) %>%
    dplyr::select(OTU, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T2' & group %in% c('Chow', 'IA')) %>%
    group_by(OTU, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'OTU', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') %>%
    left_join(taxonomy,
              by = c('OTU' = 'Feature ID'))
  
  ra_summary_T2AC <- feature_table_filt %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(OTU %in% sig_otus_T2AC) %>%
    dplyr::select(OTU, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T2' & group %in% c('Chow', 'Ad lib')) %>%
    group_by(OTU, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'OTU', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') %>%
    left_join(taxonomy,
              by = c('OTU' = 'Feature ID'))
  
  ra_summary_T2BC <- feature_table_filt %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(OTU %in% sig_otus_T2BC) %>%
    dplyr::select(OTU, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T2' & group %in% c('IA', 'Ad lib')) %>%
    group_by(OTU, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'OTU', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') %>%
    left_join(taxonomy,
              by = c('OTU' = 'Feature ID'))
  
  list(T1AB = ra_summary_T1AB,
       T1AC = ra_summary_T1AB,
       T1BC = ra_summary_T1BC,
       T2AB = ra_summary_T2AB,
       T2AC = ra_summary_T2AC,
       T2BC = ra_summary_T2BC) %>%
    writexl::write_xlsx('aldex2_group_sig_rasummary.xlsx')
}

# collapse feature table to Genus level
{
  genus_table <- feature_table %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    left_join(taxonomy, by = c('OTU' = 'Feature ID')) %>%
    separate(Taxon, 
             into = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), 
             sep = '; ') %>%
    pivot_longer(cols = starts_with('DMG'), names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample, Genus) %>%
    summarise(Count = sum(Count)) %>%
    pivot_wider(id_cols = Genus, names_from = 'Sample', values_from = 'Count', values_fill = 0) %>%
    mutate(Genus = if_else(is.na(Genus), 'Unassigned', Genus)) %>%
    column_to_rownames('Genus')
}

# filter genus table to only include genera that are > 0.1% in at least one group-timepoint combo
{
  genera_to_keep <- genus_table %>%
    data.frame() %>%
    rownames_to_column('Genus') %>%
    pivot_longer(cols = -Genus, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    dplyr::select(Genus, Sample, RA) %>%
    filter(RA >= 0.1) %>%
    filter(!Genus %in% c('g__', 'Unassigned')) %>%
    pull(Genus) %>%
    unique()
  
  genus_table_filt <- genus_table %>%
    data.frame() %>%
    rownames_to_column('Genus') %>%
    filter(Genus %in% genera_to_keep) %>%
    column_to_rownames('Genus')
}

# identify genera that are differentially abundant between timepoints within each group
{
  samples_A <- metadata %>%
    filter(group == 'Chow')
  samples_B <- metadata %>%
    filter(group == 'IA')
  samples_C <- metadata %>%
    filter(group == 'Ad lib')
  
  aldex_in_A <- genus_table_filt %>%
    dplyr::select(samples_A$sample_name)
  aldex_in_B <- genus_table_filt %>%
    dplyr::select(samples_B$sample_name)
  aldex_in_C <- genus_table_filt %>%
    dplyr::select(samples_C$sample_name)
  
  aldex_out_A <- aldex(aldex_in_A,
                       samples_A$timepoint,
                       mc.samples = 128,
                       test = "t",
                       paired.test = TRUE,
                       effect = TRUE,
                       include.sample.summary = FALSE, 
                       denom = "all",
                       verbose = FALSE) %>%
    rownames_to_column('Genus')
  
  
  aldex_out_B <- aldex(aldex_in_B,
                       samples_B$timepoint,
                       mc.samples = 128,
                       test = "t",
                       paired.test = TRUE,
                       effect = TRUE,
                       include.sample.summary = FALSE, 
                       denom = "all",
                       verbose = FALSE) %>%
    rownames_to_column('Genus')
  
  aldex_out_C <- aldex(aldex_in_C,
                       samples_C$timepoint,
                       mc.samples = 128,
                       test = "t",
                       paired.test = TRUE,
                       effect = TRUE,
                       include.sample.summary = FALSE, 
                       denom = "all",
                       verbose = FALSE)  %>%
    rownames_to_column('Genus')
  
  list(GroupA = aldex_out_A,
       GroupB = aldex_out_B,
       GroupC = aldex_out_C) %>%
    writexl::write_xlsx('aldex2_genus_timepoint_all.xlsx')
  
  list(GroupA = aldex_out_A %>%
         filter(wi.eBH <= 0.05),
       GroupB = aldex_out_B %>%
         filter(wi.eBH <= 0.05),
       GroupC = aldex_out_C %>%
         filter(wi.eBH <= 0.05)) %>%
    writexl::write_xlsx('aldex2_genus_timepoint_sig.xlsx')
  
  sig_genera_A <- aldex_out_A %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Genus)
  sig_genera_B <- aldex_out_B %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Genus)
  sig_genera_C <- aldex_out_C %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Genus)
  
  ra_summary_A <- genus_table %>%
    data.frame() %>%
    rownames_to_column('Genus') %>%
    pivot_longer(cols = -Genus, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Genus %in% sig_genera_A) %>%
    dplyr::select(Genus, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(group == 'Chow') %>%
    group_by(Genus, timepoint) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Genus', names_from = 'timepoint', names_prefix = 'RA_', values_from = 'MeanRA')
  
  ra_summary_B <- genus_table %>%
    data.frame() %>%
    rownames_to_column('Genus') %>%
    pivot_longer(cols = -Genus, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Genus %in% sig_genera_B) %>%
    dplyr::select(Genus, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(group == 'IA') %>%
    group_by(Genus, timepoint) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Genus', names_from = 'timepoint', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  ra_summary_C <- genus_table %>%
    data.frame() %>%
    rownames_to_column('Genus') %>%
    pivot_longer(cols = -Genus, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Genus %in% sig_genera_C) %>%
    dplyr::select(Genus, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(group == 'Ad lib') %>%
    group_by(Genus, timepoint) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Genus', names_from = 'timepoint', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  list(GroupA = ra_summary_A,
       GroupB = ra_summary_B,
       GroupC = ra_summary_C) %>%
    writexl::write_xlsx('aldex2_genus_timepoint_sig_rasummary.xlsx')
}

# identify genera that are differentially abundant between groups within each timepoint
{
  samples_T1AB <- metadata %>%
    filter(group %in% c('Chow', 'IA') & timepoint == 'T1')
  samples_T1AC <- metadata %>%
    filter(group %in% c('Chow', 'Ad lib') & timepoint == 'T1')
  samples_T1BC <- metadata %>%
    filter(group %in% c('IA', 'Ad lib') & timepoint == 'T1')
  samples_T2AB <- metadata %>%
    filter(group %in% c('Chow', 'IA') & timepoint == 'T2')
  samples_T2AC <- metadata %>%
    filter(group %in% c('Chow', 'Ad lib') & timepoint == 'T2')
  samples_T2BC <- metadata %>%
    filter(group %in% c('IA', 'Ad lib') & timepoint == 'T2')
  
  aldex_in_T1AB <- genus_table_filt %>%
    dplyr::select(samples_T1AB$sample_name)
  aldex_in_T1AC <- genus_table_filt %>%
    dplyr::select(samples_T1AC$sample_name)
  aldex_in_T1BC <- genus_table_filt %>%
    dplyr::select(samples_T1BC$sample_name)
  aldex_in_T2AB <- genus_table_filt %>%
    dplyr::select(samples_T1AB$sample_name)
  aldex_in_T2AC <- genus_table_filt %>%
    dplyr::select(samples_T2AC$sample_name)
  aldex_in_T2BC <- genus_table_filt %>%
    dplyr::select(samples_T2BC$sample_name)
  
  aldex_out_T1AB <- aldex(aldex_in_T1AB, samples_T1AB$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>% rownames_to_column("Genus")
  aldex_out_T1AC <- aldex(aldex_in_T1AC, samples_T1AC$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>% rownames_to_column("Genus")
  aldex_out_T1BC <- aldex(aldex_in_T1BC, samples_T1BC$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>% rownames_to_column("Genus") 
  aldex_out_T2AB <- aldex(aldex_in_T2AB, samples_T2AB$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>% rownames_to_column("Genus") 
  aldex_out_T2AC <- aldex(aldex_in_T2AC, samples_T2AC$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>% rownames_to_column("Genus") 
  aldex_out_T2BC <- aldex(aldex_in_T2BC, samples_T2BC$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>% rownames_to_column("Genus") 
  
  list(T1AB = aldex_out_T1AB,
       T1AC = aldex_out_T1AC,
       T1BC = aldex_out_T1BC,
       T2AB = aldex_out_T2AB,
       T2AC = aldex_out_T2AC,
       T2BC = aldex_out_T2BC) %>%
    writexl::write_xlsx('aldex2_genus_group_all.xlsx')
  
  list(T1AB = aldex_out_T1AB %>%
         filter(wi.eBH <= 0.05),
       T1AC = aldex_out_T1AC %>%
         filter(wi.eBH <= 0.05),
       T1BC = aldex_out_T1BC %>%
         filter(wi.eBH <= 0.05),
       T2AB = aldex_out_T2AB %>%
         filter(wi.eBH <= 0.05),
       T2AC = aldex_out_T2AC %>%
         filter(wi.eBH <= 0.05),
       T2BC = aldex_out_T2BC %>%
         filter(wi.eBH <= 0.05)) %>%
    writexl::write_xlsx('aldex2_genus_group_sig.xlsx')
  
  sig_genus_T1AB <- aldex_out_T1AB %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Genus)
  sig_genus_T1AC <- aldex_out_T1AC %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Genus)
  sig_genus_T1BC <- aldex_out_T1BC %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Genus)
  sig_genus_T2AB <- aldex_out_T2AB %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Genus)
  sig_genus_T2AC <- aldex_out_T2AC %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Genus)
  sig_genus_T2BC <- aldex_out_T2BC %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Genus)
  
  ra_summary_T1AB <- genus_table_filt %>%
    data.frame() %>%
    rownames_to_column('Genus') %>%
    pivot_longer(cols = -Genus, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Genus %in% sig_genus_T1AB) %>%
    dplyr::select(Genus, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T1' & group %in% c('Chow', 'IA')) %>%
    group_by(Genus, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Genus', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA')
  
  ra_summary_T1AC <- genus_table_filt %>%
    data.frame() %>%
    rownames_to_column('Genus') %>%
    pivot_longer(cols = -Genus, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Genus %in% sig_genus_T1AC) %>%
    dplyr::select(Genus, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T1' & group %in% c('Chow', 'Ad lib')) %>%
    group_by(Genus, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Genus', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  ra_summary_T1BC <- genus_table_filt %>%
    data.frame() %>%
    rownames_to_column('Genus') %>%
    pivot_longer(cols = -Genus, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Genus %in% sig_genus_T1BC) %>%
    dplyr::select(Genus, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T1' & group %in% c('IA', 'Ad lib')) %>%
    group_by(Genus, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Genus', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  ra_summary_T2AB <- genus_table_filt %>%
    data.frame() %>%
    rownames_to_column('Genus') %>%
    pivot_longer(cols = -Genus, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Genus %in% sig_genus_T2AB) %>%
    dplyr::select(Genus, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T2' & group %in% c('Chow', 'IA')) %>%
    group_by(Genus, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Genus', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  ra_summary_T2AC <- genus_table_filt %>%
    data.frame() %>%
    rownames_to_column('Genus') %>%
    pivot_longer(cols = -Genus, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Genus %in% sig_genus_T2AC) %>%
    dplyr::select(Genus, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T2' & group %in% c('Chow', 'Ad lib')) %>%
    group_by(Genus, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Genus', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  ra_summary_T2BC <- genus_table_filt %>%
    data.frame() %>%
    rownames_to_column('Genus') %>%
    pivot_longer(cols = -Genus, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Genus %in% sig_genus_T2BC) %>%
    dplyr::select(Genus, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T2' & group %in% c('IA', 'Ad lib')) %>%
    group_by(Genus, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Genus', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  list(T1AB = ra_summary_T1AB,
       T1AC = ra_summary_T1AB,
       T1BC = ra_summary_T1BC,
       T2AB = ra_summary_T2AB,
       T2AC = ra_summary_T2AC,
       T2BC = ra_summary_T2BC) %>%
    writexl::write_xlsx('aldex2_genus_group_sig_rasummary.xlsx')
}

# collapse feature table to Family level
{
  family_table <- feature_table %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    left_join(taxonomy, by = c('OTU' = 'Feature ID')) %>%
    separate(Taxon, 
             into = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), 
             sep = '; ') %>%
    pivot_longer(cols = starts_with('DMG'), names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample, Family) %>%
    summarise(Count = sum(Count)) %>%
    pivot_wider(id_cols = Family, names_from = 'Sample', values_from = 'Count', values_fill = 0) %>%
    mutate(Family = if_else(is.na(Family), 'Unassigned', Family)) %>%
    column_to_rownames('Family')
}

# filter family table to only include families that are > 0.1% in at least one group-timepoint combo
{
  families_to_keep <- family_table %>%
    data.frame() %>%
    rownames_to_column('Family') %>%
    pivot_longer(cols = -Family, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    dplyr::select(Family, Sample, RA) %>%
    filter(RA >= 0.1) %>%
    filter(!Family %in% c('f__', 'Unassigned')) %>%
    pull(Family) %>%
    unique()
  
  family_table_filt <- family_table %>%
    data.frame() %>%
    rownames_to_column('Family') %>%
    filter(Family %in% families_to_keep) %>%
    column_to_rownames('Family')
}

# identify families that are differentially abundant between timepoints within each group
{
  samples_A <- metadata %>%
    filter(group == 'Chow')
  samples_B <- metadata %>%
    filter(group == 'IA')
  samples_C <- metadata %>%
    filter(group == 'Ad lib')
  
  aldex_in_A <- family_table_filt %>%
    dplyr::select(samples_A$sample_name)
  aldex_in_B <- family_table_filt %>%
    dplyr::select(samples_B$sample_name)
  aldex_in_C <- family_table_filt %>%
    dplyr::select(samples_C$sample_name)
  
  aldex_out_A <- aldex(aldex_in_A,
                       samples_A$timepoint,
                       mc.samples = 128,
                       test = "t",
                       paired.test = TRUE,
                       effect = TRUE,
                       include.sample.summary = FALSE, 
                       denom = "all",
                       verbose = FALSE) %>%
    rownames_to_column('Family')
  
  
  aldex_out_B <- aldex(aldex_in_B,
                       samples_B$timepoint,
                       mc.samples = 128,
                       test = "t",
                       paired.test = TRUE,
                       effect = TRUE,
                       include.sample.summary = FALSE, 
                       denom = "all",
                       verbose = FALSE) %>%
    rownames_to_column('Family')
  
  aldex_out_C <- aldex(aldex_in_C,
                       samples_C$timepoint,
                       mc.samples = 128,
                       test = "t",
                       paired.test = TRUE,
                       effect = TRUE,
                       include.sample.summary = FALSE, 
                       denom = "all",
                       verbose = FALSE)  %>%
    rownames_to_column('Family')
  
  list(GroupA = aldex_out_A,
       GroupB = aldex_out_B,
       GroupC = aldex_out_C) %>%
    writexl::write_xlsx('aldex2_family_timepoint_all.xlsx')
  
  list(GroupA = aldex_out_A %>%
         filter(wi.eBH <= 0.05),
       GroupB = aldex_out_B %>%
         filter(wi.eBH <= 0.05),
       GroupC = aldex_out_C %>%
         filter(wi.eBH <= 0.05)) %>%
    writexl::write_xlsx('aldex2_family_timepoint_sig.xlsx')
  
  sig_families_A <- aldex_out_A %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Family)
  sig_families_B <- aldex_out_B %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Family)
  sig_families_C <- aldex_out_C %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Family)
  
  ra_summary_A <- family_table %>%
    data.frame() %>%
    rownames_to_column('Family') %>%
    pivot_longer(cols = -Family, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Family %in% sig_families_A) %>%
    dplyr::select(Family, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(group == 'Chow') %>%
    group_by(Family, timepoint) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Family', names_from = 'timepoint', names_prefix = 'RA_', values_from = 'MeanRA')
  
  ra_summary_B <- family_table %>%
    data.frame() %>%
    rownames_to_column('Family') %>%
    pivot_longer(cols = -Family, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Family %in% sig_families_B) %>%
    dplyr::select(Family, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(group == 'IA') %>%
    group_by(Family, timepoint) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Family', names_from = 'timepoint', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  ra_summary_C <- family_table %>%
    data.frame() %>%
    rownames_to_column('Family') %>%
    pivot_longer(cols = -Family, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Family %in% sig_families_C) %>%
    dplyr::select(Family, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(group == 'Ad lib') %>%
    group_by(Family, timepoint) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Family', names_from = 'timepoint', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  list(GroupA = ra_summary_A,
       GroupB = ra_summary_B,
       GroupC = ra_summary_C) %>%
    writexl::write_xlsx('aldex2_family_timepoint_sig_rasummary.xlsx')
}

# identify families that are differentially abundant between groups within each timepoint
{
  samples_T1AB <- metadata %>%
    filter(group %in% c('Chow', 'IA') & timepoint == 'T1')
  samples_T1AC <- metadata %>%
    filter(group %in% c('Chow', 'Ad lib') & timepoint == 'T1')
  samples_T1BC <- metadata %>%
    filter(group %in% c('IA', 'Ad lib') & timepoint == 'T1')
  samples_T2AB <- metadata %>%
    filter(group %in% c('Chow', 'IA') & timepoint == 'T2')
  samples_T2AC <- metadata %>%
    filter(group %in% c('Chow', 'Ad lib') & timepoint == 'T2')
  samples_T2BC <- metadata %>%
    filter(group %in% c('IA', 'Ad lib') & timepoint == 'T2')
  
  aldex_in_T1AB <- family_table_filt %>%
    dplyr::select(samples_T1AB$sample_name)
  aldex_in_T1AC <- family_table_filt %>%
    dplyr::select(samples_T1AC$sample_name)
  aldex_in_T1BC <- family_table_filt %>%
    dplyr::select(samples_T1BC$sample_name)
  aldex_in_T2AB <- family_table_filt %>%
    dplyr::select(samples_T1AB$sample_name)
  aldex_in_T2AC <- family_table_filt %>%
    dplyr::select(samples_T2AC$sample_name)
  aldex_in_T2BC <- family_table_filt %>%
    dplyr::select(samples_T2BC$sample_name)
  
  aldex_out_T1AB <- aldex(aldex_in_T1AB, samples_T1AB$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>% rownames_to_column("Family")
  aldex_out_T1AC <- aldex(aldex_in_T1AC, samples_T1AC$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>% rownames_to_column("Family")
  aldex_out_T1BC <- aldex(aldex_in_T1BC, samples_T1BC$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>% rownames_to_column("Family") 
  aldex_out_T2AB <- aldex(aldex_in_T2AB, samples_T2AB$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>% rownames_to_column("Family") 
  aldex_out_T2AC <- aldex(aldex_in_T2AC, samples_T2AC$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>% rownames_to_column("Family") 
  aldex_out_T2BC <- aldex(aldex_in_T2BC, samples_T2BC$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>% rownames_to_column("Family") 
  
  list(T1AB = aldex_out_T1AB,
       T1AC = aldex_out_T1AC,
       T1BC = aldex_out_T1BC,
       T2AB = aldex_out_T2AB,
       T2AC = aldex_out_T2AC,
       T2BC = aldex_out_T2BC) %>%
    writexl::write_xlsx('aldex2_family_group_all.xlsx')
  
  list(T1AB = aldex_out_T1AB %>%
         filter(wi.eBH <= 0.05),
       T1AC = aldex_out_T1AC %>%
         filter(wi.eBH <= 0.05),
       T1BC = aldex_out_T1BC %>%
         filter(wi.eBH <= 0.05),
       T2AB = aldex_out_T2AB %>%
         filter(wi.eBH <= 0.05),
       T2AC = aldex_out_T2AC %>%
         filter(wi.eBH <= 0.05),
       T2BC = aldex_out_T2BC %>%
         filter(wi.eBH <= 0.05)) %>%
    writexl::write_xlsx('aldex2_family_group_sig.xlsx')
  
  sig_family_T1AB <- aldex_out_T1AB %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Family)
  sig_family_T1AC <- aldex_out_T1AC %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Family)
  sig_family_T1BC <- aldex_out_T1BC %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Family)
  sig_family_T2AB <- aldex_out_T2AB %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Family)
  sig_family_T2AC <- aldex_out_T2AC %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Family)
  sig_family_T2BC <- aldex_out_T2BC %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Family)
  
  ra_summary_T1AB <- family_table_filt %>%
    data.frame() %>%
    rownames_to_column('Family') %>%
    pivot_longer(cols = -Family, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Family %in% sig_family_T1AB) %>%
    dplyr::select(Family, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T1' & group %in% c('Chow', 'IA')) %>%
    group_by(Family, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Family', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA')
  
  ra_summary_T1AC <- family_table_filt %>%
    data.frame() %>%
    rownames_to_column('Family') %>%
    pivot_longer(cols = -Family, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Family %in% sig_family_T1AC) %>%
    dplyr::select(Family, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T1' & group %in% c('Chow', 'Ad lib')) %>%
    group_by(Family, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Family', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  ra_summary_T1BC <- family_table_filt %>%
    data.frame() %>%
    rownames_to_column('Family') %>%
    pivot_longer(cols = -Family, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Family %in% sig_family_T1BC) %>%
    dplyr::select(Family, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T1' & group %in% c('IA', 'Ad lib')) %>%
    group_by(Family, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Family', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  ra_summary_T2AB <- family_table_filt %>%
    data.frame() %>%
    rownames_to_column('Family') %>%
    pivot_longer(cols = -Family, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Family %in% sig_family_T2AB) %>%
    dplyr::select(Family, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T2' & group %in% c('Chow', 'IA')) %>%
    group_by(Family, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Family', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  ra_summary_T2AC <- family_table_filt %>%
    data.frame() %>%
    rownames_to_column('Family') %>%
    pivot_longer(cols = -Family, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Family %in% sig_family_T2AC) %>%
    dplyr::select(Family, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T2' & group %in% c('Chow', 'Ad lib')) %>%
    group_by(Family, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Family', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  ra_summary_T2BC <- family_table_filt %>%
    data.frame() %>%
    rownames_to_column('Family') %>%
    pivot_longer(cols = -Family, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Family %in% sig_family_T2BC) %>%
    dplyr::select(Family, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T2' & group %in% c('IA', 'Ad lib')) %>%
    group_by(Family, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Family', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  list(T1AB = ra_summary_T1AB,
       T1AC = ra_summary_T1AB,
       T1BC = ra_summary_T1BC,
       T2AB = ra_summary_T2AB,
       T2AC = ra_summary_T2AC,
       T2BC = ra_summary_T2BC) %>%
    writexl::write_xlsx('aldex2_family_group_sig_rasummary.xlsx')
}

# collapse feature table to Phylum level
{
  phylum_table <- feature_table %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    left_join(taxonomy, by = c('OTU' = 'Feature ID')) %>%
    separate(Taxon, 
             into = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), 
             sep = '; ') %>%
    pivot_longer(cols = starts_with('DMG'), names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample, Phylum) %>%
    summarise(Count = sum(Count)) %>%
    pivot_wider(id_cols = Phylum, names_from = 'Sample', values_from = 'Count', values_fill = 0) %>%
    mutate(Phylum = if_else(is.na(Phylum), 'Unassigned', Phylum)) %>%
    column_to_rownames('Phylum')
}

# filter phylum table to only include phyla that are > 0.1% in at least one group-timepoint combo
{
  phyla_to_keep <- phylum_table %>%
    data.frame() %>%
    rownames_to_column('Phylum') %>%
    pivot_longer(cols = -Phylum, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    dplyr::select(Phylum, Sample, RA) %>%
    filter(RA >= 0.1) %>%
    filter(!Phylum %in% c('p__', 'Unassigned')) %>%
    pull(Phylum) %>%
    unique()
  
  phylum_table_filt <- phylum_table %>%
    data.frame() %>%
    rownames_to_column('Phylum') %>%
    filter(Phylum %in% phyla_to_keep) %>%
    column_to_rownames('Phylum')
}

# identify phyla that are differentially abundant between timepoints within each group
{
  samples_A <- metadata %>%
    filter(group == 'Chow')
  samples_B <- metadata %>%
    filter(group == 'IA')
  samples_C <- metadata %>%
    filter(group == 'Ad lib')
  
  aldex_in_A <- phylum_table_filt %>%
    dplyr::select(samples_A$sample_name)
  aldex_in_B <- phylum_table_filt %>%
    dplyr::select(samples_B$sample_name)
  aldex_in_C <- phylum_table_filt %>%
    dplyr::select(samples_C$sample_name)
  
  aldex_out_A <- aldex(aldex_in_A,
                       samples_A$timepoint,
                       mc.samples = 128,
                       test = "t",
                       paired.test = TRUE,
                       effect = TRUE,
                       include.sample.summary = FALSE, 
                       denom = "all",
                       verbose = FALSE) %>%
    rownames_to_column('Phylum')
  
  
  aldex_out_B <- aldex(aldex_in_B,
                       samples_B$timepoint,
                       mc.samples = 128,
                       test = "t",
                       paired.test = TRUE,
                       effect = TRUE,
                       include.sample.summary = FALSE, 
                       denom = "all",
                       verbose = FALSE) %>%
    rownames_to_column('Phylum')
  
  aldex_out_C <- aldex(aldex_in_C,
                       samples_C$timepoint,
                       mc.samples = 128,
                       test = "t",
                       paired.test = TRUE,
                       effect = TRUE,
                       include.sample.summary = FALSE, 
                       denom = "all",
                       verbose = FALSE)  %>%
    rownames_to_column('Phylum')
  
  list(GroupA = aldex_out_A,
       GroupB = aldex_out_B,
       GroupC = aldex_out_C) %>%
    writexl::write_xlsx('aldex2_phylum_timepoint_all.xlsx')
  
  list(GroupA = aldex_out_A %>%
         filter(wi.eBH <= 0.05),
       GroupB = aldex_out_B %>%
         filter(wi.eBH <= 0.05),
       GroupC = aldex_out_C %>%
         filter(wi.eBH <= 0.05)) %>%
    writexl::write_xlsx('aldex2_phylum_timepoint_sig.xlsx')
  
  sig_phyla_A <- aldex_out_A %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Phylum)
  sig_phyla_B <- aldex_out_B %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Phylum)
  sig_phyla_C <- aldex_out_C %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Phylum)
  
  ra_summary_A <- phylum_table %>%
    data.frame() %>%
    rownames_to_column('Phylum') %>%
    pivot_longer(cols = -Phylum, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Phylum %in% sig_phyla_A) %>%
    dplyr::select(Phylum, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(group == 'Chow') %>%
    group_by(Phylum, timepoint) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Phylum', names_from = 'timepoint', names_prefix = 'RA_', values_from = 'MeanRA')
  
  ra_summary_B <- phylum_table %>%
    data.frame() %>%
    rownames_to_column('Phylum') %>%
    pivot_longer(cols = -Phylum, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Phylum %in% sig_phyla_B) %>%
    dplyr::select(Phylum, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(group == 'IA') %>%
    group_by(Phylum, timepoint) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Phylum', names_from = 'timepoint', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  ra_summary_C <- phylum_table %>%
    data.frame() %>%
    rownames_to_column('Phylum') %>%
    pivot_longer(cols = -Phylum, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Phylum %in% sig_phyla_C) %>%
    dplyr::select(Phylum, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(group == 'Ad lib') %>%
    group_by(Phylum, timepoint) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Phylum', names_from = 'timepoint', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  list(GroupA = ra_summary_A,
       GroupB = ra_summary_B,
       GroupC = ra_summary_C) %>%
    writexl::write_xlsx('aldex2_phylum_timepoint_sig_rasummary.xlsx')
}

# identify phyla that are differentially abundant between groups within each timepoint
{
  samples_T1AB <- metadata %>%
    filter(group %in% c('Chow', 'IA') & timepoint == 'T1')
  samples_T1AC <- metadata %>%
    filter(group %in% c('Chow', 'Ad lib') & timepoint == 'T1')
  samples_T1BC <- metadata %>%
    filter(group %in% c('IA', 'Ad lib') & timepoint == 'T1')
  samples_T2AB <- metadata %>%
    filter(group %in% c('Chow', 'IA') & timepoint == 'T2')
  samples_T2AC <- metadata %>%
    filter(group %in% c('Chow', 'Ad lib') & timepoint == 'T2')
  samples_T2BC <- metadata %>%
    filter(group %in% c('IA', 'Ad lib') & timepoint == 'T2')
  
  aldex_in_T1AB <- phylum_table_filt %>%
    dplyr::select(samples_T1AB$sample_name)
  aldex_in_T1AC <- phylum_table_filt %>%
    dplyr::select(samples_T1AC$sample_name)
  aldex_in_T1BC <- phylum_table_filt %>%
    dplyr::select(samples_T1BC$sample_name)
  aldex_in_T2AB <- phylum_table_filt %>%
    dplyr::select(samples_T1AB$sample_name)
  aldex_in_T2AC <- phylum_table_filt %>%
    dplyr::select(samples_T2AC$sample_name)
  aldex_in_T2BC <- phylum_table_filt %>%
    dplyr::select(samples_T2BC$sample_name)
  
  aldex_out_T1AB <- aldex(aldex_in_T1AB, samples_T1AB$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>% rownames_to_column("Phylum")
  aldex_out_T1AC <- aldex(aldex_in_T1AC, samples_T1AC$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>% rownames_to_column("Phylum")
  aldex_out_T1BC <- aldex(aldex_in_T1BC, samples_T1BC$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>% rownames_to_column("Phylum") 
  aldex_out_T2AB <- aldex(aldex_in_T2AB, samples_T2AB$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>% rownames_to_column("Phylum") 
  aldex_out_T2AC <- aldex(aldex_in_T2AC, samples_T2AC$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>% rownames_to_column("Phylum") 
  aldex_out_T2BC <- aldex(aldex_in_T2BC, samples_T2BC$group, mc.samples = 128, test = "t", effect = TRUE, 
                          include.sample.summary = FALSE, denom = "all", verbose = FALSE) %>% rownames_to_column("Phylum") 
  
  list(T1AB = aldex_out_T1AB,
       T1AC = aldex_out_T1AC,
       T1BC = aldex_out_T1BC,
       T2AB = aldex_out_T2AB,
       T2AC = aldex_out_T2AC,
       T2BC = aldex_out_T2BC) %>%
    writexl::write_xlsx('aldex2_phylum_group_all.xlsx')
  
  list(T1AB = aldex_out_T1AB %>%
         filter(wi.eBH <= 0.05),
       T1AC = aldex_out_T1AC %>%
         filter(wi.eBH <= 0.05),
       T1BC = aldex_out_T1BC %>%
         filter(wi.eBH <= 0.05),
       T2AB = aldex_out_T2AB %>%
         filter(wi.eBH <= 0.05),
       T2AC = aldex_out_T2AC %>%
         filter(wi.eBH <= 0.05),
       T2BC = aldex_out_T2BC %>%
         filter(wi.eBH <= 0.05)) %>%
    writexl::write_xlsx('aldex2_phylum_group_sig.xlsx')
  
  sig_phylum_T1AB <- aldex_out_T1AB %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Phylum)
  sig_phylum_T1AC <- aldex_out_T1AC %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Phylum)
  sig_phylum_T1BC <- aldex_out_T1BC %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Phylum)
  sig_phylum_T2AB <- aldex_out_T2AB %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Phylum)
  sig_phylum_T2AC <- aldex_out_T2AC %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Phylum)
  sig_phylum_T2BC <- aldex_out_T2BC %>%
    filter(wi.eBH <= 0.05) %>%
    pull(Phylum)
  
  ra_summary_T1AB <- phylum_table_filt %>%
    data.frame() %>%
    rownames_to_column('Phylum') %>%
    pivot_longer(cols = -Phylum, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Phylum %in% sig_phylum_T1AB) %>%
    dplyr::select(Phylum, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T1' & group %in% c('Chow', 'IA')) %>%
    group_by(Phylum, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Phylum', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA')
  
  ra_summary_T1AC <- phylum_table_filt %>%
    data.frame() %>%
    rownames_to_column('Phylum') %>%
    pivot_longer(cols = -Phylum, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Phylum %in% sig_phylum_T1AC) %>%
    dplyr::select(Phylum, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T1' & group %in% c('Chow', 'Ad lib')) %>%
    group_by(Phylum, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Phylum', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  ra_summary_T1BC <- phylum_table_filt %>%
    data.frame() %>%
    rownames_to_column('Phylum') %>%
    pivot_longer(cols = -Phylum, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Phylum %in% sig_phylum_T1BC) %>%
    dplyr::select(Phylum, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T1' & group %in% c('IA', 'Ad lib')) %>%
    group_by(Phylum, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Phylum', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  ra_summary_T2AB <- phylum_table_filt %>%
    data.frame() %>%
    rownames_to_column('Phylum') %>%
    pivot_longer(cols = -Phylum, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Phylum %in% sig_phylum_T2AB) %>%
    dplyr::select(Phylum, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T2' & group %in% c('Chow', 'IA')) %>%
    group_by(Phylum, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Phylum', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  ra_summary_T2AC <- phylum_table_filt %>%
    data.frame() %>%
    rownames_to_column('Phylum') %>%
    pivot_longer(cols = -Phylum, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Phylum %in% sig_phylum_T2AC) %>%
    dplyr::select(Phylum, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T2' & group %in% c('Chow', 'Ad lib')) %>%
    group_by(Phylum, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Phylum', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  ra_summary_T2BC <- phylum_table_filt %>%
    data.frame() %>%
    rownames_to_column('Phylum') %>%
    pivot_longer(cols = -Phylum, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Phylum %in% sig_phylum_T2BC) %>%
    dplyr::select(Phylum, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    filter(timepoint == 'T2' & group %in% c('IA', 'Ad lib')) %>%
    group_by(Phylum, group) %>%
    summarise(MeanRA = mean(RA)) %>%
    pivot_wider(id_cols = 'Phylum', names_from = 'group', names_prefix = 'RA_', values_from = 'MeanRA') 
  
  list(T1AB = ra_summary_T1AB,
       T1AC = ra_summary_T1AB,
       T1BC = ra_summary_T1BC,
       T2AB = ra_summary_T2AB,
       T2AC = ra_summary_T2AC,
       T2BC = ra_summary_T2BC) %>%
    writexl::write_xlsx('aldex2_phylum_group_sig_rasummary.xlsx')
}

# plotting relative abundance of differentially abundant phyla by group and timepoint
{
  ra_diffabund <- phylum_table_filt %>%
    rownames_to_column('Phylum') %>%
    pivot_longer(cols = -Phylum, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    filter(Phylum %in% c('p__Actinobacteria',
                         'p__Bacteroidetes',
                         'p__Firmicutes')) %>%
    dplyr::select(Phylum, Sample, RA) %>%
    left_join(metadata,
              by = c('Sample' = 'sample_name')) %>%
    mutate(Phylum = str_replace(Phylum, '^p__', '')) %>%
    write_delim('ra_differentialphyla.tsv', delim = '\t')
}

# one way anova comparing alpha diversity at T1 and T2
{
  alpha_df <- read_tsv('alpha_diversity.tsv')
  
  alpha_df %>%
    filter(timepoint == 'T1') %>%
    anova_test(dv = shannon,
               between = group)
  
  alpha_df %>%
    filter(timepoint == 'T2') %>%
    anova_test(dv = shannon,
               between = group)
}

# testing for differences in alpha diversity using linear mixed effect models
{
  alpha_df <- read_tsv('alpha_diversity.tsv')
  
  alpha_df <- alpha_df %>%
    mutate(group = factor(group, levels = c('Chow', 'Ad lib', 'IA'))) %>%
    mutate(mouse = factor(mouse)) %>%
    mutate(timepoint = factor(timepoint))
  
  alpha_df %>%
    ggplot(aes(x = timepoint, y = shannon, group = group)) +
    geom_point(aes(fill = group), position = position_jitter(width = 0.1),
               size = 3, colour = 'black', shape = 21) +
    theme_bw() +
    geom_smooth(aes(colour = group), method = 'lm', se = F)
  
  model1 <- lme(shannon ~ group * timepoint,
               random = ~ 1 | mouse,
               data = alpha_df)
  
  model2 <- lme(shannon ~ group * timepoint,
               random = ~ 1 + timepoint | mouse,
               data = alpha_df)
  
  anova(model1)
  anova(model2)
  
  anova(model1, model2)
}

# import mouse behaviour data etc.
{
  behaviour_T1 <- readxl::read_xlsx('metadata_behaviour.xlsx', sheet = 'T1')
  behaviour_T2 <- readxl::read_xlsx('metadata_behaviour.xlsx', sheet = 'T2')
}

# look for relationships between microbiome features and behaviour data using correlations
{
  library('psych')
  
  in_behaviour_T2 <- behaviour_T2 %>% 
    column_to_rownames('DAMG_ID') %>% 
    select(1, 5:85)
  
  feature_table <- phyloseq::import_biom('feature-table_json.biom')
  
  otus_to_keep <- feature_table %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    dplyr::select(OTU, Sample, RA) %>%
    filter(RA >= 0.1) %>%
    pull(OTU) %>%
    unique()
  
  feature_table_filt_ra <- feature_table %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    dplyr::select(OTU, Sample, RA) %>%
    filter(OTU %in% otus_to_keep) %>%
    pivot_wider(id_cols = Sample, names_from = OTU, values_from = RA, values_fill = 0) %>%
    filter(Sample %in% rownames(in_behaviour_T2)) %>%
    column_to_rownames('Sample')
  
  corrout_all_otu <- psych::corr.test(x = feature_table_filt_ra, 
                                      y = in_behaviour_T2 %>% select(-Group), 
                                      method = 'spearman', 
                                      adjust = 'holm')
  
  in_behaviour <- in_behaviour_T2 %>% filter(Group == 'Chow')
  in_feat <- feature_table_filt_ra[rownames(in_behaviour),]
  
  corrout_chow_otu <- psych::corr.test(x = in_feat, 
                                       y = in_behaviour %>% select(-Group), 
                                       method = 'spearman', 
                                       adjust = 'holm')
  
  in_behaviour <- in_behaviour_T2 %>% filter(Group == 'IA')
  in_feat <- feature_table_filt_ra[rownames(in_behaviour),]
  
  corrout_IA_otu <- psych::corr.test(x = in_feat, 
                                     y = in_behaviour %>% select(-Group), 
                                     method = 'spearman', 
                                     adjust = 'holm')
  
  in_behaviour <- in_behaviour_T2 %>% filter(Group == 'Ad lib')
  in_feat <- feature_table_filt_ra[rownames(in_behaviour),]
  
  corrout_Adlib_otu <- psych::corr.test(x = in_feat, 
                                        y = in_behaviour %>% select(-Group), 
                                        method = 'spearman', 
                                        adjust = 'holm')
  
  p_all <- corrout_all_otu$p.adj %>%
    data.frame() %>% 
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Variable', values_to = 'P_holm') %>%
    filter(P_holm <= 0.05) %>%
    left_join(taxonomy, by = c('OTU' = 'Feature ID'))
  
  p_chow <- corrout_chow_otu$p.adj %>%
    data.frame() %>% 
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Variable', values_to = 'P_holm') %>%
    filter(P_holm <= 0.05) %>%
    left_join(taxonomy, by = c('OTU' = 'Feature ID'))
  
  p_IA <- corrout_IA_otu$p.adj %>%
    data.frame() %>% 
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Variable', values_to = 'P_holm') %>%
    filter(P_holm <= 0.05)  %>%
    left_join(taxonomy, by = c('OTU' = 'Feature ID'))
  
  p_adlib <- corrout_Adlib_otu$p.adj %>%
    data.frame() %>% 
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Variable', values_to = 'P_holm') %>%
    filter(P_holm <= 0.05)  %>%
    left_join(taxonomy, by = c('OTU' = 'Feature ID'))
  
  list(All = p_all,
       Chow = p_chow,
       IA = p_IA,
       Adlib = p_adlib) %>%
    writexl::write_xlsx('corrsig_otu.xlsx')
  
  r_all <- corrout_all_otu$r %>%
    data.frame() %>% 
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Variable', values_to = 'R') %>%
    left_join(taxonomy, by = c('OTU' = 'Feature ID'))
  
  r_chow <- corrout_chow_otu$r %>%
    data.frame() %>% 
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Variable', values_to = 'R') %>%
    left_join(taxonomy, by = c('OTU' = 'Feature ID'))
  
  r_IA <- corrout_IA_otu$r %>%
    data.frame() %>% 
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Variable', values_to = 'R') %>%
    left_join(taxonomy, by = c('OTU' = 'Feature ID'))
  
  r_adlib <- corrout_Adlib_otu$r %>%
    data.frame() %>% 
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Variable', values_to = 'R') %>%
    left_join(taxonomy, by = c('OTU' = 'Feature ID'))
  
  list(All = r_all,
       Chow = r_chow,
       IA = r_IA,
       Adlib = r_adlib) %>%
    writexl::write_xlsx('corrR_otu.xlsx')
  
  list(All = p_all %>% left_join(r_all, by = c('OTU', 'Variable')),
       Chow = p_chow %>% left_join(r_chow, by = c('OTU', 'Variable')),
       IA = p_IA %>% left_join(r_IA, by = c('OTU', 'Variable')),
       Adlib = p_adlib %>% left_join(r_adlib, by = c('OTU', 'Variable'))) %>%
    writexl::write_xlsx('corrsigR_otu.xlsx')
  
  # performing correlations at genus level
  
  genus_to_keep <- feature_table %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Sample', values_to = 'Count') %>%
    left_join(taxonomy, by = c('OTU' = 'Feature ID')) %>%
    mutate(Taxon = str_remove(Taxon, '; s__.*')) %>%
    group_by(Sample, Taxon) %>%
    summarise(Count = sum(Count)) %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    dplyr::select(Taxon, Sample, RA) %>%
    filter(RA >= 0.1) %>%
    pull(Taxon) %>%
    unique()
  
  genus_table_filt_ra <- feature_table %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Sample', values_to = 'Count') %>%
    left_join(taxonomy, by = c('OTU' = 'Feature ID')) %>%
    mutate(Taxon = str_remove(Taxon, '; s__.*')) %>%
    group_by(Sample, Taxon) %>%
    summarise(Count = sum(Count)) %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    dplyr::select(Taxon, Sample, RA) %>%
    filter(Taxon %in% genus_to_keep) %>%
    pivot_wider(id_cols = Sample, names_from = Taxon, values_from = RA, values_fill = 0) %>%
    filter(Sample %in% rownames(in_behaviour_T2)) %>%
    column_to_rownames('Sample')
  
  corrout_all_genus <- psych::corr.test(x = genus_table_filt_ra,
                                        y = in_behaviour_T2 %>% select(-Group), 
                                        method = 'spearman', 
                                        adjust = 'holm')
  
  in_behaviour <- in_behaviour_T2 %>% filter(Group == 'Chow')
  in_feat <- genus_table_filt_ra[rownames(in_behaviour),]
  
  corrout_chow_genus <- psych::corr.test(x = in_feat, 
                                       y = in_behaviour %>% select(-Group), 
                                       method = 'spearman', 
                                       adjust = 'holm')
  
  in_behaviour <- in_behaviour_T2 %>% filter(Group == 'IA')
  in_feat <- genus_table_filt_ra[rownames(in_behaviour),]
  
  corrout_IA_genus <- psych::corr.test(x = in_feat,
                                       y = in_behaviour %>% select(-Group), 
                                       method = 'spearman', 
                                       adjust = 'holm')
  
  in_behaviour <- in_behaviour_T2 %>% filter(Group == 'Ad lib')
  in_feat <- genus_table_filt_ra[rownames(in_behaviour),]
  
  corrout_Adlib_genus <- psych::corr.test(x = in_feat, 
                                        y = in_behaviour %>% select(-Group), 
                                        method = 'spearman', 
                                        adjust = 'holm')
  
  p_all <- corrout_all_genus$p.adj %>%
    data.frame() %>% 
    rownames_to_column('Genus') %>%
    pivot_longer(cols = -Genus, names_to = 'Variable', values_to = 'P_holm') %>%
    filter(P_holm <= 0.05) 
  
  p_chow <- corrout_chow_genus$p.adj %>%
    data.frame() %>% 
    rownames_to_column('Genus') %>%
    pivot_longer(cols = -Genus, names_to = 'Variable', values_to = 'P_holm') %>%
    filter(P_holm <= 0.05)
  
  p_IA <- corrout_IA_genus$p.adj %>%
    data.frame() %>% 
    rownames_to_column('Genus') %>%
    pivot_longer(cols = -Genus, names_to = 'Variable', values_to = 'P_holm') %>%
    filter(P_holm <= 0.05)  
  
  p_adlib <- corrout_Adlib_genus$p.adj %>%
    data.frame() %>% 
    rownames_to_column('Genus') %>%
    pivot_longer(cols = -Genus, names_to = 'Variable', values_to = 'P_holm') %>%
    filter(P_holm <= 0.05) 
  
  list(All = p_all,
       Chow = p_chow,
       IA = p_IA,
       Adlib = p_adlib) %>%
    writexl::write_xlsx('corrsig_genus.xlsx')
  
  r_all <- corrout_all_genus$r %>%
    data.frame() %>% 
    rownames_to_column('Genus') %>%
    pivot_longer(cols = -Genus, names_to = 'Variable', values_to = 'R') 
  
  r_chow <- corrout_chow_genus$r %>%
    data.frame() %>% 
    rownames_to_column('Genus') %>%
    pivot_longer(cols = -Genus, names_to = 'Variable', values_to = 'R') 
  
  r_IA <- corrout_IA_genus$r %>%
    data.frame() %>% 
    rownames_to_column('Genus') %>%
    pivot_longer(cols = -Genus, names_to = 'Variable', values_to = 'R') 
  
  r_adlib <- corrout_Adlib_genus$r %>%
    data.frame() %>% 
    rownames_to_column('Genus') %>%
    pivot_longer(cols = -Genus, names_to = 'Variable', values_to = 'R')
  
  list(All = r_all,
       Chow = r_chow,
       IA = r_IA,
       Adlib = r_adlib) %>%
    writexl::write_xlsx('corrR_genus.xlsx')
  
  list(All = p_all %>% left_join(r_all, by = c('Genus', 'Variable')),
       Chow = p_chow %>% left_join(r_chow, by = c('Genus', 'Variable')),
       IA = p_IA %>% left_join(r_IA, by = c('Genus', 'Variable')),
       Adlib = p_adlib %>% left_join(r_adlib, by = c('Genus', 'Variable'))) %>%
    writexl::write_xlsx('corrsigR_genus.xlsx')
  
  # performing correlations at family level
  
  family_to_keep <- feature_table %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Sample', values_to = 'Count') %>%
    left_join(taxonomy, by = c('OTU' = 'Feature ID')) %>%
    mutate(Taxon = str_remove(Taxon, '; g__.*')) %>%
    group_by(Sample, Taxon) %>%
    summarise(Count = sum(Count)) %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    dplyr::select(Taxon, Sample, RA) %>%
    filter(RA >= 0.1) %>%
    pull(Taxon) %>%
    unique()
  
  family_table_filt_ra <- feature_table %>%
    data.frame() %>%
    rownames_to_column('OTU') %>%
    pivot_longer(cols = -OTU, names_to = 'Sample', values_to = 'Count') %>%
    left_join(taxonomy, by = c('OTU' = 'Feature ID')) %>%
    mutate(Taxon = str_remove(Taxon, '; g__.*')) %>%
    group_by(Sample, Taxon) %>%
    summarise(Count = sum(Count)) %>%
    group_by(Sample) %>%
    mutate(CountSum = sum(Count)) %>%
    mutate(RA = 100*(Count/CountSum)) %>%
    dplyr::select(Taxon, Sample, RA) %>%
    filter(Taxon %in% genus_to_keep) %>%
    pivot_wider(id_cols = Sample, names_from = Taxon, values_from = RA, values_fill = 0) %>%
    filter(Sample %in% rownames(in_behaviour_T2)) %>%
    column_to_rownames('Sample')
  
  corrout_all_family <- psych::corr.test(x = family_table_filt_ra,
                                        y = in_behaviour_T2 %>% select(-Group), 
                                        method = 'spearman', 
                                        adjust = 'holm')
  
  in_behaviour <- in_behaviour_T2 %>% filter(Group == 'Chow')
  in_feat <- family_table_filt_ra[rownames(in_behaviour),]
  
  corrout_chow_family <- psych::corr.test(x = in_feat, 
                                         y = in_behaviour %>% select(-Group), 
                                         method = 'spearman', 
                                         adjust = 'holm')
  
  in_behaviour <- in_behaviour_T2 %>% filter(Group == 'IA')
  in_feat <- family_table_filt_ra[rownames(in_behaviour),]
  
  corrout_IA_family <- psych::corr.test(x = in_feat, 
                                       y = in_behaviour %>% select(-Group), 
                                       method = 'spearman', 
                                       adjust = 'holm')
  
  in_behaviour <- in_behaviour_T2 %>% filter(Group == 'Ad lib')
  in_feat <- family_table_filt_ra[rownames(in_behaviour),]
  
  corrout_Adlib_family <- psych::corr.test(x = in_feat, 
                                          y = in_behaviour %>% select(-Group), 
                                          method = 'spearman', 
                                          adjust = 'holm')
  
  p_all <- corrout_all_family$p.adj %>%
    data.frame() %>% 
    rownames_to_column('Family') %>%
    pivot_longer(cols = -Family, names_to = 'Variable', values_to = 'P_holm') %>%
    filter(P_holm <= 0.05) 
  
  p_chow <- corrout_chow_family$p.adj %>%
    data.frame() %>% 
    rownames_to_column('Family') %>%
    pivot_longer(cols = -Family, names_to = 'Variable', values_to = 'P_holm') %>%
    filter(P_holm <= 0.05)
  
  p_IA <- corrout_IA_family$p.adj %>%
    data.frame() %>% 
    rownames_to_column('Family') %>%
    pivot_longer(cols = -Family, names_to = 'Variable', values_to = 'P_holm') %>%
    filter(P_holm <= 0.05)  
  
  p_adlib <- corrout_Adlib_family$p.adj %>%
    data.frame() %>% 
    rownames_to_column('Family') %>%
    pivot_longer(cols = -Family, names_to = 'Variable', values_to = 'P_holm') %>%
    filter(P_holm <= 0.05) 
  
  list(All = p_all,
       Chow = p_chow,
       IA = p_IA,
       Adlib = p_adlib) %>%
    writexl::write_xlsx('corrsig_family.xlsx')
  
  r_all <- corrout_all_family$r %>%
    data.frame() %>% 
    rownames_to_column('Family') %>%
    pivot_longer(cols = -Family, names_to = 'Variable', values_to = 'R') 
  
  r_chow <- corrout_chow_family$r %>%
    data.frame() %>% 
    rownames_to_column('Family') %>%
    pivot_longer(cols = -Family, names_to = 'Variable', values_to = 'R') 
  
  r_IA <- corrout_IA_family$r %>%
    data.frame() %>% 
    rownames_to_column('Family') %>%
    pivot_longer(cols = -Family, names_to = 'Variable', values_to = 'R') 
  
  r_adlib <- corrout_Adlib_family$r %>%
    data.frame() %>% 
    rownames_to_column('Family') %>%
    pivot_longer(cols = -Family, names_to = 'Variable', values_to = 'R')
  
  list(All = r_all,
       Chow = r_chow,
       IA = r_IA,
       Adlib = r_adlib) %>%
    writexl::write_xlsx('corrR_family.xlsx')
  
  list(All = p_all %>% left_join(r_all, by = c('Family', 'Variable')),
       Chow = p_chow %>% left_join(r_chow, by = c('Family', 'Variable')),
       IA = p_IA %>% left_join(r_IA, by = c('Family', 'Variable')),
       Adlib = p_adlib %>% left_join(r_adlib, by = c('Family', 'Variable'))) %>%
    writexl::write_xlsx('corrsigR_family.xlsx')
  
  # combining results
  rbind(readxl::read_xlsx('corrsigR_otu.xlsx', sheet = 'Chow') %>%
          select(-c('OTU', 'Taxon.y')) %>%
          relocate(Taxon.x, .before = Variable) %>%
          rename(Taxon = Taxon.x) %>%
          mutate(Rank = 'OTU') %>% 
          mutate(Group = 'Chow'),
        readxl::read_xlsx('corrsigR_genus.xlsx', sheet = 'Chow') %>%
          mutate(Rank = 'Genus') %>% 
          mutate(Group = 'Chow') %>% 
          rename(Taxon = Genus),
        readxl::read_xlsx('corrsigR_genus.xlsx', sheet = 'IA') %>%
          mutate(Rank = 'Genus') %>% 
          mutate(Group = 'IA') %>% 
          rename(Taxon = Genus),
        readxl::read_xlsx('corrsigR_family.xlsx', sheet = 'IA') %>%
          mutate(Rank = 'Family') %>% 
          mutate(Group = 'IA') %>% 
          rename(Taxon = Family)) %>%
    write_tsv('corrsigR_combined.tsv')
  
}

# plot significant correlations
{
  corrsig_otu_chow <- readxl::read_xlsx('corrsigR_otu.xlsx', sheet = 'Chow') %>%
    rename(Taxon = OTU) %>%
    mutate(Group = 'Chow', Rank = 'OTU') %>%
    select(Taxon, Variable, P_holm, R, Group, Rank)
  corrsig_genus_chow <- readxl::read_xlsx('corrsigR_genus.xlsx', sheet = 'Chow') %>%
    mutate(Group = 'Chow', Rank = 'Genus') %>% rename(Taxon = 1)
  corrsig_genus_ia <- readxl::read_xlsx('corrsigR_genus.xlsx', sheet = 'IA') %>%
    mutate(Group = 'IA', Rank = 'Genus') %>% rename(Taxon = 1)
  corrsig_family_ia <- readxl::read_xlsx('corrsigR_family.xlsx', sheet = 'IA') %>%
    mutate(Group = 'IA', Rank = 'Family') %>% rename(Taxon = 1)
  
  corrsig <- bind_rows(corrsig_otu_chow,
                       corrsig_genus_chow,
                       corrsig_genus_ia,
                       corrsig_family_ia)
  
  
  behaviour_T2 <- readxl::read_xlsx('metadata_behaviour.xlsx', sheet = 'T2') %>% 
    column_to_rownames('DAMG_ID') %>% 
    select(unique(corrsig$Variable))
  
}

# look for relationship between microbiome profile and behaviour using envfit
{
  behaviour_T2 <- readxl::read_xlsx('metadata_behaviour.xlsx', sheet = 'T2') %>% 
    column_to_rownames('DAMG_ID') %>% 
    select(1, 5:85)
  
  bc_dist<- read_qza('CoreMetricsPhylogenetic/bray_curtis_distance_matrix.qza')$data
  bc_dist_T2 <- usedist::dist_subset(bc_dist, rownames(behaviour_T2))
  
  bc_ord_T2 <- cmdscale(bc_dist_T2, eig = T)
  
  envfit_bcT2 <- envfit(bc_ord_T2, behaviour_T2, na.rm = T)
  
  envfit_out <- bind_cols(names(envfit_bcT2$vectors$r),
                          envfit_bcT2$vectors$r,
                          envfit_bcT2$vectors$pvals) %>%
    rename(Variable = 1,
           Rsquared = 2,
           Pvalue = 3) %>%
    arrange(Pvalue)
  
  behaviour_T2_chow <- behaviour_T2 %>% 
    filter(Group == 'Chow')
  behaviour_T2_adlib <- behaviour_T2 %>% 
    filter(Group == 'Ad lib')
  behaviour_T2_ia <- behaviour_T2 %>% 
    filter(Group == 'IA')
  
  bc_dist_T2_chow <- usedist::dist_subset(bc_dist, rownames(behaviour_T2_chow))
  bc_dist_T2_adlib <- usedist::dist_subset(bc_dist, rownames(behaviour_T2_adlib))
  bc_dist_T2_ia <- usedist::dist_subset(bc_dist, rownames(behaviour_T2_ia))
  
  bc_ord_T2_chow <- cmdscale(bc_dist_T2_chow, eig = T)
  bc_ord_T2_adlib <- cmdscale(bc_dist_T2_adlib, eig = T)
  bc_ord_T2_ia <- cmdscale(bc_dist_T2_ia, eig = T)
  
  envfit_bcT2_chow <- envfit(bc_ord_T2_chow, behaviour_T2_chow, na.rm = T)
  envfit_bcT2_adlib <- envfit(bc_ord_T2_adlib, behaviour_T2_adlib, na.rm = T)
  envfit_bcT2_ia <- envfit(bc_ord_T2_ia, behaviour_T2_ia, na.rm = T)
  
  envfit_out_chow <- bind_cols(names(envfit_bcT2_chow$vectors$r),
                          envfit_bcT2_chow$vectors$r,
                          envfit_bcT2_chow$vectors$pvals) %>%
    rename(Variable = 1,
           Rsquared = 2,
           Pvalue = 3) %>%
    arrange(Pvalue)
  
  envfit_out_adlib <- bind_cols(names(envfit_bcT2_adlib$vectors$r),
                          envfit_bcT2_adlib$vectors$r,
                          envfit_bcT2_adlib$vectors$pvals) %>%
    rename(Variable = 1,
           Rsquared = 2,
           Pvalue = 3) %>%
    arrange(Pvalue)
  
  envfit_out_ia <- bind_cols(names(envfit_bcT2_ia$vectors$r),
                          envfit_bcT2_ia$vectors$r,
                          envfit_bcT2_ia$vectors$pvals) %>%
    rename(Variable = 1,
           Rsquared = 2,
           Pvalue = 3) %>%
    arrange(Pvalue)
  
  list(All = envfit_out,
       Chow = envfit_out_chow,
       Ad_Lib = envfit_out_adlib,
       IA = envfit_out_ia) %>%
    writexl::write_xlsx('envfit_T2.xlsx')
  
  arrows <- envfit_bcT2$vectors$arrows %>%
    data.frame() %>%
    rownames_to_column('Variable') %>%
    filter(Variable %in% (envfit_out %>% 
             filter(Pvalue <= 0.05) %>% 
             pull(Variable))) %>%
    mutate(VariableAlias = LETTERS[1:6])
  
  arrows %>%
    write_tsv('envfitarrows.tsv')
  
  eigs <- bc_ord_T2$eig / sum(bc_ord_T2$eig)
  
  p <- bc_ord_T2$points %>%
    data.frame() %>%
    rownames_to_column('sample_name') %>%
    left_join(metadata) %>%
    ggplot(aes(x = X1, y = X2)) +
    geom_point(aes(fill = Group_Time), shape = 21, colour = 'black', size = 3) +
    stat_ellipse(aes(colour = Group_Time), linewidth = 1, show.legend = F) +
    geom_segment(data = arrows, aes(x = 0, y = 0, xend = Dim1/2, yend = Dim2/2, group = Variable), 
                 arrow = arrow(length = unit(0.2, "cm"), type = 'closed')) +
    ggrepel::geom_label_repel(data = arrows, aes(x = Dim1/2, y = Dim2/2, label = VariableAlias), fill = NA, size = 2) +
    theme_bw() +
    labs(x = paste0('PC1 [', round(100*(eigs[1]), digits = 2), '%]'),
         y = paste0('PC2 [', round(100*(eigs[2]), digits = 2), '%]')) +
    scale_fill_manual(values = groupColours) +
    scale_colour_manual(values = groupColours) 
  
  pdf('bray_T2_biplot.pdf', width = 6, height = 4)
  print(p)
  dev.off()
  
  saveRDS(p, 'bray_T2_biplot.RDS')
  
}

# export counts and relative abundances of feature table at each taxonomic rank
{
  feature_table <- phyloseq::import_biom('feature-table_json.biom')
  taxonomy <- read_delim('taxonomy.tsv', delim = '\t') %>%
    dplyr::select(-Confidence) %>%
    rename(OTU = 1) 
  
  featuretable <- feature_table %>%
    data.frame() %>%
    rownames_to_column('OTU')
  
  featuretable_long <- featuretable %>%
    pivot_longer(cols = -OTU, names_to = 'Sample', values_to = 'Count') %>%
    group_by(Sample) %>%
    mutate(RA = 100*(Count/sum(Count))) %>%
    ungroup() %>%
    left_join(taxonomy, by = 'OTU')
  
  list(Phylum = featuretable_long %>%
         mutate(Phylum = str_remove(Taxon, '; c__.*')) %>%
         group_by(Sample, Phylum) %>%
         summarise(Count = sum(Count)) %>%
         pivot_wider(id_cols = Phylum, names_from = Sample, values_from = Count, values_fill = 0),
       Class = featuretable_long %>%
         mutate(Class = str_remove(Taxon, '; o__.*')) %>%
         group_by(Sample, Class) %>%
         summarise(Count = sum(Count)) %>%
         pivot_wider(id_cols = Class, names_from = Sample, values_from = Count, values_fill = 0),
       Order = featuretable_long %>%
         mutate(Order = str_remove(Taxon, '; f__.*')) %>%
         group_by(Sample, Order) %>%
         summarise(Count = sum(Count)) %>%
         pivot_wider(id_cols = Order, names_from = Sample, values_from = Count, values_fill = 0),
       Family = featuretable_long %>%
         mutate(Family = str_remove(Taxon, '; g__.*')) %>%
         group_by(Sample, Family) %>%
         summarise(Count = sum(Count)) %>%
         pivot_wider(id_cols = Family, names_from = Sample, values_from = Count, values_fill = 0),
       Genus = featuretable_long %>%
         mutate(Genus = str_remove(Taxon, '; s__.*')) %>%
         group_by(Sample, Genus) %>%
         summarise(Count = sum(Count)) %>%
         pivot_wider(id_cols = Genus, names_from = Sample, values_from = Count, values_fill = 0)) %>%
    writexl::write_xlsx('featurecount_ranks.xlsx')
  
  list(Phylum = featuretable_long %>%
         mutate(Phylum = str_remove(Taxon, '; c__.*')) %>%
         group_by(Sample, Phylum) %>%
         summarise(RA = sum(RA)) %>%
         pivot_wider(id_cols = Phylum, names_from = Sample, values_from = RA, values_fill = 0),
       Class = featuretable_long %>%
         mutate(Class = str_remove(Taxon, '; o__.*')) %>%
         group_by(Sample, Class) %>%
         summarise(RA = sum(RA)) %>%
         pivot_wider(id_cols = Class, names_from = Sample, values_from = RA, values_fill = 0),
       Order = featuretable_long %>%
         mutate(Order = str_remove(Taxon, '; f__.*')) %>%
         group_by(Sample, Order) %>%
         summarise(RA = sum(RA)) %>%
         pivot_wider(id_cols = Order, names_from = Sample, values_from = RA, values_fill = 0),
       Family = featuretable_long %>%
         mutate(Family = str_remove(Taxon, '; g__.*')) %>%
         group_by(Sample, Family) %>%
         summarise(RA = sum(RA)) %>%
         pivot_wider(id_cols = Family, names_from = Sample, values_from = RA, values_fill = 0),
       Genus = featuretable_long %>%
         mutate(Genus = str_remove(Taxon, '; s__.*')) %>%
         group_by(Sample, Genus) %>%
         summarise(RA = sum(RA)) %>%
         pivot_wider(id_cols = Genus, names_from = Sample, values_from = RA, values_fill = 0)) %>%
    writexl::write_xlsx('featurerelabund_ranks.xlsx')
}
}


















