library('tidyverse')
library('rstatix')
library('vegan')

# importing metadata
metadata <- readxl::read_xlsx('metadata.xlsx') %>%
  mutate(TreatmentSex = paste(Treatment, Sex, sep = '_'))

# plotting read count
reads <- readxl::read_xlsx('numreads.xlsx') %>%
  left_join(metadata)

reads %>%
  writexl::write_xlsx('readcounts.xlsx')

# importing Emu profiles
df <- read_tsv('emu-combined-species-counts.tsv')

df <- df %>%
  dplyr::select(species, starts_with('DMG')) %>%
  rename(Species = species)
df[is.na(df)] <- 0

# removing samples with low read counts
df <- df %>% 
  dplyr::select(-DMG2303277)
metadata <- metadata %>%
  filter(Sample != 'DMG2303277')

# have already seen that there is a large variation by cohort
# some samples spent a long time in the freezer, so going to try to account for this with PLSDAbatch
library(PLSDAbatch)
library(mixOmics)

df <- df %>%
  column_to_rownames('Species') %>%
  t()

df_filtering <- PreFL(df, keep.spl = 10, keep.var = 0.01)
df_filtered <- df_filtering$data.filter

dim(df)
dim(df_filtered)

df_clr <- logratio.transfo(X = df_filtered, logratio = 'CLR', offset = 1)

df_pca_before <- pca(df_clr, ncomp = 3, scale = TRUE)
Scatter_Density(object = df_pca_before, batch = metadata$Cohort, trt = metadata$Treatment, 
                batch.legend.title = 'Cohort', trt.legend.title = 'Treatment')

df_treatment_tune <- plsda(X = df_clr, Y = metadata$Treatment, ncomp = 5)
df_treatment_tune$prop_expl_var

df_batch_tune <- PLSDA_batch(X = df_clr, 
                             Y.trt = metadata$Treatment, Y.bat = metadata$Cohort,
                             ncomp.trt = 2, ncomp.bat = 5)
df_batch_tune$explained_variance.bat 

df_batch_res <- PLSDA_batch(X = df_clr, 
                            Y.trt = metadata$Treatment, Y.bat = metadata$Cohort,
                            ncomp.trt = 1, ncomp.bat = 1)
df_batch <- df_batch_res$X.nobatch

df_pca_after <- pca(df_batch, ncomp = 3, scale = TRUE)
Scatter_Density(object = df_pca_after, batch = metadata$Cohort, trt = metadata$Treatment, 
                batch.legend.title = 'Cohort', trt.legend.title = 'Treatment')




