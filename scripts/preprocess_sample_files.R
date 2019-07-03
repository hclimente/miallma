library(tidyverse)
library(readxl)

# ALL
phase1 <- read_xlsx('META_TARGET_ALL/TARGET_ALL_PhaseI_patient_table.xlsx') %>%
  filter(!is.na(Project_name))
phase2 <- read_xlsx('META_TARGET_ALL/TARGET_ALL_PhaseII_patient_table.xlsx') %>%
  filter(!is.na(Project_name))

all <- bind_rows(phase1, phase2) %>%
  filter(`DX/RL/PT` != 'DX') %>%
  mutate(tumor = 'ALL',
         type = ifelse(`DX/RL/PT` == 'PT', 'Control', 'Case'),
         sample = Patient_runID) %>%
  select(sample, tumor, type)

# AML
aml <- read_xlsx('META_TARGET_AML/TARGET_AML_patient_table.xlsx') %>%
  filter(!is.na(Project_name)) %>%
  filter(`DX/RL/PT` != 'DX') %>%
  mutate(tumor = 'AML',
         type = ifelse(`DX/RL/PT` == 'PT', 'Control', 'Case'),
         sample = Patient_runID) %>%
  select(sample, tumor, type)

bind_rows(all, aml) %>%
  write_tsv('sample_info.tsv')
