#!/usr/bin/env nextflow

params.data = "../data/"
data = params.data

// PREPARE INTOGEN DRIVERS
////////////////////////////////////////
cancerGenes = file("${data}CancerGenes_Pediatric_Adult.xlsx")
tumorDrivers = file("${data}Mutational_drivers_per_tumor_type.tsv")
matching = file("${data}intogen2tcga.csv")

process get_mutational_drivers {

	publishDir data, overwrite: true, mode: "copy"

	input:
		file cancerGenes
		file tumorDrivers
		file matching

	output:
		file "mutational_drivers.tsv" into drivers

	"""
	#!/usr/bin/env Rscript

	library(tidyverse)
	library(readxl)

	matching <- read_csv("$matching")
	drivers <- matching\$TCGA
	names(drivers) <- matching\$Intogen

	# convert intogen tumor names into tcga
	tumorDrivers <- read_tsv("$tumorDrivers", comment="#") %>%
		mutate(Tumor_type = ifelse(Tumor_type %in% names(drivers), drivers[Tumor_type], Tumor_type)) %>%
		separate_rows(Tumor_type)

	# annotate cancer genes
	read_excel("$cancerGenes") %>%
		select(Gene) %>%
		left_join(tumorDrivers, by = c('Gene' = 'geneHGNCsymbol')) %>%
		mutate(Tumor_type = ifelse(is.na(Tumor_type), '', Tumor_type)) %>%
		write_tsv("mutational_drivers.tsv")
	"""

}

process get_ppi {

  output:
    file 'BIOGRID-MV-Physical-*mitab.txt' into mitab

  """
  wget https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-MV-Physical-LATEST.mitab.zip
  unzip BIOGRID-MV-Physical-LATEST.mitab.zip
  """

}

process get_d1 {

	publishDir data, overwrite: true, mode: "copy"

	input:
		file drivers
		file mitab

	output:
		file 'd1.txt' into d1

	"""
	#!/usr/bin/env Rscript

	library(tidyverse)

	drivers <- read_tsv('$drivers', col_types = 'cc') %>%
    select(Gene) %>%
    unique %>%
    mutate(driver = TRUE)


	read_tsv('$mitab', col_types = 'ccccccccccccccc') %>%
		filter(`Taxid Interactor A` == 'taxid:9606' & `Taxid Interactor B` == 'taxid:9606') %>%
		mutate(gene_A = gsub('.+locuslink:', '', `Alt IDs Interactor A`),
					 gene_B = gsub('.+locuslink:', '', `Alt IDs Interactor B`)) %>%
		select(starts_with('gene')) %>%
		left_join(drivers, by = c('gene_A' = 'Gene')) %>%
    left_join(drivers, by = c('gene_B' = 'Gene')) %>%
    mutate(driver.x = ! is.na(driver.x), driver.y = ! is.na(driver.y)) %>%
    filter(driver.x != driver.y) %>%
    gather(key = 'pos', value = 'Symbol', starts_with('gene')) %>%
    mutate(d1 = (driver.x & pos == 'gene_B') | (driver.y & pos == 'gene_A')) %>%
    filter(d1) %>%
    select(Symbol) %>%
    unique %>%
		write_tsv('d1.txt')
	"""

}