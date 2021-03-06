#!/usr/bin/env nextflow

params.data = "../data"
data = params.data
params.out = "./"
out = params.out

network = file("$data/gencode_v27.pklz")
ddi = file("$data/ddi.tsv")
drivers = file("$data/mutational_drivers.tsv")
expression = file("$data/expression")
d1 = file("$data/d1.txt")

// COMPUTE TRANS SWITCHES
////////////////////////////////////////
tumors = ['ALL']
projects = ['TARGET_other','TARGET_phase1','TARGET_phase2_dc','TARGET_phase2_di', 'TARGET_phase2_vl']

process calculate_trans_switches {

	publishDir "${out}/${tumor}_${project}", overwrite: true, mode: "copy"

	input:
		each tumor from tumors
		each project from projects
		file "annotation_tumor.pklz" from network
		file expression

	output:
		set tumor, project,'switches_spada.tsv',"annotation_${tumor}.pklz" into switches

	"""
	spada init --name ${tumor} --annotation annotation_tumor.pklz
	spada switches --minimum-expression 0.1 \
--expression-control $expression/${tumor}_${project}_DX_isoform_tpm.tsv \
--expression-case $expression/${tumor}_${project}_RL_isoform_tpm.tsv
	mv annotation.pklz annotation_${tumor}.pklz
	"""

}

// ANALYZE SWITCHES
////////////////////////////////////////
switches .into { functional; summary }

process functional_analysis {

	memory = '15GB'
	publishDir "${out}/${tumor}_${project}", overwrite: true, mode: "copy"

	input:
		set val(tumor), val(project), file('switches.tsv'), file('annotation.pklz') from functional

	output:
		file 'pfam_analysis.tsv'
		file 'prosite_analysis.tsv'
		file 'ppi_analysis.tsv'

	"""
	spada function --switches switches.tsv
	"""

}

process summary {

	publishDir "${out}/${tumor}_${project}", overwrite: true, mode: "copy"

	input:
		set val(tumor), val(project), file('switches.tsv'), file('annotation.pklz') from summary
		file expression

	output:
		set tumor, project, 'switches_spada.tsv' into spada_switches
		file 'proteome_features.tsv'

	"""
	spada summary --switches switches.tsv --minimum-expression 0.1 \
--expression-control $expression/${tumor}_${project}_DX_isoform_tpm.tsv \
--expression-case $expression/${tumor}_${project}_RL_isoform_tpm.tsv
	"""

}

// SWITCH POST-PROCESSING
////////////////////////////////////////
process driver_annotation {

	publishDir "${out}/${tumor}_${project}", overwrite: true, mode: "copy"

	input:
		set val(tumor), val(project), file('switches.tsv') from spada_switches
		file d1
		file drivers

	output:
		file 'switches_final.tsv'

	"""
	#!/usr/bin/env Rscript

	library(tidyverse)

	drivers <- read_tsv('$drivers', col_types = 'cc') %>%
	    mutate(Driver = 'Tumor-specific driver')

	d1 <- read_tsv('$d1', col_types = 'c')

	read_tsv('switches.tsv', col_types = 'ccccciiiic') %>%
		left_join(drivers, by = c('Symbol' = 'Gene', 'Experiment' = 'Tumor_type')) %>%
		mutate(Driver = ifelse(is.na(Driver) & Symbol %in% d1\$Symbol, 'Driver interactor', Driver),
			   Driver = ifelse(is.na(Driver) & Symbol %in% drivers\$Gene, 'Foreign driver', Driver),
		       Driver = ifelse(is.na(Driver), 'No', Driver)) %>%
		write_tsv('switches_final.tsv')
	"""

}
