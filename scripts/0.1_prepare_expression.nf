#!/usr/bin/env nextflow

params.data = "../data/"
data = params.data

tx_expression = file("${data}target_phases1_2.tsv")
sample_info = file("${data}sample_info.tsv")

tumors = Channel
		.from(sample_info)
		.splitCsv(header: true)
		.map { row -> row.tumor }
		.unique()
origins = ["Control", "Case"]

process get_gtf {

	output:
		file 'gtf' into gtf

	'''
	wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_23/gencode.v23.annotation.gtf.gz
	gunzip -c gencode.v23.annotation.gtf.gz >gtf
	'''

}

process get_tx2gene {

	publishDir data, overwrite: true, mode: "copy"

	input:
		file gtf

	output:
		file "tx2gene.csv" into tx2gene

	"""
	echo "gene,transcript" >tx2gene.csv
	awk '\$3 == "transcript"' $gtf | cut -f9 | sed 's/gene_id "//' | sed 's/"; transcript_id "/,/' | sed -E 's/";.+//' >>tx2gene.csv
	"""

}

process sort_expression {

	input:
		each i from 1..10
		each j from 1..10
		file tx_expression
		file tx2gene

	output:
		file 'sorted_bit' into sortedBits

	"""
	grep ^ENST00000$i$j $tx_expression | sort >expression
	grep ,ENST00000$i$j $tx2gene | sort --key=2 -t, | cut -f1 -d, >genes
	paste genes expression | sort | cut -f2- >>sorted_bit
	"""

}

process join_expression_bits {

	input:
		file tx_expression
		file tx2gene
		file 'sorted_bit*' from sortedBits.collect()

	output:
		file 'sorted' into sortedTxExpression

	"""
	head -n1 $tx_expression >sorted
	grep ^ENSTR0000 $tx_expression | sort >expression
	grep ,ENSTR0000 $tx2gene | sort --key=2 -t, | cut -f1 -d, >genes
	paste genes expression | sort | cut -f2- >>sorted
	cat sorted_bit* >>sorted
	"""

}

process get_tumor_expression {

	publishDir "${data}/expression", overwrite: true, mode: "copy"

	input:
		each tumor from tumors
		each origin from origins
		file sample_info 
		file sortedTxExpression

	output:
		file "${tumor}_${origin}_isoform_tpm.tsv" into tumor_tx_expression

	"""
	#!/usr/bin/env Rscript

	library(tidyverse)

	tumorSamples <- read_tsv("$sample_info") %>%
		filter(tumor == "$tumor" & type == "$origin") %>%
		.\$sample

	allSamples <- read_tsv("$tx_expression", n_max = 1) %>% colnames

	cols <- which(allSamples %in% tumorSamples) %>% paste(collapse = ",") %>% paste0("-f1,", .)

	system2("cut", args = c(cols, "$tx_expression"), stdout = "tmp")
	system2('sed', args = c("'s/^sample/transcript/'", 'tmp'), stdout = "${tumor}_${origin}_isoform_tpm.tsv")
	"""

}
