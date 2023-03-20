#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test ont_rna {flames basic SIRV test} {
	file delete -force tmp/sirv
	file mkdir tmp/sirv/fastq
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/sirv/fastq/[file tail $file]
	}
	foreach file [glob -nocomplain data/SIRV-flames/*] {
		if {$file eq "data/SIRV-flames/fastq"} continue
		mklink $file tmp/sirv/[file tail $file]
	}
	# flames works directly form fastq (and we could give the fastq directory instead of the bam file
	# going via bam to fit into the the typical workflow
	cg refseq_minimap2 tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta splice
	cg map \
		-method minimap2 -preset splice -paired 0 \
		tmp/sirv/map-minimap2-sirv.bam \
		tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		tmp/sirv/sirv \
		tmp/sirv/fastq/sample1.fastq.gz tmp/sirv/fastq/sample2.fastq.gz
	exec samtools index tmp/sirv/map-minimap2-sirv.bam
	cg iso_flames -stack 1 \
		-refseq tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf \
		tmp/sirv/map-minimap2-sirv.bam
	# check vs expected
	exec diff tmp/sirv/isoform_counts-flames-fastqs-sirv.tsv data/isoform_counts-flames-fastqs-sirv.tsv
	exec diff tmp/sirv/gene_counts-flames-fastqs-sirv.tsv data/gene_counts-flames-fastqs-sirv.tsv
} {}

test ont_rna {flames empty fastq} {
	file delete -force tmp/sirv
	file mkdir tmp/sirv/fastq
	file_write tmp/sirv/fastq/empty.fastq ""
	foreach file [glob -nocomplain data/SIRV-flames/SIRV*] {
		if {$file eq "data/SIRV-flames/fastq"} continue
		mklink $file tmp/sirv/[file tail $file]
	}
	# flames works directly form fastq (and we could give the fastq directory instead of the bam file
	# going via bam to fit into the the typical workflow
	cg refseq_minimap2 tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta splice
	cg map \
		-method minimap2 -preset splice -paired 0 \
		tmp/sirv/map-minimap2-sirv.bam \
		tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		tmp/sirv/sirv \
		tmp/sirv/fastq/empty.fastq
	exec samtools index tmp/sirv/map-minimap2-sirv.bam
	cg iso_flames -stack 1 \
		-refseq tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf \
		tmp/sirv/map-minimap2-sirv.bam
	# check vs expected
	cg select -overwrite 1 -q 0 data/isoform_counts-flames-fastqs-sirv.tsv tmp/expected-isoform_counts-flames-fastqs-sirv.tsv
	cg select -overwrite 1 -q 0 data/gene_counts-flames-fastqs-sirv.tsv tmp/expected-gene_counts-flames-fastqs-sirv.tsv
	exec diff tmp/sirv/isoform_counts-flames-fastqs-sirv.tsv tmp/expected-isoform_counts-flames-fastqs-sirv.tsv
	exec diff tmp/sirv/gene_counts-flames-fastqs-sirv.tsv tmp/expected-gene_counts-flames-fastqs-sirv.tsv
} {}

test ont_rna {flair basic SIRV test} {
	file delete -force tmp/sirv
	file mkdir tmp/sirv/fastq
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/sirv/fastq/[file tail $file]
	}
	foreach file [glob -nocomplain data/SIRV-flames/*] {
		if {$file eq "data/SIRV-flames/fastq"} continue
		mklink $file tmp/sirv/[file tail $file]
	}
	cg refseq_minimap2 tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta splice
	cg map \
		-method minimap2 -preset splice -paired 0 \
		tmp/sirv/map-minimap2-sirv.bam \
		tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		tmp/sirv/sirv \
		tmp/sirv/fastq/sample1.fastq.gz tmp/sirv/fastq/sample2.fastq.gz
	exec samtools index tmp/sirv/map-minimap2-sirv.bam
	cg flair -stack 1 \
		-refseq tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf \
		tmp/sirv/map-minimap2-sirv.bam
	# check vs expected
	exec diff tmp/sirv/isoform_counts-flair-minimap2-sirv.tsv data/isoform_counts-flair-minimap2-sirv.tsv
	exec diff tmp/sirv/gene_counts-flair-minimap2-sirv.tsv data/gene_counts-flair-minimap2-sirv.tsv
} {}

test ont_rna {isoquant basic SIRV test} {
	file delete -force tmp/sirv
	file mkdir tmp/sirv/fastq
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/sirv/fastq/[file tail $file]
	}
	foreach file [glob -nocomplain data/SIRV-flames/*] {
		if {$file eq "data/SIRV-flames/fastq"} continue
		mklink $file tmp/sirv/[file tail $file]
	}
	exec samtools faidx tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta
	cg refseq_minimap2 tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta splice
	cg map \
		-method minimap2 -preset splice -paired 0 \
		tmp/sirv/map-minimap2-sirv.bam \
		tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		tmp/sirv/sirv \
		tmp/sirv/fastq/sample1.fastq.gz tmp/sirv/fastq/sample2.fastq.gz
	exec samtools index tmp/sirv/map-minimap2-sirv.bam
	cg iso_isoquant -stack 1 \
		-refseq tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf \
		tmp/sirv/map-minimap2-sirv.bam
	# check vs expected
	exec diff tmp/sirv/isoform_counts-isoquant-minimap2-sirv.tsv data/isoform_counts-isoquant-minimap2-sirv.tsv
	exec diff tmp/sirv/gene_counts-isoquant-minimap2-sirv.tsv data/gene_counts-isoquant-minimap2-sirv.tsv
} {}

test ont_rna {flames SIRV test no ref} {
	file delete -force tmp/sirv
	file mkdir tmp/sirv/fastq
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/sirv/fastq/[file tail $file]
	}
	foreach file [glob -nocomplain data/SIRV-flames/*] {
		if {$file eq "data/SIRV-flames/fastq"} continue
		mklink $file tmp/sirv/[file tail $file]
	}
	exec samtools faidx tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta
	cg refseq_minimap2 tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta splice
	cg map \
		-method minimap2 -preset splice -paired 0 \
		tmp/sirv/map-minimap2-sirv.bam \
		tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		tmp/sirv/sirv \
		tmp/sirv/fastq/sample1.fastq.gz tmp/sirv/fastq/sample2.fastq.gz
	exec samtools index tmp/sirv/map-minimap2-sirv.bam
	# does not seem to find new genes, so add one of each
	exec grep SIRV101 tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf > tmp/sirv/part_SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf
	foreach id {SIRV201N SIRV205P SIRV301P SIRV308N SIRV403N SIRV409P SIRV501P SIRV512N SIRV601P SIRV617N SIRV701N} {
		exec grep $id tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf >> tmp/sirv/part_SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf
	}
	cg iso_flames \
		-refseq tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts tmp/sirv/part_SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf \
		tmp/sirv/map-minimap2-sirv.bam
	cg gtf2tsv tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf | cg select -s - -f {* structural_category="known" counts-ref=1} > tmp/sirv/ref.tsv
	file delete tmp/sirv/multitranscript.tsv
	cg multitranscript -match . tmp/sirv/multitranscript.tsv tmp/sirv/isoform_counts-flames-fastqs-sirv.tsv tmp/sirv/ref.tsv 
	# check vs expected
	exec diff tmp/sirv/isoform_counts-flames-fastqs-sirv.tsv data/isoform_counts-flames-fastqs-noref_sirv.tsv
	exec diff tmp/sirv/gene_counts-flames-fastqs-sirv.tsv data/gene_counts-flames-fastqs-noref_sirv.tsv
} {}

test ont_rna {flair SIRV test no ref} {
	file delete -force tmp/sirv
	file mkdir tmp/sirv/fastq
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/sirv/fastq/[file tail $file]
	}
	foreach file [glob -nocomplain data/SIRV-flames/*] {
		if {$file eq "data/SIRV-flames/fastq"} continue
		mklink $file tmp/sirv/[file tail $file]
	}
	exec samtools faidx tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta
	cg refseq_minimap2 tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta splice
	cg map \
		-method minimap2 -preset splice -paired 0 \
		tmp/sirv/map-minimap2-sirv.bam \
		tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		tmp/sirv/sirv \
		tmp/sirv/fastq/sample1.fastq.gz tmp/sirv/fastq/sample2.fastq.gz
	exec samtools index tmp/sirv/map-minimap2-sirv.bam
	exec grep SIRV101 tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf > tmp/sirv/part_SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf
	foreach id {SIRV201N SIRV205P SIRV301P SIRV308N SIRV403N SIRV409P SIRV501P SIRV512N SIRV601P SIRV617N SIRV701N} {
		exec grep $id tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf >> tmp/sirv/part_SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf
	}
	cg flair -stack 1 -v 2 \
		-refseq tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts tmp/sirv/part_SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf \
		tmp/sirv/map-minimap2-sirv.bam
	cg gtf2tsv tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf | cg select -s - -f {* structural_category="known" counts-ref=1} > tmp/sirv/ref.tsv
	file delete tmp/sirv/multitranscript.tsv
	cg multitranscript -match . tmp/sirv/multitranscript.tsv tmp/sirv/isoform_counts-flair-minimap2-sirv.tsv tmp/sirv/ref.tsv 
	# check vs expected
	exec diff tmp/sirv/isoform_counts-flair-minimap2-sirv.tsv data/isoform_counts-flair-noref_sirv.tsv
	exec diff tmp/sirv/gene_counts-flair-minimap2-sirv.tsv data/gene_counts-flair-noref_sirv.tsv
} {}

test ont_rna {isoquant SIRV test no ref} {
	file delete -force tmp/sirv
	file mkdir tmp/sirv/fastq
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/sirv/fastq/[file tail $file]
	}
	foreach file [glob -nocomplain data/SIRV-flames/*] {
		if {$file eq "data/SIRV-flames/fastq"} continue
		mklink $file tmp/sirv/[file tail $file]
	}
	exec samtools faidx tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta
	cg refseq_minimap2 tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta splice
	cg map \
		-method minimap2 -preset splice -paired 0 \
		tmp/sirv/map-minimap2-sirv.bam \
		tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		tmp/sirv/sirv \
		tmp/sirv/fastq/sample1.fastq.gz tmp/sirv/fastq/sample2.fastq.gz
	exec samtools index tmp/sirv/map-minimap2-sirv.bam
	# isoquant will find the "novel" genes so only giving one (not all). It will find a few more isoforms with all genes given though
	cg gtf2tsv tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf | cg select -s - > tmp/sirv/ref.tsv
	cg select -q {$transcript in "SIRV101"} tmp/sirv/ref.tsv > tmp/sirv/part_ref.tsv
	cg iso_isoquant -stack 1 \
		-refseq tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts tmp/sirv/part_ref.tsv \
		tmp/sirv/map-minimap2-sirv.bam
	cg gtf2tsv tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf | cg select -s - -f {* structural_category="known" counts-ref=1} > tmp/sirv/ref.tsv
	file delete tmp/sirv/multitranscript.tsv
	cg multitranscript -match . tmp/sirv/multitranscript.tsv tmp/sirv/isoform_counts-isoquant-minimap2-sirv.tsv tmp/sirv/ref.tsv 
	# check vs expected
	exec diff tmp/sirv/isoform_counts-isoquant-minimap2-sirv.tsv data/isoform_counts-isoquant-noref_sirv.tsv
	exec diff tmp/sirv/gene_counts-isoquant-minimap2-sirv.tsv data/gene_counts-isoquant-noref_sirv.tsv
} {}

test ont_rna {isoquant SIRV test no ref -skipregions} {
	file delete -force tmp/sirv
	file mkdir tmp/sirv/fastq
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/sirv/fastq/[file tail $file]
	}
	foreach file [glob -nocomplain data/SIRV-flames/*] {
		if {$file eq "data/SIRV-flames/fastq"} continue
		mklink $file tmp/sirv/[file tail $file]
	}
	exec samtools faidx tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta
	cg refseq_minimap2 tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta splice
	cg map \
		-method minimap2 -preset splice -paired 0 \
		tmp/sirv/map-minimap2-sirv.bam \
		tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		tmp/sirv/sirv \
		tmp/sirv/fastq/sample1.fastq.gz tmp/sirv/fastq/sample2.fastq.gz
	exec samtools index tmp/sirv/map-minimap2-sirv.bam
	# isoquant will find the "novel" genes so only giving one (not all). It will find a few more isoforms with all genes given though
	cg gtf2tsv tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf | cg select -s - > tmp/sirv/ref.tsv
	cg select -q {$transcript in "SIRV101"} tmp/sirv/ref.tsv > tmp/sirv/part_ref.tsv
	cg iso_isoquant -stack 1 \
		-refseq tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts tmp/sirv/part_ref.tsv \
		-skipregions SIRV7 \
		tmp/sirv/map-minimap2-sirv.bam
	cg gtf2tsv tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf | cg select -s - -f {* structural_category="known" counts-ref=1} > tmp/sirv/ref.tsv
	file delete tmp/sirv/multitranscript.tsv
	cg multitranscript -match . tmp/sirv/multitranscript.tsv tmp/sirv/isoform_counts-isoquant-minimap2-sirv.tsv tmp/sirv/ref.tsv 
	# check vs expected
	exec diff tmp/sirv/isoform_counts-isoquant-minimap2-sirv.tsv data/isoform_counts-isoquant-noref_sirv.tsv
} {54a55,57
> SIRV7	1000	147946	-	1000,2993,3809,114680,147608	2675,3111,3896,114988,147946	novelt_SIRV7_1000-e1675i318e118i698e87i110784e308i32620e338	novelg_SIRV7_m_1001_147947	novelg_SIRV7_m_1001_147947	transcript				5	IsoQuant		transcript6.SIRV7.nnic	novel_gene	2526	16.50	16.50	6	5	0.00	0	0
> SIRV7	1000	147946	-	1000,2993,43028,114680,147608	2675,3111,43077,114988,147946	novelt_SIRV7_1000-e1675i318e118i39917e49i71603e308i32620e338	novelg_SIRV7_m_1001_147947	novelg_SIRV7_m_1001_147947	transcript				5	IsoQuant		transcript8.SIRV7.nnic	novel_gene	2488	23.50	23.50	13	12	0.00	0	0
> SIRV7	56033	147947	-	56033,70883,78841,114680,147608	56097,70987,78965,114960,147947	novelt_SIRV7_56033-e64i14786e104i7854e124i35715e280i32648e339	novelg_SIRV7_m_1001_147947	novelg_SIRV7_m_1001_147947	transcript				5	IsoQuant		transcript11.SIRV7.nnic	novel_gene	911	26.50	50.00	50	37	0.00	0	0
child process exited abnormally} error

test ont_rna {isoquant_sens SIRV test no ref} {
	file delete -force tmp/sirv
	file mkdir tmp/sirv/fastq
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/sirv/fastq/[file tail $file]
	}
	foreach file [glob -nocomplain data/SIRV-flames/*] {
		if {$file eq "data/SIRV-flames/fastq"} continue
		mklink $file tmp/sirv/[file tail $file]
	}
	exec samtools faidx tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta
	cg refseq_minimap2 tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta splice
	cg map \
		-method minimap2 -preset splice -paired 0 \
		tmp/sirv/map-minimap2-sirv.bam \
		tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		tmp/sirv/sirv \
		tmp/sirv/fastq/sample1.fastq.gz tmp/sirv/fastq/sample2.fastq.gz
	exec samtools index tmp/sirv/map-minimap2-sirv.bam
	# isoquant will find the "novel" genes so only giving one (not all). It will find a few more isoforms with all genes given though
	cg gtf2tsv tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf | cg select -s - > tmp/sirv/ref.tsv
	cg select -q {$transcript in "SIRV101"} tmp/sirv/ref.tsv > tmp/sirv/part_ref.tsv
	cg iso_isoquant -stack 1 \
		-preset sens \
		-refseq tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts tmp/sirv/part_ref.tsv \
		tmp/sirv/map-minimap2-sirv.bam
	cg gtf2tsv tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf | cg select -s - -f {* structural_category="known" counts-ref=1} > tmp/sirv/ref.tsv
	file delete tmp/sirv/multitranscript.tsv
	cg multitranscript -match . tmp/sirv/multitranscript.tsv tmp/sirv/isoform_counts-isoquant_sens-minimap2-sirv.tsv tmp/sirv/ref.tsv 
	# check vs expected
	exec diff tmp/sirv/isoform_counts-isoquant_sens-minimap2-sirv.tsv data/isoform_counts-isoquant_sens-noref_sirv.tsv
	exec diff tmp/sirv/gene_counts-isoquant_sens-minimap2-sirv.tsv data/gene_counts-isoquant_sens-noref_sirv.tsv
} {}

test ont_rna {isoquant_all SIRV test no ref} {
	file delete -force tmp/sirv
	file mkdir tmp/sirv/fastq
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/sirv/fastq/[file tail $file]
	}
	foreach file [glob -nocomplain data/SIRV-flames/*] {
		if {$file eq "data/SIRV-flames/fastq"} continue
		mklink $file tmp/sirv/[file tail $file]
	}
	exec samtools faidx tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta
	cg refseq_minimap2 tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta splice
	cg map \
		-method minimap2 -preset splice -paired 0 \
		tmp/sirv/map-minimap2-sirv.bam \
		tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		tmp/sirv/sirv \
		tmp/sirv/fastq/sample1.fastq.gz tmp/sirv/fastq/sample2.fastq.gz
	exec samtools index tmp/sirv/map-minimap2-sirv.bam
	# isoquant will find the "novel" genes so only giving one (not all). It will find a few more isoforms with all genes given though
	cg gtf2tsv tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf | cg select -s - > tmp/sirv/ref.tsv
	cg select -q {$transcript in "SIRV101"} tmp/sirv/ref.tsv > tmp/sirv/part_ref.tsv
	cg iso_isoquant -stack 1 \
		-preset all \
		-refseq tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts tmp/sirv/part_ref.tsv \
		tmp/sirv/map-minimap2-sirv.bam
	cg gtf2tsv tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf | cg select -s - -f {* structural_category="known" counts-ref=1} > tmp/sirv/ref.tsv
	file delete tmp/sirv/multitranscript.tsv
	cg multitranscript -match . tmp/sirv/multitranscript.tsv tmp/sirv/isoform_counts-isoquant_all-minimap2-sirv.tsv tmp/sirv/ref.tsv 
	# check vs expected
	exec diff tmp/sirv/isoform_counts-isoquant_all-minimap2-sirv.tsv data/isoform_counts-isoquant_all-noref_sirv.tsv
	exec diff tmp/sirv/gene_counts-isoquant_all-minimap2-sirv.tsv data/gene_counts-isoquant_all-noref_sirv.tsv
} {}

test ont_rna {process_project multi methods} {
	test_cleantmp
	file mkdir tmp/samples/sirv1/fastq
	file mkdir tmp/samples/sirv2/fastq
	file mkdir tmp/ref/sirv
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/samples/sirv1/fastq/[file tail $file]
		exec cg fastq2tsv $file | cg select -q {$ROW < 1000} | cg tsv2fastq | cg bgzip > tmp/samples/sirv2/fastq/[file tail $file]
	}
	mklink data/SIRV-flames/SIRV_isoforms_multi-fasta_170612a.fasta tmp/ref/sirv/genome_sirv.ifas
	mklink data/SIRV-flames/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf tmp/ref/sirv/gene_sirv.gtf
	exec samtools faidx tmp/ref/sirv/genome_sirv.ifas
	cg refseq_minimap2 tmp/ref/sirv/genome_sirv.ifas splice
	file delete tmp/compar/isoform_counts-tmp.tsv
	exec cg process_project -stack 1 -v 2 \
		-split 1 \
		-d 4 \
		-threads 2 \
		-paired 0 -clip 0 \
		-maxfastqdistr 250 \
		-aligner {minimap2_splice} \
		-removeduplicates 0 \
		-realign 0 \
		-distrreg 0 \
		-svcallers {} \
		-varcallers {} \
		-isocallers {isoquant flair flames} \
		-iso_match . \
		-reports {} \
		-dbdir tmp/ref/sirv \
		tmp \
		>& tmp/ontrna.log
	# check vs expected
	exec diff tmp/compar/isoform_counts-tmp.tsv data/ontrna/isoform_counts-tmp.tsv
	exec diff tmp/compar/gene_counts-tmp.tsv data/ontrna/gene_counts-tmp.tsv
} {}

testsummarize
