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
		-ali_keepcomments 0 \
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
		-ali_keepcomments 0 \
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

test ont_rna {flair basic SIRV test resultfile} {
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
		-ali_keepcomments 0 \
		tmp/sirv/map-minimap2-sirv.bam \
		tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		tmp/sirv/sirv \
		tmp/sirv/fastq/sample1.fastq.gz tmp/sirv/fastq/sample2.fastq.gz
	exec samtools index tmp/sirv/map-minimap2-sirv.bam
	cg flair -stack 1 \
		-refseq tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta \
		-reftranscripts tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf \
		tmp/sirv/map-minimap2-sirv.bam tmp/result/iso_count-flair-result.tsv
	# check vs expected
	cg select -overwrite 1 -f {transcript gene geneid chromosome strand begin end exonStarts exonEnds cdsStart cdsEnd exonCount type transcripttype counts-flair-result=$counts-flair-minimap2-sirv} \
		data/isoform_counts-flair-minimap2-sirv.tsv expected.tsv 
	exec diff tmp/result/iso_count-flair-result.tsv expected.tsv
	cg select -overwrite 1 -f {type	gene	gene_type	chromosome	begin	end	strand	nrtranscripts	counts-flair-result=$counts-flair-minimap2-sirv} \
		data/gene_counts-flair-minimap2-sirv.tsv expected.tsv 
	exec diff tmp/result/gene_counts-flair-result.tsv expected.tsv 
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
		-ali_keepcomments 0 \
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
		-ali_keepcomments 0 \
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
		-ali_keepcomments 0 \
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
		-ali_keepcomments 0 \
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
		-ali_keepcomments 0 \
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
> SIRV7	1000	147946	-	1000,2993,3809,114680,147608	2675,3111,3896,114988,147946	novelt_SIRV7_1000-e1675i318e118i698e87i110784e308i32620e338	novelg_SIRV7_m_1001_147946	novelg_SIRV7_m_1001_147946	transcript				5	IsoQuant		transcript1.SIRV7.nnic	novel_gene	2526	16.50	16.5	6	5	6.5	6	5	0	0	0
> SIRV7	1000	147946	-	1000,2993,43028,114680,147608	2675,3111,43077,114988,147946	novelt_SIRV7_1000-e1675i318e118i39917e49i71603e308i32620e338	novelg_SIRV7_m_1001_147946	novelg_SIRV7_m_1001_147946	transcript				5	IsoQuant		transcript3.SIRV7.nnic	novel_gene	2488	23.50	23.5	13	12	13.5	13	12	0	0	0
> SIRV7	56033	147947	-	56033,70883,78841,114680,147608	56097,70987,78965,114960,147947	novelt_SIRV7_56033-e64i14786e104i7854e124i35715e280i32648e339	novelg_SIRV7_m_56034_147947	novelg_SIRV7_m_56034_147947	transcript				5	IsoQuant		transcript6.SIRV7.nnic	novel_gene	911	29.00	29	29	19	29	29	19	0	0	0
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
		-ali_keepcomments 0 \
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

test ont_rna {flair basic SIRV test -compar joint} {
	test_cleantmp
	file mkdir tmp/samples/sirv1/fastq
	file mkdir tmp/samples/sirv2/fastq
	file mkdir tmp/ref/sirv
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/samples/sirv1/fastq/[file tail $file]
		exec cg fastq2tsv $file | cg select -q {$ROW < 1000} | cg tsv2fastq | cg bgzip > tmp/samples/sirv2/fastq/[file tail $file]
	}
	mkdir tmp/ref
	foreach file [glob -nocomplain data/SIRV-flames/*] {
		if {$file eq "data/SIRV-flames/fastq"} continue
		mklink $file tmp/ref/[file tail $file]
	}
	mklink data/SIRV-flames/SIRV_isoforms_multi-fasta_170612a.fasta tmp/ref/sirv/genome_sirv.ifas
	mklink data/SIRV-flames/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf tmp/ref/sirv/gene_sirv.gtf
	exec samtools faidx tmp/ref/sirv/genome_sirv.ifas
	cg refseq_minimap2 tmp/ref/sirv/genome_sirv.ifas splice
	file delete tmp/compar/isoform_counts-tmp.tsv
	cg map \
		-method minimap2 -preset splice -paired 0 \
		-ali_keepcomments 0 \
		tmp/samples/sirv1/map-minimap2-sirv1.bam \
		tmp/ref/sirv/genome_sirv.ifas \
		tmp/samples/sirv1/sirv1 \
		tmp/samples/sirv1/fastq/sample1.fastq.gz tmp/samples/sirv1/fastq/sample2.fastq.gz
	exec samtools index tmp/samples/sirv1/map-minimap2-sirv1.bam
	cg map \
		-method minimap2 -preset splice -paired 0 \
		-ali_keepcomments 0 \
		tmp/samples/sirv2/map-minimap2-sirv2.bam \
		tmp/ref/sirv/genome_sirv.ifas \
		tmp/samples/sirv2/sirv2 \
		tmp/samples/sirv2/fastq/sample1.fastq.gz tmp/samples/sirv2/fastq/sample2.fastq.gz
	exec samtools index tmp/samples/sirv2/map-minimap2-sirv2.bam
	cg flair -stack 1 -compar joint \
		-refseq tmp/ref/sirv/genome_sirv.ifas \
		-reftranscripts tmp/ref/sirv/gene_sirv.gtf \
		tmp
	# check vs expected
	exec diff tmp/compar/isoform_counts-flair-tmp.tsv data/isoform_counts-flair-tmp.tsv
	exec diff tmp/compar/gene_counts-flair-tmp.tsv data/gene_counts-flair-tmp.tsv
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
	cg tsvdiff tmp/samples/sirv1/isoform_counts-isoquant-sminimap2_splice-sirv1.tsv tmp/ref/sirv/gene_sirv.tsv
} {diff tmp/samples/sirv1/isoform_counts-isoquant-sminimap2_splice-sirv1.tsv tmp/ref/sirv/gene_sirv.tsv
header diff
<extrafields: geneid gene_ori category size counts_iqall-isoquant-sminimap2_splice-sirv1 counts_weighed-isoquant-sminimap2_splice-sirv1 counts_unique-isoquant-sminimap2_splice-sirv1 counts_strict-isoquant-sminimap2_splice-sirv1 counts_sweighed-isoquant-sminimap2_splice-sirv1 counts_sunique-isoquant-sminimap2_splice-sirv1 counts_sstrict-isoquant-sminimap2_splice-sirv1 counts_aweighed-isoquant-sminimap2_splice-sirv1 counts_aunique-isoquant-sminimap2_splice-sirv1 counts_astrict-isoquant-sminimap2_splice-sirv1
---
>extrafields: gene_name gene_id

child process exited abnormally} error

test ont_rna {process_project multi methods using -preset} {
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
		-preset ontr \
		-iso_joint {} \
		-d 4 \
		-threads 2 \
		-isocallers {isoquant flair flames} \
		-iso_match . \
		-reports {} \
		-dbdir tmp/ref/sirv \
		tmp \
		>& tmp/ontrna.log
	# check vs expected
	exec diff tmp/compar/isoform_counts-tmp.tsv data/ontrna/isoform_counts-tmp.tsv
	exec diff tmp/compar/gene_counts-tmp.tsv data/ontrna/gene_counts-tmp.tsv
	cg tsvdiff tmp/samples/sirv1/isoform_counts-isoquant-sminimap2_splice-sirv1.tsv tmp/ref/sirv/gene_sirv.tsv
} {diff tmp/samples/sirv1/isoform_counts-isoquant-sminimap2_splice-sirv1.tsv tmp/ref/sirv/gene_sirv.tsv
header diff
<extrafields: geneid gene_ori category size counts_iqall-isoquant-sminimap2_splice-sirv1 counts_weighed-isoquant-sminimap2_splice-sirv1 counts_unique-isoquant-sminimap2_splice-sirv1 counts_strict-isoquant-sminimap2_splice-sirv1 counts_sweighed-isoquant-sminimap2_splice-sirv1 counts_sunique-isoquant-sminimap2_splice-sirv1 counts_sstrict-isoquant-sminimap2_splice-sirv1 counts_aweighed-isoquant-sminimap2_splice-sirv1 counts_aunique-isoquant-sminimap2_splice-sirv1 counts_astrict-isoquant-sminimap2_splice-sirv1
---
>extrafields: gene_name gene_id

child process exited abnormally} error

test ont_rna {isoquant joint analysis} {
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
	exec cg process_project -stack 1 -v 2 -d 4 \
		-split 1 \
		-threads 2 \
		-paired 0 -clip 0 \
		-maxfastqdistr 250 \
		-aligner {minimap2_splice} \
		-removeduplicates 0 \
		-realign 0 \
		-distrreg 0 \
		-svcallers {} \
		-varcallers {} \
		-isocallers {isoquant} \
		-iso_match . \
		-iso_joint {isoquant} \
		-reports {} \
		-dbdir tmp/ref/sirv \
		tmp \
		>& tmp/ontrna.log
	# check vs expected
	exec diff tmp/compar/isoform_counts-isoquant_joint-tmp.tsv data/ontrna/isoform_counts-isoquant_joint-tmp.tsv
	exec diff tmp/compar/gene_counts-isoquant_joint-tmp.tsv data/ontrna/gene_counts-isoquant_joint-tmp.tsv
} {}

test ont_rna {isoquant joint analysis no ref} {
	test_cleantmp
	file mkdir tmp/samples/sirv1/fastq
	file mkdir tmp/samples/sirv2/fastq
	file mkdir tmp/ref/sirv
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/samples/sirv1/fastq/[file tail $file]
		exec cg fastq2tsv $file | cg select -q {$ROW < 1000} | cg tsv2fastq | cg bgzip > tmp/samples/sirv2/fastq/[file tail $file]
	}
	mklink data/SIRV-flames/SIRV_isoforms_multi-fasta_170612a.fasta tmp/ref/sirv/genome_sirv.ifas
	file delete tmp/ref/sirv/gene_sirv.gtf
	cg gtf2tsv data/SIRV-flames/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf | cg select -s - -q {$transcript in ""} > tmp/ref/sirv/gene_sirv.tsv
	exec samtools faidx tmp/ref/sirv/genome_sirv.ifas
	cg refseq_minimap2 tmp/ref/sirv/genome_sirv.ifas splice
	file delete tmp/compar/isoform_counts-tmp.tsv
	exec cg process_project -stack 1 -v 2 -d 4 \
		-split 1 \
		-threads 2 \
		-paired 0 -clip 0 \
		-maxfastqdistr 250 \
		-aligner {minimap2_splice} \
		-removeduplicates 0 \
		-realign 0 \
		-distrreg 0 \
		-svcallers {} \
		-varcallers {} \
		-isocallers {isoquant} \
		-iso_match . \
		-iso_joint {isoquant} \
		-reports {} \
		-dbdir tmp/ref/sirv \
		tmp \
		>& tmp/ontrna.log
	# check vs expected
	exec diff tmp/compar/isoform_counts-isoquant_joint-tmp.tsv data/ontrna/isoform_counts-isoquant_joint-tmp_noref.tsv
	exec diff tmp/compar/gene_counts-isoquant_joint-tmp.tsv data/ontrna/gene_counts-isoquant_joint-tmp_noref.tsv
} {}

test ont_rna {isoquant SIRV test overlap distrreg borders} {
	file delete -force tmp/sirv
	file mkdir tmp/sirv/fastq
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/sirv/fastq/[file tail $file]
	}
	foreach file [glob -nocomplain data/SIRV-flames/*] {
		if {$file eq "data/SIRV-flames/fastq"} continue
		mklink $file tmp/sirv/[file tail $file]
	}
	file copy tmp/sirv/SIRV_isoforms_multi-fasta_170612a.fasta tmp/sirv/genome_sirv.ifas
	exec samtools faidx tmp/sirv/genome_sirv.ifas
	cg refseq_minimap2 tmp/sirv/genome_sirv.ifas splice
	mkdir tmp/sirv/extra
	file_write tmp/sirv/extra/reg_sirv_distrg.tsv [deindent {
		chromsomoe	begin	end
		SIRV1	0	7000
		SIRV1	7000	12643
		SIRV2	0	6911
		SIRV3	0	10943
		SIRV4	0	16122
		SIRV5	0	14606
		SIRV6	0	12837
		SIRV7	0	148957
	}]\n
	cg map \
		-method minimap2 -preset splice -paired 0 \
		-ali_keepcomments 0 \
		tmp/sirv/map-minimap2-sirv.bam \
		tmp/sirv/genome_sirv.ifas \
		tmp/sirv/sirv \
		tmp/sirv/fastq/sample1.fastq.gz tmp/sirv/fastq/sample2.fastq.gz
	exec samtools index tmp/sirv/map-minimap2-sirv.bam
	cg iso_isoquant -stack 1 \
		-refseq tmp/sirv \
		-reftranscripts tmp/sirv/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf \
		-distrreg g \
		tmp/sirv/map-minimap2-sirv.bam
	# check vs expected
	exec diff tmp/sirv/isoform_counts-isoquant-minimap2-sirv.tsv data/isoform_counts-isoquant-minimap2-sirv-regoverlap.tsv
	exec diff tmp/sirv/gene_counts-isoquant-minimap2-sirv.tsv data/gene_counts-isoquant-minimap2-sirv-regoverlap.tsv
} {}

test ont_rna {isoquant SIRV test no ref overlap distrreg borders} {
	test_cleantmp
	file mkdir tmp/samples/sirv1/fastq
	file mkdir tmp/samples/sirv2/fastq
	file mkdir tmp/ref/sirv
	foreach file [glob -nocomplain data/SIRV-flames/fastq/*] {
		mklink $file tmp/samples/sirv1/fastq/[file tail $file]
		exec cg fastq2tsv $file | cg select -q {$ROW < 1000} | cg tsv2fastq | cg bgzip > tmp/samples/sirv2/fastq/[file tail $file]
	}
	mklink data/SIRV-flames/SIRV_isoforms_multi-fasta_170612a.fasta tmp/ref/sirv/genome_sirv.ifas
	exec samtools faidx tmp/ref/sirv/genome_sirv.ifas
	cg refseq_minimap2 tmp/ref/sirv/genome_sirv.ifas splice
	file delete tmp/ref/sirv/gene_sirv.gtf
	cg gtf2tsv data/SIRV-flames/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf | cg select -s - -q {$transcript in ""} > tmp/ref/sirv/gene_sirv.tsv
	mkdir tmp/ref/sirv/extra
	file_write tmp/ref/sirv/extra/reg_sirv_distrg.tsv [deindent {
		chromsomoe	begin	end
		SIRV1	0	7000
		SIRV1	7000	12643
		SIRV2	0	6911
		SIRV3	0	10943
		SIRV4	0	16122
		SIRV5	0	14606
		SIRV6	0	12837
		SIRV7	0	148957
	}]\n
	file delete tmp/compar/isoform_counts-tmp.tsv
	exec cg process_project -stack 1 -v 2 -d 4 \
		-split 1 \
		-threads 2 \
		-paired 0 -clip 0 \
		-maxfastqdistr 250 \
		-aligner {minimap2_splice} \
		-removeduplicates 0 \
		-realign 0 \
		-svcallers {} \
		-varcallers {} \
		-isocallers {isoquant} \
		-iso_match . \
		-iso_joint {isoquant} \
		-reports {} \
		-dbdir tmp/ref/sirv \
		-distrreg g \
		tmp \
		>& tmp/ontrna.log
	# check vs expected
	exec diff tmp/compar/isoform_counts-isoquant_joint-tmp.tsv data/ontrna/isoform_counts-isoquant_joint-tmp_noref-regoverlap.tsv
	exec diff tmp/compar/gene_counts-isoquant_joint-tmp.tsv data/ontrna/gene_counts-isoquant_joint-tmp_noref-regoverlap.tsv
} {}

test ont_rna {isoquant gene_name_check} {
	array set checkisosa {ENSG00000104228.13 {{27284885,27289161,27290155,27294079,27298463,27310800, 27288127,27289280,27290178,27294310,27298559,27311272,} {27285808,27288035,27289161,27290155,27294079, 27286137,27288127,27289280,27290178,27294135,} {27287700,27289161,27294079,27310800, 27288127,27289280,27294310,27311238,}} ENSG00000015592.17 {{27235307,27239970,27241053,27241676,27242396,27243710,27258350, 27236905,27240162,27241262,27241757,27242492,27243801,27258404,} {27235884,27239225,27239970,27241053,27242396,27258350, 27236905,27239279,27240162,27241262,27242492,27258401,} {27236302,27239970,27241053,27242396,27243710,27258350, 27236905,27240162,27241262,27242492,27243801,27258409,} {27236474,27239225,27239970,27241053,27241676,27242396,27243710,27258350, 27236905,27239279,27240162,27241262,27241757,27242492,27243801,27258420,} {27239339,27241053,27242396,27243710,27258350, 27240162,27241262,27242492,27243801,27258404,} {27239689,27241053,27241676,27242396,27243710,27258350, 27240162,27241262,27241757,27242492,27243801,27258404,} {27241167,27241676,27242109,27243710,27258350, 27241262,27241757,27242492,27243801,27258398,}} ENSG00000253888.2 {{27171913,27184179,27206454,27210041, 27172042,27184387,27206519,27210783,}}}
	array set checka {chr8,+ {{27171913 27210783 ENSG00000253888.2 ENSG00000253888}} chr8,- {{27235307 27258420 ENSG00000015592.17} {27284885 27311272 ENSG00000104228.13 ENSG00000104228}}}
	foreach {name chr strand begin end expected} {
		first chr8 - 27236303 27241544 ENSG00000015592.17
		second chr8 - 27238303 27311238 ENSG00000104228.13
		novel chr8 - 27171000 27171913 novelg_chr8_m_27171000_27171913
	} {
		set result [gene_name_check $chr $strand $begin $end checka checkisosa]
		if {$result ne $expected} {
			error "Wrong call for $name: $result instead of $expected"
		}
	}
	# test nested
	array set checkisosa {nested {{2000,3500, 2500,4000,}} enclosing {{1000,4500, 1500,5000,} }}
	unset -nocomplain checka
	array set checka {chr8,+ {{2000 4000 nested nested} {1000 5000 enclosing enclosing}} }
	foreach {name chr strand begin end expected} {
		nested chr8 + 2100 3000 nested
		nestedoutside chr8 + 1900 4100 nested
		enclosingstart chr8 + 1400 3000 enclosing
		enclosingend chr8 + 1900 4600 enclosing
	} {
		set result [gene_name_check $chr $strand $begin $end checka checkisosa]
		if {$result ne $expected} {
			error "Wrong call for $name: $result instead of $expected"
		}
	}
} {}

testsummarize
