#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"
# small (limited) process_sample and process_project tests using only data directly distributed with genomecomb

source tools.tcl

proc splitfastqs {fastq1 pre1 post1 {perfile 20}} {
	set num 1
	catch {close $o1 ; close $f1}
	set f1 [gzopen $fastq1]
	set o1 [open $pre1${num}$post1 w]
	set lines 0
	while {[gets $f1 line1] != -1} {
		puts $o1 $line1
		incr lines
		if {$lines == $perfile} {
			close $o1
			if {[eof $f1]} break
			incr num
			set o1 [open $pre1${num}$post1 w]
			set lines 0
		}
	}
	close $o1; close $f1
	if {$lines == 0} {file delete $pre1${num}$post1}
}

test process_sample {bwa distrreg} {
	test_cleantmp
	file mkdir tmp/NA19240m/fastq
	file copy -force data/seq_R1.fq.gz data/seq_R2.fq.gz tmp/NA19240m/fastq
	exec cg process_sample {*}$::dopts -threads 1 -aligners bwa -distrreg chr \
		-dbdir $::refseqdir/hg19/genome_hg19.ifas tmp/NA19240m > tmp/NA19240m.startuplog 2> tmp/NA19240m.startuperror
	# chr21:42730799-42762826
	set genomecombversion [cg version]
	file_write tmp/expected_varall-gatk-rdsbwa-NA19240m.tsv.analysisinfo [subst [deindent {
		sample	clipping	clipping_version	clipping_cg_version	aligner	aligner_version	reference	aligner_paired	aligner_sort	aligner_sort_version	sammerge	sammerge_version	sammerge_sort	sammerge_mergesort	bamclean	bamclean_version	removeduplicates	removeduplicates_version	realign	realign_version	analysis	varcaller	varcaller_version	varcaller_cg_version	varcaller_region
		gatk-rdsbwa-NA19240m	fastq-mcf	1.1.2-537 adapted	$genomecombversion	bwa	0.7.15-r1140	hg19	1	gnusort	8.31	genomecomb	$genomecombversion	coordinate	1	genomecomb	$genomecombversion	samtools	1.15 (using htslib 1.15)	gatk	3.8-1-0-gf15c1c3ef	gatk-rdsbwa-NA19240m	gatk	3.8-1-0-gf15c1c3ef	$genomecombversion	sreg-cov5-rdsbwa-NA19240m.tsv.zst
	}]]\n
	cg tsvdiff -q 1 -x fastq -x *.bai -x *.crai -x *.zsti \
		-x projectinfo.tsv -x *.analysisinfo -x *.stats.zst \
		-x info_analysis.tsv -x *.submitting -x *.finished -x fastqc* -x *.index -x log_jobs \
		tmp/NA19240m data/NA19240m
	cg tsvdiff tmp/NA19240m/varall-gatk-rdsbwa-NA19240m.tsv.analysisinfo tmp/expected_varall-gatk-rdsbwa-NA19240m.tsv.analysisinfo
} {}

test process_sample {bwa distrreg ubam source} {
	test_cleantmp
	file mkdir tmp/NA19240m/ubam
	exec samtools import data/seq_R1.fq.gz data/seq_R2.fq.gz > tmp/NA19240m/ubam/seq.bam
	exec cg process_sample {*}$::dopts -clip 0 -threads 1 -aligners bwa -distrreg chr \
		-dbdir $::refseqdir/hg19/genome_hg19.ifas tmp/NA19240m > tmp/NA19240m.startuplog 2> tmp/NA19240m.startuperror
	# chr21:42730799-42762826
	set genomecombversion [cg version]
	file_write tmp/expected_varall-gatk-rdsbwa-NA19240m.tsv.analysisinfo [subst [deindent {
		sample	clipping	clipping_version	clipping_cg_version	aligner	aligner_version	reference	aligner_paired	aligner_sort	aligner_sort_version	sammerge	sammerge_version	sammerge_sort	sammerge_mergesort	bamclean	bamclean_version	removeduplicates	removeduplicates_version	realign	realign_version	analysis	varcaller	varcaller_version	varcaller_cg_version	varcaller_region
		gatk-rdsbwa-NA19240m	fastq-mcf	1.1.2-537 adapted	$genomecombversion	bwa	0.7.15-r1140	hg19	1	gnusort	8.31	genomecomb	$genomecombversion	coordinate	1	genomecomb	$genomecombversion	samtools	1.15 (using htslib 1.15)	gatk	3.8-1-0-gf15c1c3ef	gatk-rdsbwa-NA19240m	gatk	3.8-1-0-gf15c1c3ef	$genomecombversion	sreg-cov5-rdsbwa-NA19240m.tsv.zst
	}]]\n
	cg tsvdiff -q 1 -x fastq -x *.bai -x *.crai -x *.zsti \
		-x projectinfo.tsv -x *.analysisinfo -x *.stats.zst \
		-x info_analysis.tsv* -x *.submitting -x *.finished -x fastqc* -x *.index -x log_jobs \
		-x fastq_stats* -x fastx* -x report_fastq_* -x ubam \
		tmp/NA19240m data/NA19240m
} {}

test process_sample {bwa distrreg cram} {
	test_cleantmp
	file mkdir tmp/NA19240m/fastq
	file copy -force data/seq_R1.fq.gz data/seq_R2.fq.gz tmp/NA19240m/fastq
	exec cg process_sample -threads 1 -aligners bwa -aliformat cram -distrreg chr \
		-dbdir $::refseqdir/hg19/genome_hg19.ifas tmp/NA19240m >& tmp/NA19240m.startuplog
	# chr21:42730799-42762826
	set genomecombversion [cg version]
	file_write tmp/expected_varall-gatk-rdsbwa-NA19240m.tsv.analysisinfo [deindent [subst {
		sample	clipping	clipping_version	clipping_cg_version	aligner	aligner_version	reference	aligner_paired	aligner_sort	aligner_sort_version	sammerge	sammerge_version	sammerge_sort	sammerge_mergesort	bamclean	bamclean_version	removeduplicates	removeduplicates_version	realign	realign_version	analysis	varcaller	varcaller_version	varcaller_cg_version	varcaller_region
		gatk-rdsbwa-NA19240m	fastq-mcf	1.1.2-537 adapted	$genomecombversion	bwa	0.7.15-r1140	hg19	1	gnusort	8.31	genomecomb	$genomecombversion	coordinate	1	genomecomb	$genomecombversion	samtools	1.15 (using htslib 1.15)	gatk	3.8-1-0-gf15c1c3ef	gatk-rdsbwa-NA19240m	gatk	3.8-1-0-gf15c1c3ef	$genomecombversion	sreg-cov5-rdsbwa-NA19240m.tsv.zst
	}]]\n
	set result {}
	lappend result [tsvdiff -q 1 -x fastq -x *.bai -x *.crai -x *.zsti \
		-x *.cram -x *.sam -x *.stats.zst \
		-x projectinfo.tsv -x *.analysisinfo \
		-x info_analysis.tsv -x *.submitting -x *.finished -x fastqc* -x *.index -x log_jobs \
		tmp/NA19240m data/NA19240m]
	# do cram diff separate as there may be differences in the "other" fields
	cg sam2tsv tmp/NA19240m/map-rdsbwa-NA19240m.cram | cg select -rf other > tmp/temp.tsv
	cg sam2tsv data/NA19240m/map-rdsbwa-NA19240m.sam | cg select -rf other > tmp/expected.tsv
	lappend result [tsvdiff tmp/temp.tsv tmp/expected.tsv]
	lappend result [tsvdiff tmp/NA19240m/varall-gatk-rdsbwa-NA19240m.tsv.analysisinfo tmp/expected_varall-gatk-rdsbwa-NA19240m.tsv.analysisinfo]
	join [list_remove $result {}] \n
} {}

test process_sample {bwa distrreg -removeduplicates 0 cram} {
	test_cleantmp
	file mkdir tmp/NA19240m/fastq
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp/NA19240m/fastq
	cg process_sample -aligners bwa -varcallers {} -removeduplicates 0 -aliformat cram -distrreg chr -dbdir $::refseqdir/hg19/genome_hg19.ifas tmp/NA19240m
	# chr21:42730799-42762826
	dbdir $::refseqdir/hg19
	exec samtools sort --no-PG -O sam tmp/NA19240m/map-rsbwa-NA19240m.cram > tmp/ali.sam
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} tmp/ali.sam | cg select -s {chromosome begin end qname} > tmp/ali.sam.tsv
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} data/bwa.sam | cg select -s {chromosome begin end qname} > tmp/bwa.sam.tsv
	set result {}
	lappend result [tsvdiff tmp/ali.sam.tsv tmp/bwa.sam.tsv]
	join [list_remove $result {}] \n
} {}

test process_sample {bwa distrreg -removeduplicates picard bam} {
	test_cleantmp
	file mkdir tmp/NA19240m/fastq
	file copy data/seq_R1.fq.gz data/seq_R2.fq.gz tmp/NA19240m/fastq
	cg process_sample -stack 1 -aligners bwa -varcallers {} -removeduplicates picard -distrreg chr -dbdir $::refseqdir/hg19/genome_hg19.ifas tmp/NA19240m
	# chr21:42730799-42762826
	exec samtools sort --no-PG -O sam tmp/NA19240m/map-rdsbwa-NA19240m.bam > tmp/ali.sam
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} tmp/ali.sam | cg select -rf other -s {chromosome begin end qname} > tmp/ali.sam.tsv
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} data/NA19240m/map-rdsbwa-NA19240m.sam | cg select -rf other -s {chromosome begin end qname} > tmp/bwa.sam.tsv
	set result {}
	lappend result [tsvdiff tmp/ali.sam.tsv tmp/bwa.sam.tsv]
	join [list_remove $result {}] \n
} {}

test process_sample {map bwa distrreg mutiple fastq} {
	test_cleantmp
	file mkdir tmp/NA19240m/fastq
	# make split up fastq
	splitfastqs data/seq_R1.fq.gz tmp/NA19240m/fastq/seq _R1.fq 20
	splitfastqs data/seq_R2.fq.gz tmp/NA19240m/fastq/seq _R2.fq 20
	file delete tmp/NA19240m/map-rdsbwa-NA19240m.bam
	cg process_sample -clip 0 -aligners bwa -varcallers {} -distrreg chr -varcallers {} \
		-dbdir $::refseqdir/hg19/genome_hg19.ifas tmp/NA19240m >@ stdout 2>@ stderr
	# chr21:42730799-42762826
	exec samtools sort --no-PG -O sam tmp/NA19240m/map-rdsbwa-NA19240m.bam > tmp/ali.sam
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} tmp/ali.sam | cg select -rf {duplicate other properpair mapquality XS MQ} -s {chromosome begin end qname} > tmp/ali.sam.tsv
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} data/bwa.sam | cg select -rf {duplicate other properpair mapquality XS MQ} -s {chromosome begin end qname} > tmp/bwa.sam.tsv
	set result {}
	lappend result [tsvdiff tmp/ali.sam.tsv tmp/bwa.sam.tsv]
	join [list_remove $result {}] \n
} {}

test process_sample {map bwa distrreg mutiple fastq -maxfastqdistr 2} {
	test_cleantmp
	file mkdir tmp/NA19240m/fastq
	# make split up fastq
	splitfastqs data/seq_R1.fq.gz tmp/NA19240m/fastq/seq _R1.fq 20
	splitfastqs data/seq_R2.fq.gz tmp/NA19240m/fastq/seq _R2.fq 20
	file delete tmp/NA19240m/map-rdsbwa-NA19240m.bam
	cg process_sample -stack 1 -clip 0 -maxfastqdistr 2 -aligners bwa -varcallers {} -distrreg chr -varcallers {} -dbdir $::refseqdir/hg19/genome_hg19.ifas tmp/NA19240m >@ stdout 2>@ stderr
	# chr21:42730799-42762826
	exec samtools sort --no-PG -O sam tmp/NA19240m/map-rdsbwa-NA19240m.bam > tmp/ali.sam
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} tmp/ali.sam | cg select -rf {duplicate other properpair mapquality XS MQ} -s {chromosome begin end qname} > tmp/ali.sam.tsv
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} data/bwa.sam | cg select -rf {duplicate other properpair mapquality XS MQ} -s {chromosome begin end qname} > tmp/bwa.sam.tsv
	set result {}
	lappend result [tsvdiff tmp/ali.sam.tsv tmp/bwa.sam.tsv]
	join [list_remove $result {}] \n
} {}

test process_project {process_project bwa distrreg multiple fastq -maxfastqdistr 2} {
	test_cleantmp
	file mkdir tmp/samples/NA19240m/fastq
	# make split up fastq
	splitfastqs data/seq_R1.fq.gz tmp/samples/NA19240m/fastq/seq _R1.fq 20
	splitfastqs data/seq_R2.fq.gz tmp/samples/NA19240m/fastq/seq _R2.fq 20
	cg process_project -stack 1 -clip 0 -maxfastqdistr 2 -aligners bwa -varcallers {} -distrreg chr -varcallers {} -dbdir $::refseqdir/hg19/genome_hg19.ifas tmp >@ stdout 2>@ stderr
	# chr21:42730799-42762826
	exec samtools sort --no-PG -O sam tmp/samples/NA19240m/map-rdsbwa-NA19240m.bam > tmp/ali.sam
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} tmp/ali.sam | cg select -rf {duplicate other properpair mapquality XS MQ} -s {chromosome begin end qname} > tmp/ali.sam.tsv
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} data/bwa.sam | cg select -rf {duplicate other properpair mapquality XS MQ} -s {chromosome begin end qname} > tmp/bwa.sam.tsv
	set result {}
	lappend result [tsvdiff tmp/ali.sam.tsv tmp/bwa.sam.tsv]
	join [list_remove $result {}] \n
} {}

test process_project {process_project include msamples directory (analyse, but not include in compar)} {
	test_cleantmp
	cg project_addsample -transfer rel tmp NA19240m data/seq_R1.fq.gz data/seq_R2.fq.gz
	cg project_addsample -transfer rel tmp msample data/seq_R1.fq.gz data/seq_R2.fq.gz
	mkdir tmp/msamples
	file rename tmp/samples/msample tmp/msamples/msample
	cg process_project -stack 1 \
		-clip 0 -maxfastqdistr 2 -aligners bwa -varcallers bcf \
		-distrreg chr -dbdir $::refseqdir/hg19 \
		tmp >& tmp/startup.log
	# chr21:42730799-42762826
	set result {}
	exec samtools sort --no-PG -O sam tmp/samples/NA19240m/map-rdsbwa-NA19240m.bam > tmp/ali.sam
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} tmp/ali.sam | cg select -rf {duplicate other properpair mapquality XS MQ} -s {chromosome begin end qname} > tmp/ali.sam.tsv
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} data/bwa.sam | cg select -rf {duplicate other properpair mapquality XS MQ} -s {chromosome begin end qname} > tmp/bwa.sam.tsv
	lappend result [tsvdiff tmp/ali.sam.tsv tmp/bwa.sam.tsv]
	exec samtools sort --no-PG -O sam tmp/msamples/msample/map-rdsbwa-msample.bam > tmp/ali2.sam
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} tmp/ali2.sam | cg select -rf {duplicate other properpair mapquality XS MQ NM RG} -s {chromosome begin end qname} > tmp/ali2.sam.tsv
	exec cg sam2tsv -fields {AS XS MQ MC ms MD RG NM XA} data/bwa.sam | cg select -rf {duplicate other properpair mapquality XS MQ NM RG} -s {chromosome begin end qname} > tmp/bwa2.sam.tsv
	lappend result [tsvdiff tmp/ali2.sam.tsv tmp/bwa2.sam.tsv]
	lappend result [cg select -n tmp/compar/annot_compar-tmp.tsv.zst]
	join [list_remove $result {}] \n
} NA19240m

testsummarize

