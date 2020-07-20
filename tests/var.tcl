#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test var {var_gatk basic} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_gatk {*}$::dopts tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	cg unzip tmp/var-gatk-bwa.tsv.zst 
	exec diff -I {#GATKCommandLine.UnifiedGenotyper} tmp/var-gatk-bwa.tsv data/var-gatk-bwa.tsv
	string_change [cg covered tmp/sreg-gatk-bwa.tsv.zst] [list \n\n \n]
} {chromosome	bases
chr21	1047
chr22	156
total	1203}

test var {var distrreg gatk} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var {*}$::dopts -method gatk -distrreg 1 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	cg tsvdiff tmp/var-gatk-bwa.tsv.zst data/var-gatk-bwa.tsv
	string_change [cg covered tmp/sreg-gatk-bwa.tsv.zst] [list \n\n \n]
} {chromosome	bases
chr21	1047
chr22	156
total	1203}

test var {var distrreg gatk -d 3} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	exec cg var {*}$::dopts -d 3 -method gatk -distrreg 1 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	cg tsvdiff tmp/var-gatk-bwa.tsv.zst data/var-gatk-bwa.tsv
	string_change [cg covered tmp/sreg-gatk-bwa.tsv.zst] [list \n\n \n]
} {chromosome	bases
chr21	1047
chr22	156
total	1203}

test var {var distrreg gatk sequencedgenome (part)} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg select -q {$chromosome regexp "chr2."} [gzfile $::refseqdir/hg19/extra/reg_hg19_sequencedgenome.tsv] tmp/distrreg.tsv
	cg var {*}$::dopts -method gatk -distrreg tmp/distrreg.tsv tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	cg tsvdiff tmp/var-gatk-bwa.tsv.zst data/var-gatk-bwa.tsv
	string_change [cg covered tmp/sreg-gatk-bwa.tsv.zst] [list \n\n \n]
} {chromosome	bases
chr21	1047
chr22	156
total	1203}

test var {var distrreg gatk size} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var {*}$::dopts -method gatk -distrreg 30000000 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	cg tsvdiff tmp/var-gatk-bwa.tsv.zst data/var-gatk-bwa.tsv
	string_change [cg covered tmp/sreg-gatk-bwa.tsv.zst] [list \n\n \n]
} {chromosome	bases
chr21	1047
chr22	156
total	1203}

test var {var distrreg gatk result exists already} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	foreach file {
		tmp/sreg-cov3-bwa.tsv.zst
		tmp/varall-gatk-bwa.tsv.zst tmp/varall-gatk-bwa.tsv.zst.zsti tmp/varall-gatk-bwa.tsv.analysisinfo
		tmp/var-gatk-bwa.tsv.zst tmp/var-gatk-bwa.tsv.zst.zsti tmp/var-gatk-bwa.tsv.analysisinfo
		tmp/sreg-gatk-bwa.tsv.zst tmp/sreg-gatk-bwa.tsv.zst.zsti tmp/sreg-gatk-bwa.tsv.analysisinfo
		tmp/reg_cluster-gatk-bwa.tsv.zst
	} {
		after 100
		file_write $file {}
	}
	cg var {*}$::dopts -v 2 -distrreg chr -method gatk tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	catch {exec grep "targets already completed or running" tmp/logerror} temp
	if {$temp eq ""} {error "not skipped"}
	llength [glob -nocomplain tmp/*.old]
} 0

test var {var_sam basic} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_sam {*}$::dopts tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	cg tsvdiff tmp/var-sam-bwa.tsv.zst data/var-sam-bwa.tsv
	string_change [cg covered tmp/sreg-sam-bwa.tsv.zst] [list \n\n \n]
} {chromosome	bases
chr21	1047
chr22	156
total	1203}

test var {var_bcf basic} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_bcf {*}$::dopts -mincoverage 5 -minqual 12 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	cg tsvdiff tmp/var-bcf-bwa.tsv.zst data/var-bcf-bwa.tsv
	string_change [cg covered tmp/sreg-bcf-bwa.tsv.zst] [list \n\n \n]
} {chromosome	bases
chr21	1049
chr22	156
total	1205}

test var {var_bcf distrreg} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var {*}$::dopts -method bcf -distrreg 1 -mincoverage 5 -minqual 12 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	cg tsvdiff tmp/var-bcf-bwa.tsv.zst data/var-bcf-bwa.tsv
	string_change [cg covered tmp/sreg-bcf-bwa.tsv.zst] [list \n\n \n]
} {chromosome	bases
chr21	1049
chr22	156
total	1205}

test var {var distrreg sam} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var {*}$::dopts -method sam -distrreg 1 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	cg tsvdiff tmp/var-sam-bwa.tsv.zst data/var-sam-bwa.tsv
	string_change [cg covered tmp/sreg-sam-bwa.tsv.zst] [list \n\n \n]
} {chromosome	bases
chr21	1047
chr22	156
total	1203}

test var {var distrreg sam -d 3} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var {*}$::dopts -d 3 -method sam -distrreg 1 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	cg tsvdiff tmp/var-sam-bwa.tsv.zst data/var-sam-bwa.tsv
	string_change [cg covered tmp/sreg-sam-bwa.tsv.zst] [list \n\n \n]
} {chromosome	bases
chr21	1047
chr22	156
total	1203}

test var {var sam result exists already} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var {*}$::dopts -method sam tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas
	cg var {*}$::dopts -method sam tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	catch {exec grep "targets already completed or running" tmp/logerror} temp
	if {$temp eq ""} {error "not skipped"}
	llength [glob -nocomplain tmp/*.old]
} 0

test var {var_freebayes basic} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	exec cg var_freebayes {*}$::dopts tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	cg tsvdiff tmp/var-freebayes-bwa.tsv.zst data/var-freebayes-bwa.tsv
} {}

test var {var distrreg freebayes} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var {*}$::dopts -method freebayes -distrreg 1 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	cg tsvdiff tmp/var-freebayes-bwa.tsv.zst data/var-freebayes-bwa.tsv
} {}

test var {var_gatkh -ERC BP_RESOLUTION} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	# use low quality settings because test bwa.bam has low coverage
	cg var_gatkh {*}$::dopts -mincoverage 5 -mingenoqual 12 -ERC BP_RESOLUTION tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	cg tsvdiff tmp/var-gatkh-bwa.tsv.zst data/var-gatkh-bwa.tsv
	cg tsvdiff tmp/varall-gatkh-bwa.gvcf.gz data/varall-gatkh-bwa.gvcf.gz
	string_change [cg covered tmp/sreg-gatkh-bwa.tsv.zst] [list \n\n \n]
} {chromosome	bases
chr21	1019
chr22	142
total	1161}

test var {var_gatkh basic} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	# use low quality settings because test bwa.bam has low coverage
	cg var_gatkh {*}$::dopts -mincoverage 5 -mingenoqual 12 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	cg select -overwrite 1 -rf MIN_DP tmp/var-gatkh-bwa.tsv.zst tmp/test.tsv.zst
	cg tsvdiff tmp/test.tsv.zst data/var-gatkh-bwa.tsv
	cg tsvdiff tmp/varall-gatkh-bwa.gvcf.gz data/varall-gatkh-bwa.bgvcf
	string_change [cg covered tmp/sreg-gatkh-bwa.tsv.zst] [list \n\n \n]
} {chromosome	bases
chr21	1001
chr22	148
total	1149}

test var {var distrreg gatkh} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var {*}$::dopts -ERC BP_RESOLUTION -cleanup 0 -distrreg 1 -method gatkh -mincoverage 5 -mingenoqual 12 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	cg tsvdiff tmp/var-gatkh-bwa.tsv.zst data/var-gatkh-bwa.tsv
	cg tsvdiff tmp/varall-gatkh-bwa.gvcf.gz data/varall-gatkh-bwa.gvcf.gz
	string_change [cg covered tmp/sreg-gatkh-bwa.tsv.zst] [list \n\n \n]
} {chromosome	bases
chr21	1019
chr22	142
total	1161}

test var {var_strelka basic} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	cg var_strelka -stack 1 -v 2 -datatype exome -mincoverage 5 -mingenoqual 12 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log 2> tmp/logerror
	cg tsvdiff tmp/var-strelka-bwa.tsv.zst data/var-strelka-bwa.tsv
	cg tsvdiff tmp/varall-strelka-bwa.gvcf.gz data/varall-strelka-bwa.gvcf.gz
	string_change [cg covered tmp/sreg-strelka-bwa.tsv.zst] [list \n\n \n]
} {chromosome	bases
chr21	610
chr22	81
total	691}

test var {var distrreg strelka} {
	test_cleantmp
	file copy data/bwa.bam data/bwa.bam.bai tmp
	exec cg var {*}$::dopts -d 2 -datatype exome -mincoverage 5 -distrreg 1 -method strelka -mingenoqual 12 tmp/bwa.bam $::refseqdir/hg19/genome_hg19.ifas > tmp/log.txt 2> tmp/log.err
	cg tsvdiff tmp/var-strelka-bwa.tsv.zst data/var-strelka-bwa.tsv
	cg tsvdiff tmp/varall-strelka-bwa.gvcf.gz data/varall-strelka-bwa.gvcf.gz
} {}

test var {var_longshot basic} {
	cd $::smalltestdir
	file delete -force tmp/longshot
	file mkdir tmp/longshot
	foreach file [glob ori/longshot_example_data/pacbio_reads_30x.bam* ori/longshot_example_data/genome.fa* ori/longshot_example_data/ground_truth_variants.*] {
		mklink $file tmp/longshot/[file tail $file]
	}
	cg var_longshot {*}$::dopts tmp/longshot/pacbio_reads_30x.bam tmp/longshot/genome.fa
	cg vcf2tsv tmp/longshot/ground_truth_variants.vcf tmp/longshot/ground_truth_variants.tsv
	cg multicompar tmp/longshot/compar.tsv tmp/longshot/var-longshot-pacbio_reads_30x.tsv.zst tmp/longshot/ground_truth_variants.tsv
	cg tsvdiff -x *.log -x *.finished  -x *.zsti \
		tmp/longshot expected/longshot
	cg select -g {zyg-longshot-pacbio_reads_30x * zyg-ground_truth_variants *} tmp/longshot/compar.tsv
} {zyg-longshot-pacbio_reads_30x	zyg-ground_truth_variants	count
?	m	2
?	t	9
m	m	228
t	?	1
t	t	475}

test var {var_longshot distrreg} {
	cd $::smalltestdir
	file delete -force tmp/longshot_d
	file mkdir tmp/longshot_d
	foreach file [glob ori/longshot_example_data/pacbio_reads_30x.bam* ori/longshot_example_data/genome.fa* ori/longshot_example_data/ground_truth_variants.*] {
		mklink $file tmp/longshot_d/[file tail $file]
	}
	cg var -method longshot -distrreg 1 {*}$::dopts tmp/longshot_d/pacbio_reads_30x.bam tmp/longshot_d/genome.fa
	cg vcf2tsv tmp/longshot_d/ground_truth_variants.vcf tmp/longshot_d/ground_truth_variants.tsv
	cg multicompar tmp/longshot_d/compar.tsv tmp/longshot_d/var-longshot-pacbio_reads_30x.tsv.zst tmp/longshot_d/ground_truth_variants.tsv
	cg tsvdiff -brief 1 -x *.log -x *.finished -x *.zsti \
		tmp/longshot_d expected/longshot_d
	cg select -g {zyg-longshot-pacbio_reads_30x * zyg-ground_truth_variants *} \
		tmp/longshot_d/compar.tsv
} {zyg-longshot-pacbio_reads_30x	zyg-ground_truth_variants	count
?	m	2
?	t	9
m	m	228
t	?	1
t	t	475}

test var {var_longshot distrreg with contig1_ and empty contig and contig not in bam} {
	cd $::smalltestdir
	set files [glob -nocomplain tmp/longshot_d2/*]
	if {[llength $files]} {
		file delete -force {*}$files
	}
	file mkdir tmp/longshot_d2
	set c [file_read ori/longshot_example_data/genome.fa]
	set c [string_change $c {>contig1 >contig1_random1 >contig2 >contig1_random2}]
	append c >contig4\nAAAAAAAAAA\n>contig5\nAAAAAAAAAA
	file_write tmp/longshot_d2/genome.fa $c\n
	exec samtools faidx tmp/longshot_d2/genome.fa
	mklink ori/longshot_example_data/renamed_pacbio_reads_30x.bam tmp/longshot_d2/pacbio_reads_30x.bam
	exec samtools index tmp/longshot_d2/pacbio_reads_30x.bam
	set c [file_read ori/longshot_example_data/ground_truth_variants.tsv]
	set c [string_change $c {contig1 contig1_random1 contig2 contig1_random2}]
	file_write tmp/longshot_d2/ground_truth_variants.tsv $c\n
	#
	cg var -method longshot -distrreg 1 {*}$::dopts tmp/longshot_d2/pacbio_reads_30x.bam tmp/longshot_d2/genome.fa
	cg multicompar tmp/longshot_d2/compar.tsv tmp/longshot_d2/var-longshot-pacbio_reads_30x.tsv.zst tmp/longshot_d2/ground_truth_variants.tsv
	cg tsvdiff -brief 1 -x *.log -x *.bai -x *.finished -x *.zsti \
		tmp/longshot_d2 expected/longshot_d2
	cg select -g {zyg-longshot-pacbio_reads_30x * zyg-ground_truth_variants *} tmp/longshot_d2/compar.tsv
} {zyg-longshot-pacbio_reads_30x	zyg-ground_truth_variants	count
?	m	2
?	t	9
m	m	228
t	?	1
t	t	475}

test var {var_longshot distrreg with -hap_bam 1 option} {
	cd $::smalltestdir
	file delete -force tmp/longshot_hapbam
	file mkdir tmp/longshot_hapbam
	foreach file [glob ori/longshot_example_data/pacbio_reads_30x.bam* ori/longshot_example_data/genome.fa* ori/longshot_example_data/ground_truth_variants.*] {
		mklink $file tmp/longshot_hapbam/[file tail $file]
	}
	cg var {*}$::dopts -method longshot -distrreg 1 -hap_bam 1 tmp/longshot_hapbam/pacbio_reads_30x.bam tmp/longshot_hapbam/genome.fa
	cg vcf2tsv tmp/longshot_hapbam/ground_truth_variants.vcf tmp/longshot_hapbam/ground_truth_variants.tsv
	cg multicompar tmp/longshot_hapbam/compar.tsv tmp/longshot_hapbam/var-longshot-pacbio_reads_30x.tsv.zst tmp/longshot_hapbam/ground_truth_variants.tsv
	cg tsvdiff -brief 1 -x *.log -x *.finished -x *.zsti \
		tmp/longshot_hapbam expected/longshot_hapbam
	set temp [string trim [exec cg sam2tsv -fields HP tmp/longshot_hapbam/map-hlongshot-pacbio_reads_30x.bam | cg select -q {$HP == 1} | cg regjoin | cg  covered]]
	if {$temp ne [string trim [deindent {
		chromosome	bases
		contig1	185941
		contig2	179545
		contig3	176005
		
		total	541491
	}]]} {
		error "error in haplotyped bams"
	}
	set expected [string trim [deindent {
		zyg-longshot-pacbio_reads_30x	zyg-ground_truth_variants	count
		?	m	2
		?	t	9
		m	m	228
		t	?	1
		t	t	475
	}]]
	set result [exec cg select -g {zyg-longshot-pacbio_reads_30x * zyg-ground_truth_variants *} tmp/longshot_hapbam/compar.tsv]
	if {$result ne $expected} {
		error "wrong result"
	}
} {}

testsummarize

