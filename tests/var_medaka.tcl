#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

proc make_smallgiabonttest {dir} {
	file delete -force $dir
	file mkdir $dir
	set regions {chr1:2547867-2568902 chr6:32152328-32167543 chr10:975157-1000215}
	set oridir ori/nanopore-human-pangenomics_regions/HG002
	exec samtools view -b \
		$oridir/map-sminimap2-regions_HG002_hg38.bam \
		{*}$regions \
		> $dir/map-sminimap2-pHG002_hg38.bam
	exec samtools index $dir/map-sminimap2-pHG002_hg38.bam
	file_write $dir/regions.tsv chromosome\tbegin\tend\n[string_change [join $regions \n] [list : \t - \t]]\n
	cg regcommon $dir/regions.tsv $oridir/regions_HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed.tsv.zst > $dir/sreg-truth_HG002_hg38.tsv
	cg regselect $oridir/regions_HG002_GRCh38_1_22_v4.2.1_benchmark.tsv.zst \
		$dir/sreg-truth_HG002_hg38.tsv > $dir/var-truth_HG002_hg38.tsv
}

test var {var_medaka basic giab ont} {
	cd $::smalltestdir
	set workdir tmp/medaka_sgiab
	make_smallgiabonttest $workdir
	cg var_medaka {*}$::dopts \
		$workdir/map-sminimap2-pHG002_hg38.bam $::refseqdir/hg38 >& tmp/medaka_sgiab.log
	file delete $workdir/compar.tsv
	cg multicompar -reannot 1 $workdir/compar.tsv \
		$workdir/var-medaka-sminimap2-pHG002_hg38.tsv.zst \
		$workdir/var-truth_HG002_hg38.tsv
	cg benchmarkvars -refcurve_cutoffs {{} 10 20 30 40 50 60} $workdir/compar.tsv truth_HG002_hg38 $workdir/benchmark.tsv
	cg tsvdiff -q 1 -x *.log -x *.finished  -x *.zsti \
		-x compar.tsv.reannot \
		-ignorefields {varcaller_cg_version} \
		$workdir expected/[file tail $workdir]
	list [cg select -g chromosome $workdir/compar.tsv] [cg select -g {zyg-medaka-sminimap2-pHG002_hg38 * zyg-truth_HG002_hg38 *} $workdir/compar.tsv]
} {{chromosome	count
1	10071
6	8518
10	4349} {zyg-medaka-sminimap2-pHG002_hg38	zyg-truth_HG002_hg38	count
c	u	908
m	m	24
m	r	1
m	u	686
r	c	1
r	t	1
t	c	1
t	r	12
t	t	59
t	u	21245}}

test var {var_medaka basic} {
	cd $::smalltestdir
	file delete -force tmp/medaka
	file mkdir tmp/medaka
	mklink ori/pepperdeepvariant_example_data/HG002_ONT_50x_2_GRCh38.chr20.quickstart.bam tmp/medaka/test.bam
	mklink ori/pepperdeepvariant_example_data/HG002_ONT_50x_2_GRCh38.chr20.quickstart.bam.bai tmp/medaka/test.bam.bai
	cg vcf2tsv ori/pepperdeepvariant_example_data/HG002_GRCh38_1_22_v4.2.1_benchmark.quickstart.vcf.gz tmp/medaka/var-truth.tsv
	cg bed2tsv ori/pepperdeepvariant_example_data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.quickstart.bed tmp/medaka/sreg-truth.tsv
	#
	cg var_medaka {*}$::dopts tmp/medaka/test.bam $::refseqdir/hg38 \
		>& tmp/medaka.log
	file delete tmp/medaka/compar.tsv
	cg multicompar -reannot 1 tmp/medaka/compar.tsv tmp/medaka/var-medaka-test.tsv.zst tmp/medaka/var-truth.tsv
	cg benchmarkvars -refcurve_cutoffs {{} 10 20 30 40 50 60} tmp/medaka/compar.tsv truth tmp/medaka/benchmark.tsv
	set result {}
	lappend result [tsvdiff -q 1 \
		-x *.log -x *.finished  -x *.zsti -x *.submitting -x *.tsv.reannot -x *.tbi \
		-ignorefields {varcaller_cg_version sammerge_version} \
		tmp/medaka expected/medaka]
	lappend result [cg select -g chromosome tmp/medaka/compar.tsv]
	lappend result [cg select -g {zyg-medaka-test * zyg-truth *} tmp/medaka/compar.tsv]
	join [list_remove $result {}] \n
} {chromosome	count
20	10014
zyg-medaka-test	zyg-truth	count
c	u	220
m	m	36
m	r	1
m	u	212
t	r	1
t	t	12
t	u	9532}

test var {var_medaka distrreg} {
	cd $::smalltestdir
	file delete -force tmp/medaka_distrreg
	file mkdir tmp/medaka_distrreg
	mklink ori/pepperdeepvariant_example_data/HG002_ONT_50x_2_GRCh38.chr20.quickstart.bam tmp/medaka_distrreg/test.bam
	mklink ori/pepperdeepvariant_example_data/HG002_ONT_50x_2_GRCh38.chr20.quickstart.bam.bai tmp/medaka_distrreg/test.bam.bai
	cg vcf2tsv ori/pepperdeepvariant_example_data/HG002_GRCh38_1_22_v4.2.1_benchmark.quickstart.vcf.gz tmp/medaka_distrreg/var-truth.tsv
	cg bed2tsv ori/pepperdeepvariant_example_data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.quickstart.bed tmp/medaka_distrreg/sreg-truth.tsv
	#
	cg var -method medaka -distrreg 1 {*}$::dopts tmp/medaka_distrreg/test.bam $::refseqdir/hg38 \
		>& tmp/medaka_distrreg.log
	file delete tmp/medaka_distrreg/compar.tsv
	cg multicompar -reannot 1 tmp/medaka_distrreg/compar.tsv tmp/medaka_distrreg/var-medaka-test.tsv.zst tmp/medaka/var-truth.tsv
	cg benchmarkvars -refcurve_cutoffs {{} 10 20 30 40 50 60} tmp/medaka_distrreg/compar.tsv truth tmp/medaka_distrreg/benchmark.tsv
	set result {}
	lappend result [tsvdiff -q 1 \
		-x *.old -x *.log -x *.finished  -x *.zsti -x *.submitting -x *.tsv.reannot \
		-x *.tbi \
		-ignorefields {varcaller_cg_version sammerge_version} \
		tmp/medaka_distrreg expected/medaka_distrreg]
	lappend result [cg select -g chromosome tmp/medaka_distrreg/compar.tsv]
	lappend result [cg select -g {zyg-medaka-test * zyg-truth *} tmp/medaka_distrreg/compar.tsv]
	join [list_remove $result {}] \n
} {chromosome	count
20	10014
zyg-medaka-test	zyg-truth	count
c	u	220
m	m	36
m	r	1
m	u	212
t	r	1
t	t	12
t	u	9532}

test var {var -method medaka -regionfile} {
	# not all methods support -regionfile (e.g. medaka)
	# in this case the -regionfile parameter is ignored, and the entire genome is run
	cd $::smalltestdir
	file delete -force tmp/medaka_reg
	file mkdir tmp/medaka_reg
	mklink ori/pepperdeepvariant_example_data/HG002_ONT_50x_2_GRCh38.chr20.quickstart.bam tmp/medaka_reg/test.bam
	mklink ori/pepperdeepvariant_example_data/HG002_ONT_50x_2_GRCh38.chr20.quickstart.bam.bai tmp/medaka_reg/test.bam.bai
	cg vcf2tsv ori/pepperdeepvariant_example_data/HG002_GRCh38_1_22_v4.2.1_benchmark.quickstart.vcf.gz tmp/medaka_reg/var-truth.tsv
	cg bed2tsv ori/pepperdeepvariant_example_data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.quickstart.bed tmp/medaka_reg/sreg-truth.tsv
	file_write tmp/medaka_reg/targets.tsv [string trim [deindent {
		chromosome	begin	end
		chr20	831400	831600
		chr20	840000	841000
		chr20	1005000	1009000
	}]]\n
	#
	cg var -stack 1 -v 2 -method medaka \
		-regionfile tmp/medaka_reg/targets.tsv \
		{*}$::dopts \
		tmp/medaka_reg/test.bam $::refseqdir/hg38 \
		>& tmp/medaka_reg.log
	file delete tmp/medaka_reg/compar.tsv
	cg multicompar -reannot 1 \
		tmp/medaka_reg/compar.tsv tmp/medaka_reg/var-medaka-test.tsv.zst tmp/medaka_reg/var-truth.tsv
	set result {}
	lappend result [tsvdiff -q 1 \
		-x *.log -x *.finished  -x *.zsti -x *.submitting -x *.tsv.reannot -x *.tbi \
		-ignorefields {varcaller_cg_version sammerge_version} \
		tmp/medaka_reg expected/medaka_reg]
	lappend result [cg select -g chromosome tmp/medaka_reg/compar.tsv]
	lappend result [cg select -g {zyg-medaka-test * zyg-truth *} tmp/medaka_reg/compar.tsv]
	join [list_remove $result {}] \n
} {chromosome	count
20	10014
zyg-medaka-test	zyg-truth	count
c	u	220
m	m	36
m	r	1
m	u	212
t	r	1
t	t	12
t	u	9532}

test var {var -method medaka -regionfile -distreg} {
	cd $::smalltestdir
	file delete -force tmp/medaka_reg_distrreg
	file mkdir tmp/medaka_reg_distrreg
	mklink ori/pepperdeepvariant_example_data/HG002_ONT_50x_2_GRCh38.chr20.quickstart.bam tmp/medaka_reg_distrreg/test.bam
	mklink ori/pepperdeepvariant_example_data/HG002_ONT_50x_2_GRCh38.chr20.quickstart.bam.bai tmp/medaka_reg_distrreg/test.bam.bai
	cg vcf2tsv ori/pepperdeepvariant_example_data/HG002_GRCh38_1_22_v4.2.1_benchmark.quickstart.vcf.gz tmp/medaka_reg_distrreg/var-truth.tsv
	cg bed2tsv ori/pepperdeepvariant_example_data/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.quickstart.bed tmp/medaka_reg_distrreg/sreg-truth.tsv
	file_write tmp/medaka_reg_distrreg/targets.tsv [string trim [deindent {
		chromosome	begin	end
		chr20	831400	831600
		chr20	840000	841000
		chr20	1005000	1009000
	}]]\n
	#
	cg var -method medaka \
		-regionfile tmp/medaka_reg_distrreg/targets.tsv \
		-distrreg regionfile \
		{*}$::dopts \
		tmp/medaka_reg_distrreg/test.bam \
		$::refseqdir/hg38 \
		>& tmp/medaka_reg_distrreg.log
	file delete tmp/medaka_reg_distrreg/compar.tsv
	cg multicompar -reannot 1 \
		tmp/medaka_reg_distrreg/compar.tsv tmp/medaka_reg_distrreg/var-medaka-test.tsv.zst tmp/medaka_reg_distrreg/var-truth.tsv
	set result {}
	lappend result [tsvdiff -q 1 \
		-x *.log -x *.finished  -x *.zsti -x *.submitting -x *.tbi -x *.tsv.reannot \
		-ignorefields {varcaller_cg_version sammerge_version} \
		tmp/medaka_reg_distrreg expected/medaka_reg_distrreg]
	lappend result [cg select -g chromosome tmp/medaka_reg_distrreg/compar.tsv]
	lappend result [cg select -g {zyg-medaka-test * zyg-truth *} tmp/medaka_reg_distrreg/compar.tsv]
	join [list_remove $result {}] \n
} {chromosome	count
20	71
zyg-medaka-test	zyg-truth	count
m	m	6
m	u	23
u	m	30
u	t	12}

testsummarize
