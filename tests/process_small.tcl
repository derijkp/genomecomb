#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl
set keepdir [pwd]
set dopts [get argv ""]

# tests
# =====

test process_small {mastr mastr_mx2} {
	cd $::bigtestdir	
	file delete -force tmp/mastr_mx2
#	file copy ori/mastr_mx2 tmp/
	file mkdir tmp/mastr_mx2
	foreach sample [glob ori/mastr_mx2.start/samples/*] {
		set sdest tmp/mastr_mx2/[string range [file tail $sample] 0 end-3]/fastq
		file mkdir $sdest
		foreach file [glob -nocomplain $sample/ori/*] {
			file copy $file $sdest
		}
	}
	file delete -force tmp/wgs2.mastr
	file mkdir tmp/wgs2.mastr
	file copy ori/wgs2.mastr/amplicons-wgs2.tsv tmp/wgs2.mastr
	# file copy ori/mastr_mx2/demultiplex_stats.tsv tmp/mastr_mx2
	# if you want to see output while running
	 cg process_mastr --stack 1 --verbose 2 {*}$::dopts -split 1 tmp/wgs2.mastr tmp/mastr_mx2 refseqtest/hg19 2>@ stderr >@ stdout
	# no output while running
	# cg process_mastr --stack 1 --verbose 2 -split 1 tmp/wgs2.mastr tmp/mastr_mx2 refseqtest/hg19
	# check vs expected
	foreach sample {blanco2_8485 ceph1333_02_34_7220 ceph1347_02_34_7149 ceph1347_02_34_8446} {
		checkdiff -y --suppress-common-lines tmp/mastr_mx2/$sample/crsbwa-$sample.hsmetrics expected/mastr_mx2/$sample/crsbwa-$sample.hsmetrics | grep -v -E "Started on|net.sf.picard.analysis.directed.CalculateHsMetrics BAIT_INT"
	}
	checkdiff -qr -x *log_jobs -x *hsmetrics -x *.bam -x *.bai -x colinfo -x mastr_mx2.html tmp/mastr_mx2 expected/mastr_mx2
	checkdiff -y --suppress-common-lines tmp/mastr_mx2/mastr_mx2.html expected/mastr_mx2/mastr_mx2.html | grep -v -E {HistogramID|htmlwidget-|^<!|^<h2>20}
	# could have used this, but previous is faster
	# cg tsvdiff -q 1 -x log_jobs -x hsmetrics -x colinfo -x mastr_mx2.html tmp/mastr_mx2 expected/mastr_mx2 2> temp
} {}

test process_small {process_sample exome yri mx2} {
	cd $::bigtestdir
	file delete -force tmp/one_exome_yri_mx2
	file mkdir tmp/one_exome_yri_mx2/samples
	foreach sample {
		 NA19240mx2
	} {
		file mkdir tmp/one_exome_yri_mx2/samples/$sample/fastq
		file copy {*}[glob ori/exomes_yri_mx2.start/samples/$sample/ori/*.fq.gz] tmp/one_exome_yri_mx2/samples/$sample/fastq
	}
	cg process_sample --stack 1 --verbose 2 {*}$::dopts -split 1 -dbdir refseqtest/hg19 tmp/one_exome_yri_mx2/samples/NA19240mx2 2>@ stderr >@ stdout
	# cg process_sample --stack 1 --verbose 2 -d status -split 1 -dbdir refseqtest/hg19 tmp/one_exome_yri_mx2/samples/NA19240mx2 | less -S
	# check vs expected
	checkdiff -y --suppress-common-lines tmp/one_exome_yri_mx2/samples/NA19240mx2/map-dsbwa-NA19240mx2.bam.dupmetrics expected/one_exome_yri_mx2/samples/NA19240mx2/map-dsbwa-NA19240mx2.bam.dupmetrics | grep -v "Started on" | grep -v "net.sf.picard.sam.MarkDuplicates INPUT"
	checkdiff -qr -x *log_jobs -x *.bam -x *.bai -x colinfo -x *_fastqc -x *bam.dupmetrics tmp/one_exome_yri_mx2/samples/NA19240mx2 expected/one_exome_yri_mx2/samples/NA19240mx2
} {}

test process_small {process_illumina exomes yri mx2} {
	cd $::bigtestdir
	file delete -force tmp/exomes_yri_mx2
	file mkdir tmp/exomes_yri_mx2/samples
	foreach sample {
		NA19238mx2  NA19239mx2  NA19240mx2
	} {
		file mkdir tmp/exomes_yri_mx2/samples/$sample/fastq
		file copy {*}[glob ori/exomes_yri_mx2.start/samples/$sample/ori/*.fq.gz] tmp/exomes_yri_mx2/samples/$sample/fastq
	}
	# cg process_illumina --stack 1 --verbose 2 -d 2 -split 1 -dbdir refseqtest/hg19 tests/yri_exome
	cg process_illumina --stack 1 --verbose 2 {*}$::dopts -split 1 -dbdir refseqtest/hg19 tmp/exomes_yri_mx2 2>@ stderr >@ stdout
	# check vs expected
	checkdiff -y --suppress-common-lines tmp/exomes_yri_mx2/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics expected/exomes_yri_mx2/samples/NA19238mx2/map-dsbwa-NA19238mx2.bam.dupmetrics | grep -v "Started on"
	checkdiff -qr -x *log_jobs -x *.bam -x *.bai -x colinfo -x *_fastqc -x *bam.dupmetrics tmp/exomes_yri_mx2 expected/exomes_yri_mx2
	# could have used this, but previous is faster
	# cg tsvdiff -q 1 -x log_jobs -x colinfo -x _fastqc -x bam.dupmetrics tmp/exomes_yri_mx2 expected/exomes_yri_mx2
} {}

test process_small {process_sample genome yri mx2} {
	cd $::bigtestdir
	set ref $::bigtestdir/refseqtest/hg19
	file delete -force tmp/one_genome_yri_mx2
	file mkdir tmp/one_genome_yri_mx2/samples/NA19240cgmx2
	mklink ori/genomes_yri_mx2.start/samples/NA19240cgmx2/ori tmp/one_genome_yri_mx2/samples/NA19240cgmx2/ori
	cg process_sample --stack 1 --verbose 2 {*}$::dopts -split 1 -dbdir $ref tmp/one_genome_yri_mx2/samples/NA19240cgmx2 2>@ stderr >@ stdout
	# check vs expected
	checkdiff -y --suppress-common-lines tmp/one_genome_yri_mx2/samples/NA19240cgmx2/summary-NA19240cgmx2.txt expected/genomes_yri_mx2/samples/NA19240cgmx2/summary-NA19240cgmx2.txt | grep -v "finished.*finished"
	checkdiff -qr -x info_analysis.tsv -x log_jobs -x summary-NA19240cgmx2.txt tmp/one_genome_yri_mx2/samples/NA19240cgmx2 expected/genomes_yri_mx2/samples/NA19240cgmx2
} {}

test process_small {genomes yri mx2} {
	cd $::bigtestdir	
	set dest tmp/genomes_yri_mx2
	file delete -force tmp/genomes_yri_mx2
	file mkdir tmp/genomes_yri_mx2
	foreach sample {
		NA19238cgmx2 NA19239cgmx2 NA19240cgmx2
		NA19240ilmx2
	} {
		file mkdir tmp/genomes_yri_mx2/samples/$sample
		mklink ori/genomes_yri_mx2.start/samples/$sample/ori tmp/genomes_yri_mx2/samples/$sample/ori
	}
	# cg process_project --stack 1 --verbose 2 -d 2 -split 1 -dbdir /complgen/refseq/testdb2/hg19 tmp/genomes_yri_mx2
	cg process_project --stack 1 --verbose 2 {*}$::dopts -split 1 -dbdir refseqtest/hg19 tmp/genomes_yri_mx2 2>@ stderr >@ stdout
	# check vs expected
	foreach cgsample {NA19238cgmx2 NA19239cgmx2 NA19240cgmx2} {
		checkdiff -y --suppress-common-lines tmp/genomes_yri_mx2/samples/$cgsample/summary-$cgsample.txt expected/genomes_yri_mx2/samples/$cgsample/summary-$cgsample.txt | grep -v "finished.*finished"
	}
	checkdiff -y --suppress-common-lines tmp/genomes_yri_mx2/samples/NA19240ilmx2/map-dsbwa-NA19240ilmx2.bam.dupmetrics expected/genomes_yri_mx2/samples/NA19240ilmx2/map-dsbwa-NA19240ilmx2.bam.dupmetrics | grep -v "Started on"
	checkdiff -qr -x *log_jobs -x *.bam -x *.bai -x *_fastqc -x summary-* -x *dupmetrics -x colinfo tmp/genomes_yri_mx2 expected/genomes_yri_mx2
	# cg tsvdiff -q 1 -x log_jobs -x _fastqc -x summary- -x dupmetrics -x colinfo tmp/genomes_yri_mx2 expected/genomes_yri_mx2 2> temp
	# cg tsvdiff -f 'chromosome begin end type ref alt zyg-*' -d kdiff3 tmp/genomes_yri_mx2/compar/annot_compar-genomes_yri_mx2.tsv expected/genomes_yri_mx2/compar/annot_compar-genomes_yri_mx2.tsv
} {}

testsummarize

if 0 {

# create test data
# ================

proc extractfromfastq {fastq result names} {
	unset -nocomplain a
	foreach name $names {
		set a([string range $name 0 end-2]) 1
	}
	set f [gzopen $fastq]
	set o [open [gzroot $result] w]
	while 1 {
		if {[gets $f name] == -1} break
		set seq $name
		append seq \n[gets $f]
		append seq \n[gets $f]
		append seq \n[gets $f]
		regsub {/[12]$} $name {} name
		if {[info exists a([lindex $name 0])]} {puts $o $seq}
	}
	close $o
	gzclose $f
	exec gzip -f [file root $result]
}

# mastr
# ------

cd /data/genomecomb.testdata
set src expected/mastr_116068_116083
set dest ori/mastr_mx2.start/samples
foreach sample {blanco2_8485 ceph1331_01_34_8452 ceph1331_02_34_8455 ceph1347_02_34_8446 ceph1347_02_34_7149 ceph1333_02_34_7220} {
	puts "---------- $sample ----------"
	set sdest $dest/${sample}mx2
	set bamfile $sdest/map-rsbwa-${sample}mx2.bam
	set fastqs [ssort -natural [glob $src/$sample*/fastq/*]]
	set refseq [glob /data/genomecomb.testdata/refseqtest/hg19/genome_*.ifas]
	file mkdir $sdest
	map_bwa_job $bamfile.pre $refseq $fastqs sample 1
	exec samtools view -F 0x100 -b $bamfile.pre "chr1:175087565-175087840" > $bamfile.temp1
	exec samtools view -F 0x100 -b $bamfile.pre "chr21:42732949-42781869" > $bamfile.temp2
	exec samtools view -F 0x100 -b $bamfile.pre "chr22:41920849-41921805" > $bamfile.temp3
	exec samtools merge -f $bamfile $bamfile.temp1 $bamfile.temp2 $bamfile.temp3
	file delete $bamfile.temp1 $bamfile.temp2 $bamfile.temp3 $bamfile.pre $bamfile.pre.bai
	file mkdir $sdest/ori
	file mkdir $sdest/tmp
	cg bam2fastq $bamfile $sdest/tmp/${sample}mx2_R1.fq.gz $sdest/tmp/${sample}mx2_R2.fq.gz
	set names [split [exec zcat $sdest/tmp/${sample}mx2_R1.fq.gz | grep ^@] \n]
	foreach fastq $fastqs result [list $sdest/ori/${sample}mx2_R1.fq.gz $sdest/ori/${sample}mx2_R2.fq.gz] {
		puts $result
		extractfromfastq $fastq $result $names
	}
}

# test
set samples {}
foreach s {blanco2_8485 ceph1333_02_34_7220 ceph1347_02_34_7149 ceph1347_02_34_8446} {
	lappend samples gatk-crsbwa-$s sam-crsbwa-$s
}
cg select -q {region("chr1:175087565-175087840") or region("chr21:42732949-42781869") or region("chr22:41921049-41921405")} expected/mastr_116068_116083/compar/annot_compar-mastr_116068_116083.tsv \
	| cg select -ssamples $samples \
	| cg select -q {scount($sequenced eq "v") > 0} > expected.tsv
cg select -ssamples $samples tmp/mastr_mx2/compar/annot_compar-mastr_mx2.tsv test.tsv
# cg tsvdiff -f 'chromosome begin end type ref alt zyg-*' test.tsv expected.tsv

# exomes
# ------
# extract small part of exome (TNN: chr1:175087565-175087840 , MX2 gene : chr21:42732949-42781869 and ACO2 chr22:41921149-41921305)
cd /data/genomecomb.testdata
set src ori/exomes_yri_chr2122.start/samples
set dest ori/exomes_yri_mx2.start/samples
foreach sample {NA19238 NA19239 NA19240} {
	puts $sample
	set sdest $dest/${sample}mx2
	set bamfile $sdest/map-rdsbwa-${sample}mx2.bam
	set fastqs [ssort -natural [glob $src/$sample*/fastq/*]]
	file mkdir $sdest
	exec samtools view -F 0x100 -b [glob $src/$sample*/*.bam] "chr21:42732949-42781869" > $bamfile.temp1
	exec samtools view -F 0x100 -b [glob $src/$sample*/*.bam] "chr22:41921049-41923951" > $bamfile.temp2
	exec samtools merge -f $bamfile $bamfile.temp1 $bamfile.temp2
	file delete $bamfile.temp1 $bamfile.temp2
	file mkdir $sdest/tmp
	cg bam2fastq $bamfile $sdest/tmp/${sample}mx2_R1.fq.gz $sdest/tmp/${sample}mx2_R2.fq.gz
	set names [split [exec zcat $sdest/tmp/${sample}mx2_R1.fq.gz | grep {^@.*/1$}] \n]
	foreach fastq $fastqs result [list $sdest/ori/${sample}mx2_R1.fq.gz $sdest/ori/${sample}mx2_R2.fq.gz] {
		puts $result
		extractfromfastq $fastq $result $names
	}
}

# test
cg select -q {region("chr21:42732949-42781869") or region("chr22:41921049-41923951")} expected/exomes_yri_chr2122/compar/annot_compar-exomes_yri_chr2122.tsv \
	| cg select -f {chromosome begin end type ref alt {*mx2=$*chr2122} homopolymer rmsk simpleRepeat snp135_name snp135_freq} > expected.tsv \
	| cg select -q {scount($sequenced eq "v") > 0} > expected.tsv
cg select tmp/exomes_yri_mx2/compar/annot_compar-exomes_yri_mx2.tsv test.tsv
ktdiff test.tsv expected.tsv
cg tsvdiff -t xl test.tsv expected.tsv

# illumina genome
# ---------------
cd /data/genomecomb.testdata
set src ori/genomes_yritrio_chr2122.start/samples/testNA19240chr21il.ori/NA19240_GAIIx_100_chr21.bam
set dest ori/genomes_yri_mx2.start/samples
set sample NA19240il
puts $sample
set sdest $dest/${sample}mx2
file mkdir $sdest
set bamfile $sdest/map-rsbwa-${sample}mx2.bam
exec samtools view -b $src "chr21:42732949-42781869" > $bamfile
file mkdir $sdest/ori
cg bam2fastq $bamfile $sdest/ori/${sample}mx2_R1.fq.gz $sdest/ori/${sample}mx2_R2.fq.gz

# cgi genomes
# -----------
set src ori/genomes_yri_chr2122.start/samples
set dest ori/genomes_yri_mx2.start/samples
foreach sample {NA19238 NA19239 NA19240} {
	puts $sample
	set ssample $src/test${sample}chr2122cg.ori/ASM
	set sdest $dest/${sample}cgmx2/ori/ASM
	file mkdir $sdest
	catch {file copy -force $ssample/CNV $sdest}
	catch {file copy -force $ssample/SV $sdest}
	file mkdir $sdest/REF
	set reffile [glob $ssample/REF/coverageRefScore-chr21-*-ASM.tsv.bz2]
	exec cg select -q {$offset >= 42732949 and $offset < 42781869} $reffile | bzip2 > $sdest/REF/[file tail $reffile]
	set reffile [glob $ssample/REF/coverageRefScore-chr22-*-ASM.tsv.bz2]
	exec cg select -q {$offset >= 41921049 and $offset < 41921405} $reffile | bzip2 > $sdest/REF/[file tail $reffile]
	#
	set varfile [glob $ssample/var-*.tsv.bz2]
	exec cg select -q {region("chr21:42732949-42781869") or region("chr22:41921049-41921405")} $varfile | bzip2 > $sdest/[file tail $varfile]
	set genefile [glob $ssample/gene-*.tsv.bz2]
	exec cg select -q {region("chr21:42732949-42781869") or region("chr22:41921049-41921405")} $genefile | bzip2 > $sdest/[file tail $genefile]
}

set samples {cg-cg-NA19238cgmx2 cg-cg-NA19239cgmx2 cg-cg-NA19240cgmx2 sam-rdsbwa-NA19240ilmx2 gatk-rdsbwa-NA19240ilmx2}
cg select -q {region("chr21:42732949-42781869") or region("chr22:41921049-41921405")} expected/genomes_yri_chr2122/compar/annot_compar-genomes_yri_chr2122.tsv \
	| cg select -f {chromosome begin end type ref alt {*-**cgmx2=$*-test**chr2122cg} {*-**ilmx2=$*-test**chr21il} homopolymer rmsk simpleRepeat snp135_name snp135_freq} > expected.tsv \
	| cg select -ssamples $samples \
	| cg select -q {scount($sequenced eq "v") > 0} > expected.tsv
cg select tmp/genomes_yri_mx2/compar/annot_compar-genomes_yri_mx2.tsv test.tsv
ktdiff test.tsv expected.tsv
cg tsvdiff -d kdiff3 -t xl test.tsv expected.tsv

}
