# create small test data
# ======================

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
#set src ori/genomes_yritrio_chr2122.start/samples/testNA19240chr21il.ori/NA19240_GAIIx_100_chr21.bam
set src tmp/genomes_yri_chr2122/samples/testNA19240chr21il/map-rdsbwa-testNA19240chr21il.bam
set dest ori/genomes_yri_mx2.start/samples
set sample NA19240il
puts $sample
set sdest $dest/${sample}mx2
file mkdir $sdest
set bamfile $sdest/map-rsbwa-${sample}mx2.bam
exec samtools view -b $src "chr21:42732949-42781869" > $bamfile
samtools index $bamfile

file mkdir $sdest/ori
cg bam2fastq $bamfile $sdest/ori/${sample}mx2_R1.fq.gz $sdest/ori/${sample}mx2_R2.fq.gz


cd /data/genomecomb.testdata
set src tmp/genomes_yri_chr2122/samples/testNA19240chr21il/map-rdsbwa-testNA19240chr21il.bam
set bamfile tmp/map-rsbwa-NA19240ilmx2.bam
exec samtools view -b $src "chr21:42732949-42781869" > $bamfile

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
