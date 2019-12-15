#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test realign {realign_gatk basic} {
	exec samtools view -b data/bwa.sam > tmp/bwa.bam
	cg realign_gatk -stack 1 tmp/bwa.bam tmp/ratest.bam $::refseqdir/hg19
	exec samtools view tmp/ratest.bam > tmp/ratest.sam
	catch {exec diff tmp/ratest.sam data/ratest-gatk.sam}
} 0

test realign {realign_gatk pipe} {
	exec samtools view -b data/bwa.sam > tmp/bwa.bam
	cg realign_gatk -stack 1 -refseq $::refseqdir/hg19 < tmp/bwa.bam > tmp/ratest.bam
	exec samtools view tmp/ratest.bam > tmp/ratest.sam
	catch {exec diff tmp/ratest.sam data/ratest-gatk.sam}
} 0

test realign {realign -method gatk pipe} {
	exec samtools view -b data/bwa.sam > tmp/bwa.bam
	cg realign -method gatk -stack 1 -refseq $::refseqdir/hg19 < tmp/bwa.bam > tmp/ratest.bam
	exec samtools view tmp/ratest.bam > tmp/ratest.sam
	catch {exec diff tmp/ratest.sam data/ratest-gatk.sam}
} 0

test realign {realign_abra basic} {
	exec samtools view -b data/bwa.sam > tmp/bwa.bam
	cg realign_abra -stack 1 tmp/bwa.bam tmp/ratest.bam $::refseqdir/hg19
	cg sam2tsv tmp/ratest.bam tmp/ratest.tsv
	cg sam2tsv data/ratest-abra.sam tmp/expected.tsv
	catch {exec diff tmp/ratest.tsv tmp/expected.tsv}
} 0

test realign {realign_abra pipe} {
	exec samtools view -b data/bwa.sam > tmp/bwa.bam
	cg realign_abra -stack 1 -refseq $::refseqdir/hg19 < tmp/bwa.bam > tmp/ratest.bam
	cg sam2tsv tmp/ratest.bam tmp/ratest.tsv
	cg sam2tsv data/ratest-abra.sam tmp/expected.tsv
	catch {exec diff tmp/ratest.tsv tmp/expected.tsv}
} 0

test realign {realign_abra from compressed sam} {
	file copy data/bwa.sam tmp/bwa.sam
	cg zst tmp/bwa.sam
	cg realign_abra -stack 1 tmp/bwa.sam.zst tmp/ratest.bam $::refseqdir/hg19
	cg sam2tsv tmp/ratest.bam tmp/ratest.tsv
	cg sam2tsv data/ratest-abra.sam tmp/expected.tsv
	catch {exec diff tmp/ratest.tsv tmp/expected.tsv}
} 0

test realign {realign_srma basic} {
	exec samtools view -b data/bwa.sam > tmp/bwa.bam
	cg realign_srma -stack 1 tmp/bwa.bam tmp/ratest.bam $::refseqdir/hg19
	cg sam2tsv tmp/ratest.bam tmp/ratest.tsv
	cg sam2tsv data/ratest-srma.sam tmp/expected.tsv
	catch {exec diff tmp/ratest.tsv tmp/expected.tsv}
} 0

test realign {realign_srma pipe} {
	exec samtools view -b data/bwa.sam > tmp/bwa.bam
	cg realign_srma -stack 1 -refseq $::refseqdir/hg19 < tmp/bwa.bam > tmp/ratest.bam
	cg sam2tsv tmp/ratest.bam tmp/ratest.tsv
	cg sam2tsv data/ratest-srma.sam tmp/expected.tsv
	catch {exec diff tmp/ratest.tsv tmp/expected.tsv}
} 0

test markdup {bam_markduplicates picard} {
	exec samtools sort data/bwa.sam > tmp/sbwa.bam
	cg bam_markduplicates -method picard tmp/sbwa.bam tmp/result.bam
	exec cg sam2tsv tmp/result.bam | cg select -f {qname chromosome begin end duplicate} > tmp/result.tsv
	exec cg sam2tsv data/dsbwa.sam | cg select -f {qname chromosome begin end duplicate} > tmp/dsbwa.tsv
	catch {exec cg tsvdiff tmp/result.tsv tmp/dsbwa.tsv}
} 0

test markdup {bam_markduplicates picard pipe} {
	exec samtools sort data/bwa.sam > tmp/sbwa.bam
	exec cg bam_markduplicates -stack 1 -method picard < tmp/sbwa.bam > tmp/result.bam
	exec cg sam2tsv tmp/result.bam | cg select -f {qname chromosome begin end duplicate} > tmp/result.tsv
	exec cg sam2tsv data/dsbwa.sam | cg select -f {qname chromosome begin end duplicate} > tmp/dsbwa.tsv
	catch {exec cg tsvdiff tmp/result.tsv tmp/dsbwa.tsv}
} 0

test markdup {bam_markduplicates samtools} {
	exec samtools sort data/bwa.sam > tmp/sbwa.bam
	cg bam_markduplicates -method samtools tmp/sbwa.bam tmp/result.bam
	exec cg tsvdiff tmp/result.bam data/dsbwa.sam
} {}

test markdup {bam_markduplicates samtools pipe} {
	exec samtools sort data/bwa.sam > tmp/sbwa.bam
	cg bam_markduplicates -method samtools < tmp/sbwa.bam > tmp/result.bam
	exec cg tsvdiff tmp/result.bam data/dsbwa.sam
} {}

test markdup {bam_markduplicates samtools pipe -compressionlevel} {
	exec samtools sort data/bwa.sam > tmp/sbwa.bam
	cg bam_markduplicates -method samtools -compressionlevel 1 < tmp/sbwa.bam > tmp/result.bam
	exec cg tsvdiff tmp/result.bam data/dsbwa.sam
} {}

test markdup {bam_markduplicates samtools pipe to cram} {
	exec samtools sort data/bwa.sam > tmp/sbwa.bam
	cg bam_markduplicates -method samtools -outputformat cram -refseq $::refseqdir/hg19 < tmp/sbwa.bam > tmp/result.cram
	exec samtools view -h -b -T [refseq $::refseqdir/hg19] tmp/result.cram > tmp/result.bam
	exec cg sam2tsv tmp/result.bam | cg select -f {qname chromosome begin end duplicate} > tmp/result.tsv
	exec cg sam2tsv data/dsbwa.sam | cg select -f {qname chromosome begin end duplicate} > tmp/dsbwa.tsv
	catch {exec cg tsvdiff tmp/result.tsv tmp/dsbwa.tsv}
} 0

test markdup {bam_markduplicates biobambam} {
	exec samtools sort data/bwa.sam > tmp/sbwa.bam
	cg bam_markduplicates -method biobambam tmp/sbwa.bam tmp/result.bam
	exec cg tsvdiff tmp/result.bam data/dsbwa.sam
} {}

test markdup {bam_markduplicates biobambam pipe} {
	exec samtools sort data/bwa.sam > tmp/sbwa.bam
	cg bam_markduplicates -method biobambam < tmp/sbwa.bam > tmp/result.bam
	exec cg tsvdiff tmp/result.bam data/dsbwa.sam
} {}

test markdup {bam_clean} {
	test_cleantmp
	file copy -force -- data/bwa.sam tmp/bwa.sam
	cg bam_clean -stack 1 -keep 1 -refseq $::refseqdir/hg19 -sort 1 -removeduplicates 1 -realign 1 tmp/bwa.sam
	set c [string trim [file_read tmp/rdsbwa.bam.analysisinfo]]
	if {$c ne [string trim [deindent {
		removeduplicates	removeduplicates_version	realign	realign_version
		samtools	1.9 (using htslib 1.9)	gatk	3.7-0-gcfedb67
	}]]} {
		error "error in analysisinfo"
	}
	exec cg sam2tsv tmp/rdsbwa.bam | cg select -f {qname chromosome begin end duplicate} > tmp/result.tsv
	exec cg sam2tsv data/dsbwa.sam | cg select -f {qname chromosome begin end duplicate} > tmp/sbwa.tsv
	exec cg tsvdiff tmp/result.tsv tmp/sbwa.tsv
} {diff tmp/result.tsv tmp/sbwa.tsv
header
  qname	chromosome	begin	end	duplicate
147d146
< SRR792091.1611898	chr21	42779842	42779942	1
148a148
> SRR792091.1611898	chr21	42779842	42779942	1
155d154
< SRR792091.1611898	chr21	42779960	42780060	1
156a156
> SRR792091.1611898	chr21	42779960	42780060	1
child process exited abnormally} error

testsummarize
