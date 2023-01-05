package require Extral
catch {tk appname test}

set bigtestdir $env(HOME)/genomecomb.testdata
set smalltestdir $env(HOME)/genomecomb.smalltestdata
set publicdir $env(HOME)/public
if {![info exists refseqdir]} {
	set refseqdir $smalltestdir/refseqtest
}

if {[info exists argv]} {
	set pos [lsearch $argv -smalltestdir]
	if {$pos != -1} {
		set smalltestdir [lindex $argv [expr {$pos + 1}]]
		set argv [lreplace $argv $pos [expr {$pos + 1}]]
	}
	set dopts $argv
} else {
	set dopts {--stack 1 --verbose 2}
}

package require pkgtools
namespace import -force pkgtools::*
package require Extral

set test_cleantmp 1

# pkgtools::testleak 100

set keeppath $::env(PATH)
set script [info script] ; if {$script eq ""} {set script ./t}
set appdir [file dir [file dir [file normalize $script]]]
if {![info exists basetestdir]} {
	# set basetestdir [file dir [file normalize $script]]
	set basetestdir [file dir $appdir]/tests_genomecomb
}
file mkdir $basetestdir
cd $basetestdir
set testdir $basetestdir

puts "Using as testdir: $basetestdir"

if {![file exists $basetestdir/data]} {
	exec ln -sf $appdir/tests/data .
}
if {![file exists $basetestdir/tmp]} {
	file mkdir $basetestdir/tmp
}

lappend auto_path $appdir/lib $appdir/lib-exp
source $appdir/lib/file.tcl ; pwd
package require genomecomb
if {![info exists ::genomecombdir]} {genomecombenv}

proc testdir {args} {
	set testname [join [lrange $args 0 1] __]
	set testname [string_change $testname {{ } _ : _ / _ \\ _ \; _ * _ ? _ \} _ \{ _ \n _ \t _ \[ _ \] _ ( _ ) _}]
	return $::basetestdir/$testname
}

proc test {args} {
	set numargs [llength $args]
	if {$numargs == 0} {
		set ::testdir $::appdir/tests
		cd $::appdir/tests
		puts "testdir set back to [pwd]"
		return
	} elseif {$numargs == 2 || ($numargs == 3 && [string trim [lindex $args 2]] eq "")} {
		set ::testdir [testdir {*}$args]
		file mkdir $::testdir
		file mkdir $::testdir/tmp
		puts "testdir set to [pwd]"
		cd $::testdir
		return
	} elseif {$numargs < 4} {
		error "wrong # parameters, format is : test group description script expected ..."
		group description script expected args
	}
	set ::testdir [testdir {*}$args]
	file mkdir $::testdir
	file mkdir $::testdir/tmp
	cd $::testdir
	exec ln -sf $::appdir/tests/data .
	if {[get ::test_cleantmp 1]} {test_cleantmp}
	catch {job_init}
	set description [lindex $args 1]
	append description " ($::testdir)"
	lset args 1 $description
	pkgtools::test {*}$args
	set ::testdir $::appdir/tests
	cd $::appdir/tests
	return {}
}

proc tsvdiff {args} {
	if {[catch {cg tsvdiff {*}$args} result]} {
		regsub "child process exited abnormally\n?" $result {} result
		return $result
	} else {
		return ""
	}
}

proc checkdiff args {
	global e
	set err [catch {exec diff {*}$args} e]
	if {$err && $e ne {child process exited abnormally}} {
		# set pos 0; foreach v $args {if {[string index $v 0] ne "-"} break; incr pos}
		set temp [lindex [split $args |] 0]
		return "Files differ: [lrange $temp end-1 end]"
	} else {
		return ""
	}
}

proc test_cleantmp {} {
	foreach file [list_remove [glob -nocomplain $::testdir/tmp/* $::testdir/tmp/.*] $::testdir/tmp/.. $::testdir/tmp/.] {
		catch {file attributes $file -permissions ugo+xw}
		catch {file delete -force $file}
	}
	foreach file [list_remove [glob -nocomplain tmp/* tmp/.*] tmp/.. tmp/.] {
		catch {file attributes $file -permissions ugo+xw}
		catch {file delete -force $file}
	}
	cg indexclean
}

proc write_tab {file data {comment {}}} {
	set data [split [string trim $data] \n]
	set f [open $file w]
	if {$comment ne ""} {puts $f $comment}
	foreach line $data {
		puts $f [join $line \t]
	}
	close $f
}

proc write_deindent {file data} {
	file_write $file [deindent $data]\n
}

proc diff_tab {file data {comment {}}} {
	set tempfile [tempfile]
	set data [split [string trim $data] \n]
	set f [open $tempfile w]
	if {$comment ne ""} {puts $f $comment}
	foreach line $data {
		puts $f [join $line \t]
	}
	close $f
}

proc test_genomecombdir {} {
	set expdir tmp/test
	file mkdir $expdir/samples
	file mkdir $expdir/compar
	set samplesdir $expdir/samples
	foreach n {1 2 3} {
		file mkdir $samplesdir/sample$n
	}
	cg splitalleles data/var_annot.tsv > $samplesdir/sample1/var-sample1.tsv
	cg splitalleles data/var_annot2.tsv > $samplesdir/sample2/var-sample2.tsv
	cg splitalleles data/var_annot2seq.tsv > $samplesdir/sample3/prevar-sample3.tsv
	cg select -f {sequenced *} $samplesdir/sample3/prevar-sample3.tsv $samplesdir/sample3/var-sample3.tsv
	file copy data/sreg-annot1.tsv $samplesdir/sample1/sreg-sample1.tsv
	file copy data/sreg-annot2.tsv $samplesdir/sample2/sreg-sample2.tsv
	file copy data/sreg-annot2.tsv $samplesdir/sample3/sreg-sample3.tsv
	cg multicompar -reannot 1 -split 1 $expdir/compar/compar-test.tsv {*}[bsort [glob $samplesdir/*/var-*.tsv]]
	file delete $expdir/compar/compar-test.tsv.reannot
	file delete $expdir/compar/compar-test.tsv.old
	# exec diff $expdir/compar/compar-test.tsv data/expected-multicompar-split-reannot.tsv
}

proc file_regsub {exp subSpec file resultfile} {
	set f [open $file]
	set o [open $resultfile w]
	while {[gets $f line] != -1} {
		regsub -all $exp $line $subSpec line
		puts $o $line
	}
	close $o
	close $f
}

set sam_header {
	@HD	VN:1.4	GO:none	SO:coordinate
	@SQ	SN:chr1	LN:249250621
	@SQ	SN:chr2	LN:243199373
	@SQ	SN:chr3	LN:198022430
	@SQ	SN:chr4	LN:191154276
	@SQ	SN:chr5	LN:180915260
	@SQ	SN:chr6	LN:171115067
	@SQ	SN:chr7	LN:159138663
	@SQ	SN:chr8	LN:146364022
	@SQ	SN:chr9	LN:141213431
	@SQ	SN:chr10	LN:135534747
	@SQ	SN:chr11	LN:135006516
	@SQ	SN:chr12	LN:133851895
	@SQ	SN:chr13	LN:115169878
	@SQ	SN:chr14	LN:107349540
	@SQ	SN:chr15	LN:102531392
	@SQ	SN:chr16	LN:90354753
	@SQ	SN:chr17	LN:81195210
	@SQ	SN:chr18	LN:78077248
	@SQ	SN:chr19	LN:59128983
	@SQ	SN:chr20	LN:63025520
	@SQ	SN:chr21	LN:48129895
	@SQ	SN:chr22	LN:51304566
	@SQ	SN:chrM	LN:16571
	@SQ	SN:chrX	LN:155270560
	@SQ	SN:chrY	LN:59373566
	@RG	ID:sample1	PL:illumina	PU:sample1	LB:solexa-123	SM:sample1
	@PG	ID:GATK IndelRealigner	VN:2.4-9-g532efad	CL:knownAlleles=[] targetIntervals=test.intervals LODThresholdForCleaning=5.0 consensusDeterminationModel=USE_READS entropyThreshold=0.15 maxReadsInMemory=150000 maxIsizeForMovement=3000 maxPositionalMoveAllowed=200 maxConsensuses=30 maxReadsForConsensuses=120 maxReadsForRealignment=20000 noOriginalAlignmentTags=false nWayOut=null generate_nWayOut_md5s=false check_early=false noPGTag=false keepPGTags=false indelsFileForDebugging=null statisticsFileForDebugging=null SNPsFileForDebugging=null
}

proc write_sam {file data {namebase A}} {
	# creates a sam file based on data
	# sequences are all one base (given in data or A) and qualities are all -
	set tempfile [tempfile]
	set o [open $tempfile w]
	set num 1
	foreach line [split [string trim $data] \n] {
		set base A
		foreach {chr1 pos1 cigar1 seq1 chr2 pos2 cigar2 seq2 base} $line break
		if {$base eq ""} {set base A}
		if {[isint $seq1]} {
			set size1 $seq1
			set seq1 [string_fill $base $size1]
		} else {
			set size1 [string length $seq1]
		}
		set qual1 [string_fill - $size1]
		if {$seq2 ne ""} {
			if {[isint $seq2]} {
				set size2 $seq2
				set seq2 [string_fill $base $size2]
			} else {
				set size2 [string length $seq2]
			}
			set qual2 [string_fill - $size2]
			if {$chr2 eq $chr1} {set c2 =} else {set c2 $chr2}
			set tlen [expr {$pos2+$size2-$pos1}]
			set flags 99
		} else {
			set c2 *
			set pos2 0
			set tlen 0
			set flags 16
		}
		puts $o [join [list $namebase$num $flags $chr1 $pos1 60 $cigar1 $c2 $pos2 $tlen $seq1 $qual1 RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25] \t]
		if {$seq2 ne ""} {
			if {$chr2 eq $chr1} {set c1 =} else {set c1 $chr1}
			puts $o [join [list $namebase$num 147 $chr2 $pos2 60 $cigar2 $c1 $pos1 -[expr {$pos2+$size2-$pos1}] $seq2 $qual2 RG:Z:sample1	MQ:i:60	AS:i:241	XS:i:25] \t]
		}
		incr num
	}
	close $o
	set o [open $file w]
	puts $o [deindent $::sam_header]
	close $o
	exec gnusort8 -t \t -N -s -k3,3 -k4,4 -k1,1 -k2,2 $tempfile >> $file
}

proc write_vcf {file data {extracomment {}} {extrainfo {}}} {
	set data [split [string trim $data] \n]
	set f [open $file w]
	if {$extrainfo ne ""} {
		regsub -all {\n\t*} [string trim $extrainfo] "\n\t\t" temp
		set extrainfo \n\t\t$temp
	}
	puts $f [deindent [subst {
		##fileformat=VCFv4.0
		##fileDate=20090805
		##source=myImputationProgramV3.1
		##reference=1000GenomesPilot-NCBI36
		##phasing=partial
		##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
		##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
		##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
		##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
		##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
		##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">$extrainfo
		##FILTER=<ID=q10,Description="Quality below 10">
		##FILTER=<ID=s50,Description="Less than 50% of samples have data">
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
		##FORMAT=<ID=TE,Number=A,Type=Integer,Description="test for alt alleles in the order listed">
		##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
		##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
		##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
	}]]
	if {$extracomment ne ""} {puts -nonewline $f $extracomment}
	set header [lindex $data 0]
	set data [lrange $data 1 end]
	puts $f \#[join $header \t]
	foreach line $data {
		puts $f [join $line \t]
	}
	close $f
}

proc bgvcfheader {} {
	deindent {
		##fileformat=VCFv4.2
		##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">
		##FILTER=<ID=LowQual,Description="Low quality">
		##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
		##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
		##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
		##FORMAT=<ID=GQX,Number=1,Type=Integer,Description="Empirically calibrated genotype quality score for variant sites, otherwise minimum of {Genotype quality assuming variant position,Genotype quality assuming non-variant position}">
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
		##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another">
		##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
		##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
		##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">
		##GATKCommandLine=<ID=HaplotypeCaller,CommandLine="HaplotypeCaller  --annotation-group StandardAnnotation --annotation-group AS_StandardAnnotation --emit-ref-confidence GVCF --annotate-with-num-discovered-alleles true --output varall-gatkh-rdsbwa-NA19240chr2122.gvcf.gz.temp.gz --intervals /tmp/tempExtral.10600-UMRUbdTebwJFtLcIEVGB/_Extral_temp_1.bed --input /data/genomecomb.testdata/tmp/exomes_yri_chr2122/samples/NA19240chr2122/map-rdsbwa-NA19240chr2122.bam --reference /complgen/refseq/hg19/genome_hg19.fa  --disable-tool-default-annotations false --gvcf-gq-bands 1 --gvcf-gq-bands 2 --gvcf-gq-bands 3 --gvcf-gq-bands 4 --gvcf-gq-bands 5 --gvcf-gq-bands 6 --gvcf-gq-bands 7 --gvcf-gq-bands 8 --gvcf-gq-bands 9 --gvcf-gq-bands 10 --gvcf-gq-bands 11 --gvcf-gq-bands 12 --gvcf-gq-bands 13 --gvcf-gq-bands 14 --gvcf-gq-bands 15 --gvcf-gq-bands 16 --gvcf-gq-bands 17 --gvcf-gq-bands 18 --gvcf-gq-bands 19 --gvcf-gq-bands 20 --gvcf-gq-bands 21 --gvcf-gq-bands 22 --gvcf-gq-bands 23 --gvcf-gq-bands 24 --gvcf-gq-bands 25 --gvcf-gq-bands 26 --gvcf-gq-bands 27 --gvcf-gq-bands 28 --gvcf-gq-bands 29 --gvcf-gq-bands 30 --gvcf-gq-bands 31 --gvcf-gq-bands 32 --gvcf-gq-bands 33 --gvcf-gq-bands 34 --gvcf-gq-bands 35 --gvcf-gq-bands 36 --gvcf-gq-bands 37 --gvcf-gq-bands 38 --gvcf-gq-bands 39 --gvcf-gq-bands 40 --gvcf-gq-bands 41 --gvcf-gq-bands 42 --gvcf-gq-bands 43 --gvcf-gq-bands 44 --gvcf-gq-bands 45 --gvcf-gq-bands 46 --gvcf-gq-bands 47 --gvcf-gq-bands 48 --gvcf-gq-bands 49 --gvcf-gq-bands 50 --gvcf-gq-bands 51 --gvcf-gq-bands 52 --gvcf-gq-bands 53 --gvcf-gq-bands 54 --gvcf-gq-bands 55 --gvcf-gq-bands 56 --gvcf-gq-bands 57 --gvcf-gq-bands 58 --gvcf-gq-bands 59 --gvcf-gq-bands 60 --gvcf-gq-bands 70 --gvcf-gq-bands 80 --gvcf-gq-bands 90 --gvcf-gq-bands 99 --indel-size-to-eliminate-in-ref-model 10 --use-alleles-trigger false --disable-optimizations false --just-determine-active-regions false --dont-genotype false --dont-trim-active-regions false --max-disc-ar-extension 25 --max-gga-ar-extension 300 --padding-around-indels 150 --padding-around-snps 20 --kmer-size 10 --kmer-size 25 --dont-increase-kmer-sizes-for-cycles false --allow-non-unique-kmers-in-ref false --num-pruning-samples 1 --recover-dangling-heads false --do-not-recover-dangling-branches false --min-dangling-branch-length 4 --consensus false --max-num-haplotypes-in-population 128 --error-correct-kmers false --min-pruning 2 --debug-graph-transformations false --kmer-length-for-read-error-correction 25 --min-observations-for-kmer-to-be-solid 20 --likelihood-calculation-engine PairHMM --base-quality-score-threshold 18 --pair-hmm-gap-continuation-penalty 10 --pair-hmm-implementation FASTEST_AVAILABLE --pcr-indel-model CONSERVATIVE --phred-scaled-global-read-mismapping-rate 45 --native-pair-hmm-threads 4 --native-pair-hmm-use-double-precision false --debug false --use-filtered-reads-for-annotations false --bam-writer-type CALLED_HAPLOTYPES --dont-use-soft-clipped-bases false --capture-assembly-failure-bam false --error-correct-reads false --do-not-run-physical-phasing false --min-base-quality-score 10 --smith-waterman JAVA --use-new-qual-calculator false --heterozygosity 0.001 --indel-heterozygosity 1.25E-4 --heterozygosity-stdev 0.01 --standard-min-confidence-threshold-for-calling 10.0 --max-alternate-alleles 6 --max-genotype-count 1024 --sample-ploidy 2 --genotyping-mode DISCOVERY --contamination-fraction-to-filter 0.0 --output-mode EMIT_VARIANTS_ONLY --all-site-pls false --min-assembly-region-size 50 --max-assembly-region-size 300 --assembly-region-padding 100 --max-reads-per-alignment-start 50 --active-probability-threshold 0.002 --max-prob-propagation-distance 50 --interval-set-rule UNION --interval-padding 0 --interval-exclusion-padding 0 --interval-merging-rule ALL --read-validation-stringency SILENT --seconds-between-progress-updates 10.0 --disable-sequence-dictionary-validation false --create-output-bam-index true --create-output-bam-md5 false --create-output-variant-index true --create-output-variant-md5 false --lenient false --add-output-sam-program-record true --add-output-vcf-command-line true --cloud-prefetch-buffer 40 --cloud-index-prefetch-buffer -1 --disable-bam-index-caching false --help false --version false --showHidden false --verbosity INFO --QUIET false --use-jdk-deflater false --use-jdk-inflater false --gcs-max-retries 20 --disable-tool-default-read-filters false --minimum-mapping-quality 20",Version=4.0.3.0,Date="April 19, 2018 5:09:59 PM CEST">
		##GVCFBlock0-1=minGQ=0(inclusive),maxGQ=1(exclusive)
		##GVCFBlock1-2=minGQ=1(inclusive),maxGQ=2(exclusive)
		##GVCFBlock10-11=minGQ=10(inclusive),maxGQ=11(exclusive)
		##GVCFBlock11-12=minGQ=11(inclusive),maxGQ=12(exclusive)
		##GVCFBlock12-13=minGQ=12(inclusive),maxGQ=13(exclusive)
		##GVCFBlock13-14=minGQ=13(inclusive),maxGQ=14(exclusive)
		##GVCFBlock14-15=minGQ=14(inclusive),maxGQ=15(exclusive)
		##GVCFBlock15-16=minGQ=15(inclusive),maxGQ=16(exclusive)
		##GVCFBlock16-17=minGQ=16(inclusive),maxGQ=17(exclusive)
		##GVCFBlock17-18=minGQ=17(inclusive),maxGQ=18(exclusive)
		##GVCFBlock18-19=minGQ=18(inclusive),maxGQ=19(exclusive)
		##GVCFBlock19-20=minGQ=19(inclusive),maxGQ=20(exclusive)
		##GVCFBlock2-3=minGQ=2(inclusive),maxGQ=3(exclusive)
		##GVCFBlock20-21=minGQ=20(inclusive),maxGQ=21(exclusive)
		##GVCFBlock21-22=minGQ=21(inclusive),maxGQ=22(exclusive)
		##GVCFBlock22-23=minGQ=22(inclusive),maxGQ=23(exclusive)
		##GVCFBlock23-24=minGQ=23(inclusive),maxGQ=24(exclusive)
		##GVCFBlock24-25=minGQ=24(inclusive),maxGQ=25(exclusive)
		##GVCFBlock25-26=minGQ=25(inclusive),maxGQ=26(exclusive)
		##GVCFBlock26-27=minGQ=26(inclusive),maxGQ=27(exclusive)
		##GVCFBlock27-28=minGQ=27(inclusive),maxGQ=28(exclusive)
		##GVCFBlock28-29=minGQ=28(inclusive),maxGQ=29(exclusive)
		##GVCFBlock29-30=minGQ=29(inclusive),maxGQ=30(exclusive)
		##GVCFBlock3-4=minGQ=3(inclusive),maxGQ=4(exclusive)
		##GVCFBlock30-31=minGQ=30(inclusive),maxGQ=31(exclusive)
		##GVCFBlock31-32=minGQ=31(inclusive),maxGQ=32(exclusive)
		##GVCFBlock32-33=minGQ=32(inclusive),maxGQ=33(exclusive)
		##GVCFBlock33-34=minGQ=33(inclusive),maxGQ=34(exclusive)
		##GVCFBlock34-35=minGQ=34(inclusive),maxGQ=35(exclusive)
		##GVCFBlock35-36=minGQ=35(inclusive),maxGQ=36(exclusive)
		##GVCFBlock36-37=minGQ=36(inclusive),maxGQ=37(exclusive)
		##GVCFBlock37-38=minGQ=37(inclusive),maxGQ=38(exclusive)
		##GVCFBlock38-39=minGQ=38(inclusive),maxGQ=39(exclusive)
		##GVCFBlock39-40=minGQ=39(inclusive),maxGQ=40(exclusive)
		##GVCFBlock4-5=minGQ=4(inclusive),maxGQ=5(exclusive)
		##GVCFBlock40-41=minGQ=40(inclusive),maxGQ=41(exclusive)
		##GVCFBlock41-42=minGQ=41(inclusive),maxGQ=42(exclusive)
		##GVCFBlock42-43=minGQ=42(inclusive),maxGQ=43(exclusive)
		##GVCFBlock43-44=minGQ=43(inclusive),maxGQ=44(exclusive)
		##GVCFBlock44-45=minGQ=44(inclusive),maxGQ=45(exclusive)
		##GVCFBlock45-46=minGQ=45(inclusive),maxGQ=46(exclusive)
		##GVCFBlock46-47=minGQ=46(inclusive),maxGQ=47(exclusive)
		##GVCFBlock47-48=minGQ=47(inclusive),maxGQ=48(exclusive)
		##GVCFBlock48-49=minGQ=48(inclusive),maxGQ=49(exclusive)
		##GVCFBlock49-50=minGQ=49(inclusive),maxGQ=50(exclusive)
		##GVCFBlock5-6=minGQ=5(inclusive),maxGQ=6(exclusive)
		##GVCFBlock50-51=minGQ=50(inclusive),maxGQ=51(exclusive)
		##GVCFBlock51-52=minGQ=51(inclusive),maxGQ=52(exclusive)
		##GVCFBlock52-53=minGQ=52(inclusive),maxGQ=53(exclusive)
		##GVCFBlock53-54=minGQ=53(inclusive),maxGQ=54(exclusive)
		##GVCFBlock54-55=minGQ=54(inclusive),maxGQ=55(exclusive)
		##GVCFBlock55-56=minGQ=55(inclusive),maxGQ=56(exclusive)
		##GVCFBlock56-57=minGQ=56(inclusive),maxGQ=57(exclusive)
		##GVCFBlock57-58=minGQ=57(inclusive),maxGQ=58(exclusive)
		##GVCFBlock58-59=minGQ=58(inclusive),maxGQ=59(exclusive)
		##GVCFBlock59-60=minGQ=59(inclusive),maxGQ=60(exclusive)
		##GVCFBlock6-7=minGQ=6(inclusive),maxGQ=7(exclusive)
		##GVCFBlock60-70=minGQ=60(inclusive),maxGQ=70(exclusive)
		##GVCFBlock7-8=minGQ=7(inclusive),maxGQ=8(exclusive)
		##GVCFBlock70-80=minGQ=70(inclusive),maxGQ=80(exclusive)
		##GVCFBlock8-9=minGQ=8(inclusive),maxGQ=9(exclusive)
		##GVCFBlock80-90=minGQ=80(inclusive),maxGQ=90(exclusive)
		##GVCFBlock9-10=minGQ=9(inclusive),maxGQ=10(exclusive)
		##GVCFBlock90-99=minGQ=90(inclusive),maxGQ=99(exclusive)
		##GVCFBlock99-100=minGQ=99(inclusive),maxGQ=100(exclusive)
		##INFO=<ID=AS_InbreedingCoeff,Number=A,Type=Float,Description="allele specific heterozygosity as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation; relate to inbreeding coefficient">
		##INFO=<ID=AS_QD,Number=A,Type=Float,Description="Allele-specific Variant Confidence/Quality by Depth">
		##INFO=<ID=AS_RAW_BaseQRankSum,Number=1,Type=String,Description="raw data for allele specific rank sum test of base qualities">
		##INFO=<ID=AS_RAW_MQ,Number=1,Type=String,Description="Allele-specfic raw data for RMS Mapping Quality">
		##INFO=<ID=AS_RAW_MQRankSum,Number=1,Type=String,Description="Allele-specfic raw data for Mapping Quality Rank Sum">
		##INFO=<ID=AS_RAW_ReadPosRankSum,Number=1,Type=String,Description="allele specific raw data for rank sum test of read position bias">
		##INFO=<ID=AS_SB_TABLE,Number=1,Type=String,Description="Allele-specific forward/reverse read counts for strand bias tests">
		##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities">
		##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
		##INFO=<ID=DS,Number=0,Type=Flag,Description="Were any of the samples downsampled?">
		##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
		##INFO=<ID=ExcessHet,Number=1,Type=Float,Description="Phred-scaled p-value for exact test of excess heterozygosity">
		##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description="Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation">
		##INFO=<ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
		##INFO=<ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
		##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
		##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities">
		##INFO=<ID=NDA,Number=1,Type=Integer,Description="Number of alternate alleles discovered (but not necessarily genotyped) at this site">
		##INFO=<ID=RAW_MQ,Number=1,Type=Float,Description="Raw data for RMS Mapping Quality">
		##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias">
		##contig=<ID=chr1,length=249250621>
		##contig=<ID=chr2,length=243199373>
		##contig=<ID=chr3,length=198022430>
		##contig=<ID=chr4,length=191154276>
		##contig=<ID=chr5,length=180915260>
		##contig=<ID=chr6,length=171115067>
		##contig=<ID=chr7,length=159138663>
		##contig=<ID=chr8,length=146364022>
		##contig=<ID=chr9,length=141213431>
		##contig=<ID=chr10,length=135534747>
		##contig=<ID=chr11,length=135006516>
		##contig=<ID=chr12,length=133851895>
		##contig=<ID=chr13,length=115169878>
		##contig=<ID=chr14,length=107349540>
		##contig=<ID=chr15,length=102531392>
		##contig=<ID=chr16,length=90354753>
		##contig=<ID=chr17,length=81195210>
		##contig=<ID=chr18,length=78077248>
		##contig=<ID=chr19,length=59128983>
		##contig=<ID=chr20,length=63025520>
		##contig=<ID=chr21,length=48129895>
		##contig=<ID=chr22,length=51304566>
		##contig=<ID=chrM,length=16571>
		##contig=<ID=chrX,length=155270560>
		##contig=<ID=chrY,length=59373566>
		##source=HaplotypeCaller
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA19240chr2122
	}
}

proc diffanalysisinfo {file1 file2} {
	set file1 [lindex [glob $file1] 0]
	set file2 [lindex [glob $file2] 0]
	catch {exec cg tsvdiff -t xl $file1 $file2 | grep -v _version} temp
	set len [llength [split $temp \n]]
	if {$len != 3 && $len != 1} {return "files differ: $file1 $file2"}
	return ""
}

proc diffhtmlreport {file1 file2 {error 0}} {
	set file1 [lindex [glob $file1] 0]
	set file2 [lindex [glob $file2] 0]
	if {[catch {
		exec diff -r -y --suppress-common-lines $file1 $file2 \
		| grep -v {The full "sequenced genome" region in} \
		| grep -v {This report was generated by} \
		| grep -v {Made by genomecomb} \
		| grep -v {diff -r -y --suppress-common-lines}
	} msg]} {
		if {![string match {child process exited abnormally} $msg]} {
			if {$error} {error $msg} else {return $msg}
		}
	}
	return ""
}

proc cramdiff {file1 file2} {
	set temp1 tmp/[file root [file tail $file1]].tsv
	set temp2 tmp/[file root [file tail $file2]].tsv
	exec cg sam2tsv $file1 | cg select -f {qname chromosome begin end duplicate} > $temp1
	exec cg sam2tsv $file2 | cg select -f {qname chromosome begin end duplicate} > $temp2
	cg tsvdiff $temp1 $temp2
}

# remove tmp if it is a unexisting link
if {![file exists tmp]} {catch {file delete tmp}}
file mkdir tmp
set dbopt {}
set dboptt {}
