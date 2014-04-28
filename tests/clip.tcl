#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

proc write_sam {file data} {
	set tempfile [tempfile]
	set o [open $tempfile w]
	set num 1
	foreach line [split [string trim $data] \n] {
		set base A
		foreach {chr1 pos1 cigar1 size1 chr2 pos2 cigar2 size2 base} $line break
		if {$base eq ""} {set base A}
		set seq [string_fill $base $size1]
		set qual [string_fill - $size1]
		if {$chr2 eq $chr1} {set c2 =} else {set c2 $chr2}
		puts $o [join [list A$num 99 $chr1 $pos1 60 $cigar1 $c2 $pos2 [expr {$pos2+$size2-$pos1}] $seq $qual RG:Z:sample1	NM:i:4	MQ:i:60	AS:i:241	XS:i:25] \t]
		set seq [string_fill $base $size2]
		set qual [string_fill - $size2]
		if {$chr2 eq $chr1} {set c1 =} else {set c1 $chr1}
		puts $o [join [list A$num 147 $chr2 $pos2 60 $cigar2 $c1 $pos1 -[expr {$pos2+$size2-$pos1}] $seq $qual RG:Z:sample1	NM:i:4	MQ:i:60	AS:i:241	XS:i:25] \t]
		incr num
	}
	close $o
	set o [open $file w]
	foreach line [split [string trim {
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
		@PG	ID:GATK IndelRealigner	VN:2.4-9-g532efad	{CL:knownAlleles=[] targetIntervals=test.intervals LODThresholdForCleaning=5.0 consensusDeterminationModel=USE_READS entropyThreshold=0.15 maxReadsInMemory=150000 maxIsizeForMovement=3000 maxPositionalMoveAllowed=200 maxConsensuses=30 maxReadsForConsensuses=120 maxReadsForRealignment=20000 noOriginalAlignmentTags=false nWayOut=null generate_nWayOut_md5s=false check_early=false noPGTag=false keepPGTags=false indelsFileForDebugging=null statisticsFileForDebugging=null SNPsFileForDebugging=null}
	}] \n] {
		puts $o [join $line \t]
	}
	close $o
	exec gnusort8 -t \t -N -s -k3,3 -k4,4 -k1,1 $tempfile >> $file
}

if 0 {
	exec samtools view -hSb tmp/temp.sam > tmp/temp.bam
	exec samtools index tmp/temp.bam
	exec igv tmp/temp.bam &
}

test sam_clipamplicons {basic} {
	write_sam tmp/temp.sam {
		chr1	100	20M	20	chr1	121	20M	20
		chr1	100	2M1X17M	20	chr1	121	17M1X2M	20
		chr1	100	2S18M	20	chr1	121	2S16M2S	20
		chr1	100	2S2M2D16M	20	chr1	120	2S16M2D2M	20
		chr1	100	2S2M2I14M	20	chr1	123	2S14M2I2M	20
		chr1	102	4S2M2I12M	20	chr1	111	30M	30
		chr2	50	20M	20	chr2	60	20M	20
		chr2	100	50M	50	chr2	100	50M	50
	}
	write_tab tmp/samplicons.tsv {
		chromosome outer_begin begin end outer_end
		chr1 100 108 132 141
		chr2 100 108 142 150
	}
	cg sam_clipamplicons tmp/samplicons.tsv tmp/temp.sam tmp/out.sam
	cg select -sh /dev/null -f {qname flag rname pos seq qual} tmp/out.sam > tmp/result.tsv
	write_tab tmp/expected.tsv {
		A1	99	chr1	100	NNNNNNNNAAAAAAAAAAAA	!!!!!!!!------------
		A2	99	chr1	100	NNNNNNNNAAAAAAAAAAAA	!!!!!!!!------------
		A3	99	chr1	100	NNNNNNNNNNAAAAAAAAAA	!!!!!!!!!!----------
		A4	99	chr1	100	NNNNNNNNAAAAAAAAAAAA	!!!!!!!!------------
		A5	99	chr1	100	NNNNNNNNNNNNAAAAAAAA	!!!!!!!!!!!!--------
		A6	99	chr1	102	NNNNNNNNNNNNAAAAAAAA	!!!!!!!!!!!!--------
		A6	147	chr1	111	AAAAAAAAAAAAAAAAAAAAANNNNNNNNN	---------------------!!!!!!!!!
		A4	147	chr1	120	AAAAAAAAAAAAAANNNNNN	--------------!!!!!!
		A1	147	chr1	121	AAAAAAAAAAANNNNNNNNN	-----------!!!!!!!!!
		A2	147	chr1	121	AAAAAAAAAAANNNNNNNNN	-----------!!!!!!!!!
		A3	147	chr1	121	AAAAAAAAAAAAANNNNNNN	-------------!!!!!!!
		A5	147	chr1	123	AAAAAAAAAAANNNNNNNNN	-----------!!!!!!!!!
		A7	99	chr2	50	AAAAAAAAAAAAAAAAAAAA	--------------------
		A7	147	chr2	60	AAAAAAAAAAAAAAAAAAAA	--------------------
		A8	99	chr2	100	NNNNNNNNAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNN	!!!!!!!!----------------------------------!!!!!!!!
		A8	147	chr2	100	NNNNNNNNAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNN	!!!!!!!!----------------------------------!!!!!!!!
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test sam_clipamplicons {-f} {
	write_sam tmp/temp.sam {
		chr1	85	20M	20	chr1	100	20M	20	W
		chr1	99	20M	20	chr1	120	15M	15	A
		chr1	110	25M	25	chr1	122	20M	20	A
		chr1	100	20M	20	chr1	120	20M	20	A
		chr1	101	20M	20	chr1	119	20M	20	A
		chr1	120	20M	20	chr1	140	20M	20	C
		chr1	120	40M	40	chr1	125	35M	35	C
		chr1	130	20M	20	chr1	150	20M	20	T
		chr1	160	20M	20	chr1	200	20M	20	G
	}
	write_tab tmp/samplicons.tsv {
		chromosome outer_begin begin end outer_end
		chr1 100 110 130 140
		chr1 120 130 150 160
		chr1 130 140 160 170
	}
	cg sam_clipamplicons tmp/samplicons.tsv tmp/temp.sam tmp/out.sam
	cg select -sh /dev/null -f {qname flag rname pos seq qual} tmp/out.sam > tmp/result.tsv
	write_tab tmp/expected.tsv {
	}
	exec diff tmp/result.tsv tmp/expected.tsv
} {A1	99	chr1	85	NNNNNNNNNNNNNNNNNNNN	!!!!!!!!!!!!!!!!!!!!
A2	99	chr1	99	NNNNNNNNNNNAAAAAAAAA	!!!!!!!!!!!---------
A1	147	chr1	100	NNNNNNNNNNWWWWWWWWWW	!!!!!!!!!!----------
A4	99	chr1	100	NNNNNNNNNNAAAAAAAAAA	!!!!!!!!!!----------
A5	99	chr1	101	NNNNNNNNNAAAAAAAAAAA	!!!!!!!!!-----------
A3	99	chr1	110	AAAAAAAAAAAAAAAAAAAANNNNN	--------------------!!!!!
A5	147	chr1	119	AAAAAAAAAAANNNNNNNNN	-----------!!!!!!!!!
A2	147	chr1	120	AAAAAAAAAANNNNN	----------!!!!!
A4	147	chr1	120	AAAAAAAAAANNNNNNNNNN	----------!!!!!!!!!!
A6	99	chr1	120	NNNNNNNNNNCCCCCCCCCC	!!!!!!!!!!----------
A7	99	chr1	120	NNNNNNNNNNCCCCCCCCCCCCCCCCCCCCNNNNNNNNNN	!!!!!!!!!!--------------------!!!!!!!!!!
A3	147	chr1	122	AAAAAAAANNNNNNNNNNNN	--------!!!!!!!!!!!!
A7	147	chr1	125	NNNNNCCCCCCCCCCCCCCCCCCCCNNNNNNNNNN	!!!!!--------------------!!!!!!!!!!
A8	99	chr1	130	NNNNNNNNNNTTTTTTTTTT	!!!!!!!!!!----------
A6	147	chr1	140	CCCCCCCCCCNNNNNNNNNN	----------!!!!!!!!!!
A8	147	chr1	150	TTTTTTTTTTNNNNNNNNNN	----------!!!!!!!!!!
A9	99	chr1	160	GGGGGGGGGGGGGGGGGGGG	--------------------
A9	147	chr1	200	GGGGGGGGGGGGGGGGGGGG	--------------------}

testsummarize
