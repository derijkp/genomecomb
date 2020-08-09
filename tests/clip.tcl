#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

if 0 {
	exec samtools view -hSb tmp/temp.sam > tmp/temp.bam
	exec samtools index tmp/temp.bam
	exec igv tmp/temp.bam &
}

proc checksam {file testfile {refseq {}}} {
	set ext [file extension $file]
	set error [catch {exec samtools index $file} msg]
	if {$ext eq ".sam" && (!$error 
		|| (
			![regexp "is in a format that cannot be usefully indexed" $msg]
			&& ![regexp "SAM file .* not BGZF compressed" $msg]
		)
	)} {
		error "$file is not a sam file"
	} elseif {$ext eq ".bam" && ![file exists $file.bai]} {
		error "$file is not a bam file"
	} elseif {$ext eq ".cram" && ![file exists $file.crai]} {
		error "$file is not a cram file"
	}
	if {$ext eq ".bam"} {
		set tempfile [tempfile]
		exec samtools view --no-PG -h $file > $tempfile
		set file $tempfile
	} elseif {[file extension $file] eq ".cram"} {
		set refseq [refseq $refseq]
		set tempfile [tempfile]
		exec samtools view --no-PG -h -T $refseq $file > $tempfile
		set file $tempfile
	}
	cg select -sh /dev/null -f {qname flag rname pos seq qual} $file > tmp/result.tsv
	exec diff tmp/result.tsv $testfile
}

test fastq_clipadapters {fastq_clipadapters single} {
	file copy -force {*}[glob data/seq_R1.fq.gz] tmp/
	cg fastq_clipadapters tmp/seq_R1.fq.gz
	# nothing is actually clipped here, so there should be no difference
	cg tsvdiff tmp/seq_R1.fq.gz tmp/fastq.clipped/seq_R1.clipped.fastq.gz
} {}

test fastq_clipadapters {fastq_clipadapters paired} {
	file copy -force {*}[glob data/seq_R*.fq.gz] tmp/
	cg fastq_clipadapters -paired 1 tmp/seq_R1.fq.gz tmp/seq_R2.fq.gz
	# nothing is actually clipped here, so there should be no difference
	cg tsvdiff tmp/seq_R1.fq.gz tmp/fastq.clipped/seq_R1.clipped.fastq.gz
	cg tsvdiff tmp/seq_R2.fq.gz tmp/fastq.clipped/seq_R2.clipped.fastq.gz
} {}

test fastq_clipadapters {fastq_clipadapters paired, multiple fastqs} {
	cg fastq_split -parts 4 data/seq_R1.fq.gz tmp/seq_R1.fq.gz
	cg fastq_split -parts 4 data/seq_R2.fq.gz tmp/seq_R2.fq.gz
	cg fastq_clipadapters -paired 1 {*}[glob tmp/*.fq.gz]
	# nothing is actually clipped here, so there should be no difference
	foreach file [glob tmp/*.fq.gz] {
		cg tsvdiff $file tmp/fastq.clipped/[file root [file tail [gzroot $file]]].clipped.fastq.gz
	}
} {}

test sam_clipamplicons {basic} {
	write_sam tmp/temp.sam {
		chr1	100	20M	20	chr1	121	20M	20
		chr1	100	2M1X17M	20	chr1	121	17M1X2M	20
		chr1	100	2S18M	20	chr1	121	2S16M2S	20
		chr1	100	2S2M2D16M	20	chr1	120	2S16M2D2M	20
		chr1	100	2S2M2I14M	20	chr1	123	2S14M2I2M	20
		chr1	102	4S2M2I12M	20	chr1	111	30M	30
		chr1	107	20M	20	chr1	112	20M	20
		chr1	108	20M	20	chr1	113	20M	20
		chr2	50	20M	20	chr2	60	20M	20
		chr2	100	50M	50	chr2	100	50M	50
	}
	write_tab tmp/samplicons.tsv {
		chromosome outer_begin begin end outer_end
		chr1 99 107 131 140
		chr2 99 107 141 149
	}
	write_tab tmp/expected.tsv {
		A1	99	chr1	100	NNNNNNNNAAAAAAAAAAAA	!!!!!!!!------------
		A2	99	chr1	100	NNNNNNNNAAAAAAAAAAAA	!!!!!!!!------------
		A3	99	chr1	100	NNNNNNNNNNAAAAAAAAAA	!!!!!!!!!!----------
		A4	99	chr1	100	NNNNNNNNAAAAAAAAAAAA	!!!!!!!!------------
		A5	99	chr1	100	NNNNNNNNNNNNAAAAAAAA	!!!!!!!!!!!!--------
		A6	99	chr1	102	NNNNNNNNNNNNAAAAAAAA	!!!!!!!!!!!!--------
		A7	99	chr1	107	NAAAAAAAAAAAAAAAAAAA	!-------------------
		A8	99	chr1	108	AAAAAAAAAAAAAAAAAAAA	--------------------
		A6	147	chr1	111	AAAAAAAAAAAAAAAAAAAAANNNNNNNNN	---------------------!!!!!!!!!
		A7	147	chr1	112	AAAAAAAAAAAAAAAAAAAA	--------------------
		A8	147	chr1	113	AAAAAAAAAAAAAAAAAAAN	-------------------!
		A4	147	chr1	120	AAAAAAAAAAAAAANNNNNN	--------------!!!!!!
		A1	147	chr1	121	AAAAAAAAAAANNNNNNNNN	-----------!!!!!!!!!
		A2	147	chr1	121	AAAAAAAAAAANNNNNNNNN	-----------!!!!!!!!!
		A3	147	chr1	121	AAAAAAAAAAAAANNNNNNN	-------------!!!!!!!
		A5	147	chr1	123	AAAAAAAAAAANNNNNNNNN	-----------!!!!!!!!!
		A9	99	chr2	50	AAAAAAAAAAAAAAAAAAAA	--------------------
		A9	147	chr2	60	AAAAAAAAAAAAAAAAAAAA	--------------------
		A10	99	chr2	100	NNNNNNNNAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNN	!!!!!!!!----------------------------------!!!!!!!!
		A10	147	chr2	100	NNNNNNNNAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNN	!!!!!!!!----------------------------------!!!!!!!!
	}
	cg sam_clipamplicons tmp/samplicons.tsv tmp/temp.sam tmp/out.sam
	cg select -sh /dev/null -f {qname flag rname pos seq qual} tmp/out.sam > tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test sam_clipamplicons {basic 2} {
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
		chr1 99 109 129 139
		chr1 119 129 149 159
		chr1 129 139 159 169
	}
	write_tab tmp/expected.tsv {
		A1	99	chr1	85	NNNNNNNNNNNNNNNNNNNN	!!!!!!!!!!!!!!!!!!!!
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
		A9	147	chr1	200	GGGGGGGGGGGGGGGGGGGG	--------------------
	}
	cg sam_clipamplicons tmp/samplicons.tsv tmp/temp.sam tmp/out.sam
	cg select -sh /dev/null -f {qname flag rname pos seq qual} tmp/out.sam > tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test sam_clipamplicons {pipes and formats} {
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
		chr1 99 109 129 139
		chr1 119 129 149 159
		chr1 129 139 159 169
	}
	write_tab tmp/expected.tsv {
		A1	99	chr1	85	NNNNNNNNNNNNNNNNNNNN	!!!!!!!!!!!!!!!!!!!!
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
		A9	147	chr1	200	GGGGGGGGGGGGGGGGGGGG	--------------------
	}
	cg sam_clipamplicons tmp/samplicons.tsv < tmp/temp.sam > tmp/out.sam
	checksam tmp/out.sam tmp/expected.tsv
	# bam input and output
	exec samtools view --no-PG -b tmp/temp.sam > tmp/temp.bam
	cg sam_clipamplicons tmp/samplicons.tsv tmp/temp.bam tmp/out.bam
	checksam tmp/out.bam tmp/expected.tsv
	cg sam_clipamplicons -outputformat bam tmp/samplicons.tsv tmp/temp.bam > tmp/out.bam
	checksam tmp/out.bam tmp/expected.tsv
	# cram in and output via pipe
	exec samtools view --no-PG -h -C -T [refseq $::refseqdir/hg19] tmp/temp.sam > tmp/temp.cram
	cg sam_clipamplicons -stack 1 -refseq $::refseqdir/hg19 tmp/samplicons.tsv tmp/temp.cram tmp/outd.cram
	checksam tmp/outd.cram tmp/expected.tsv $::refseqdir/hg19
	cg sam_clipamplicons -stack 1 -refseq $::refseqdir/hg19 -inputformat cram -outputformat cram tmp/samplicons.tsv < tmp/temp.cram > tmp/out.cram
	checksam tmp/out.cram tmp/expected.tsv $::refseqdir/hg19
	cg sam_clipamplicons -refseq $::refseqdir/hg19 tmp/samplicons.tsv tmp/temp.sam tmp/outsd.cram
	checksam tmp/outsd.cram tmp/expected.tsv $::refseqdir/hg19
} {}

test sam_clipamplicons {skip chromosome in amplicons} {
	write_sam tmp/temp.sam {
		chr2	50	20M	20	chr2	60	20M	20
		chr2	100	50M	50	chr2	100	50M	50
	}
	write_tab tmp/samplicons.tsv {
		chromosome outer_begin begin end outer_end
		chr1 99 107 131 140
		chr2 99 107 141 149
	}
	write_tab tmp/expected.tsv {
		A1	99	chr2	50	AAAAAAAAAAAAAAAAAAAA	--------------------
		A1	147	chr2	60	AAAAAAAAAAAAAAAAAAAA	--------------------
		A2	99	chr2	100	NNNNNNNNAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNN	!!!!!!!!----------------------------------!!!!!!!!
		A2	147	chr2	100	NNNNNNNNAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNN	!!!!!!!!----------------------------------!!!!!!!!
	}
	cg sam_clipamplicons tmp/samplicons.tsv tmp/temp.sam tmp/out.sam
	cg select -sh /dev/null -f {qname flag rname pos seq qual} tmp/out.sam > tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test sam_clipamplicons {skip chromosome in sam} {
	write_sam tmp/temp.sam {
		chr1	100	20M	20	chr1	121	20M	20
		chr2	50	20M	20	chr2	60	20M	20
		chr2	100	50M	50	chr2	100	50M	50
	}
	write_tab tmp/samplicons.tsv {
		chromosome outer_begin begin end outer_end
		chr2 99 107 141 149
	}
	write_tab tmp/expected.tsv {
		A1	99	chr1	100	AAAAAAAAAAAAAAAAAAAA	--------------------
		A1	147	chr1	121	AAAAAAAAAAAAAAAAAAAA	--------------------
		A2	99	chr2	50	AAAAAAAAAAAAAAAAAAAA	--------------------
		A2	147	chr2	60	AAAAAAAAAAAAAAAAAAAA	--------------------
		A3	99	chr2	100	NNNNNNNNAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNN	!!!!!!!!----------------------------------!!!!!!!!
		A3	147	chr2	100	NNNNNNNNAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNN	!!!!!!!!----------------------------------!!!!!!!!
	}
	cg sam_clipamplicons tmp/samplicons.tsv tmp/temp.sam tmp/out.sam
	cg select -sh /dev/null -f {qname flag rname pos seq qual} tmp/out.sam > tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test sam_clipamplicons {completely in primer} {
	write_sam tmp/temp.sam {
		chr2	100	8M	8	chr2	142	8M	8
	}
	write_tab tmp/samplicons.tsv {
		chromosome outer_begin begin end outer_end
		chr2 99 107 141 149
	}
	write_tab tmp/expected.tsv {
		A1	99	chr2	100	NNNNNNNN	!!!!!!!!
		A1	147	chr2	142	NNNNNNNN	!!!!!!!!
	}
	cg sam_clipamplicons tmp/samplicons.tsv tmp/temp.sam tmp/out.sam
	cg select -sh /dev/null -f {qname flag rname pos seq qual} tmp/out.sam > tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test sam_clipamplicons {completely in primer wrong end} {
	write_sam tmp/temp.sam {
		chr2	100	8M	8	chr2	100	8M	8
		chr2	142	8M	8	chr2	142	8M	8
	}
	write_tab tmp/samplicons.tsv {
		chromosome outer_begin begin end outer_end
		chr2 99 107 141 149
	}
	write_tab tmp/expected.tsv {
		A1	99	chr2	100	NNNNNNNN	!!!!!!!!
		A1	147	chr2	100	NNNNNNNN	!!!!!!!!
		A2	99	chr2	142	NNNNNNNN	!!!!!!!!
		A2	147	chr2	142	NNNNNNNN	!!!!!!!!
	}
	cg sam_clipamplicons tmp/samplicons.tsv tmp/temp.sam tmp/out.sam
	cg select -sh /dev/null -f {qname flag rname pos seq qual} tmp/out.sam > tmp/result.tsv
	exec diff tmp/result.tsv tmp/expected.tsv
} {}

test sam_ampliconscount {basic} {
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
		chr1 99 109 129 139
		chr1 119 129 149 159
		chr1 129 139 159 169
	}
	write_tab tmp/expected.tsv {
		chromosome	begin	end	name	count
		chr1	109	129	{}	10
		chr1	129	149	{}	4
		chr1	139	159	{}	3
	}
	cg sam_ampliconscount tmp/samplicons.tsv tmp/temp.sam tmp/count.tsv
	exec diff tmp/count.tsv tmp/expected.tsv
} {}

testsummarize
