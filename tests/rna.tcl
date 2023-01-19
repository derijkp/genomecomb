#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" "$@"

source tools.tcl

test multitranscript {basic} {
	file_write tmp/t1.tsv [string trim [deindent {
		chromosome	begin	end	name	gene	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	source	gene_id	transcript_id	exon_number	count-t1
		chr3	999	2000	transcript1	gene1	-			1	999,	2000,	source2	gene1	transcript1	0	1
		chr4	999	4000	transcript2-1	gene2	+			2	999,2999,	1500,4000,	source2	gene2	transcript2-1	0,1	2
		chr4	999	5000	transcript2	gene2	+			3	999,2999,4499,	1500,4000,5000,	source2	gene2	transcript2-2	0,1,2	3
	}]]\n
	file_write tmp/t2.tsv [string trim [deindent {
		chromosome	begin	end	name	gene	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	source	gene_id	transcript_id	exon_number	count-t2
		chr3	999	2000	transcript1	gene1	-			1	999,	2000,	source2	gene1	transcript1	0	4
		chr4	999	4000	transcript2-1	gene2	+			3	999,1999,3499,	1500,3000,4000,	source2	gene2	transcript2-2	0,1,2	5
		chr4	999	4000	transcript2	gene2	+			2	999,2999,	1500,4000,	source2	gene2	transcript2-1	0,1	6
		chr4	1999	4000	transcript2-2	gene2	+			2	1999,2999,	2500,4000,	source2	gene2	transcript2-1	0,1	7
	}]]\n
	file_write tmp/expected.tsv [string trim [deindent {
		chromosome	begin	end	strand	exonStarts	exonEnds	cdsStart	cdsEnd	transcript	gene	geneid	category	name	exonCount	source	exon_number	count-t1	count-t2
		chr3	999	2000	-	999	2000			transcript1	gene1	gene1		transcript1	1	source2	0	1	4
		chr4	999	4000	+	999,1999,3499	1500,3000,4000			transcript2-2	gene2	gene2		transcript2-1	3	source2	0,1,2	0.0	5
		chr4	999	4000	+	999,2999	1500,4000			transcript2-1	gene2	gene2		transcript2-1	2	source2	0,1	2	6
		chr4	999	5000	+	999,2999,4499	1500,4000,5000			transcript2-2	gene2	gene2		transcript2	3	source2	0,1,2	3	0.0
		chr4	1999	4000	+	1999,2999	2500,4000			transcript2-1	gene2	gene2		transcript2-2	2	source2	0,1	0.0	7
	}]]\n
	file delete tmp/test.tsv
	cg multitranscript -stack 1 -exact 1 tmp/test.tsv tmp/t1.tsv tmp/t2.tsv
	exec diff tmp/test.tsv tmp/expected.tsv
} {}

test multitranscript {single overlaps (mixed match)} {
	file_write tmp/t1.tsv [string trim [deindent {
		chromosome	begin	end	name	gene	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	source	gene_id	transcript_id	exon_number	count-t1
		chr3	999	2000	known1	gene1	-			1	999	2000	source2	gene1	known1	0	1
		chr3	1500	2500	transcript1-2	gene1	-			1	1500	2500	source2	gene1	transcript1-2	0	1.5
		chr4	999	4000	transcript2-1	gene2	+			2	999,2999	1500,4000	source2	gene2	transcript2-1	0,1	2
		chr4	999	5000	transcript2	gene2	+			3	999,2999,4499	1500,4000,5000	source2	gene2	transcript2-2	0,1,2	3
	}]]\n
	file_write tmp/t2.tsv [string trim [deindent {
		chromosome	begin	end	name	gene	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	source	gene_id	transcript_id	exon_number	count-t2
		chr3	999	1990	transcript1	gene1	-			1	999	2000	source2	gene1	transcript1	0	4
		chr3	2000	2500	transcript1-2	gene1	-			1	2000	2500	source2	gene1	transcript1-2	0	4.5
		chr4	999	4000	transcript2-1	gene2	+			3	999,1999,3499	1500,3000,4000	source2	gene2	transcript2-2	0,1,2	5
		chr4	999	4000	transcript2	gene2	+			2	999,2999	1500,4000	source2	gene2	transcript2-1	0,1	6
		chr4	1999	4000	known2-2	gene2	+			2	1999,2999	2500,4000	source2	gene2	known2-2	0,1	7
		chr4	3990	5000	transcript2-3	gene2	+			2	3990	5000	source2	gene2	transcript2-3	0,1	8
		chr4	5500	6000	transcript3	gene3	+			2	5500	6000	source2	gene3	transcript3	0,1	8
	}]]\n
	file_write tmp/expected.tsv [string trim [deindent {
		chromosome	begin	end	strand	exonStarts	exonEnds	cdsStart	cdsEnd	transcript	gene	geneid	category	name	exonCount	source	exon_number	count-t1	count-t2
		chr3	999	2000	-	999	2000			known1	gene1	gene1		known1	1	source2	0	1	4
		chr3	1500	2500	-	1500	2500			transcript1-2	gene1	gene1		transcript1-2	1	source2	0	1.5	4.5
		chr4	999	4000	+	999,1999,3499	1500,3000,4000			transcript2-2	gene2	gene2		transcript2-1	3	source2	0,1,2	0.0	5
		chr4	999	4000	+	999,2999	1500,4000			transcript2-1	gene2	gene2		transcript2-1	2	source2	0,1	2	6
		chr4	999	5000	+	999,2999,4499	1500,4000,5000			transcript2-2	gene2	gene2		transcript2	3	source2	0,1,2	3	0.0
		chr4	1999	4000	+	1999,2999	2500,4000			known2-2	gene2	gene2		known2-2	2	source2	0,1	0.0	7
		chr4	3990	5000	+	3990	5000			transcript2-3	gene2	gene2		transcript2-3	2	source2	0,1	0.0	8
		chr4	5500	6000	+	5500	6000			transcript3	gene3	gene3		transcript3	2	source2	0,1	0.0	8
	}]]\n
	file delete tmp/test.tsv
	cg multitranscript -match transcript tmp/test.tsv tmp/t1.tsv tmp/t2.tsv
	exec diff tmp/test.tsv tmp/expected.tsv
} {}

test multitranscript {single overlaps (mixed cat)} {
	file_write tmp/t1.tsv [string trim [deindent {
		chromosome	begin	end	category	name	gene	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	source	gene_id	transcript_id	exon_number	count-t1
		chr3	999	2000	known	known1	gene1	-			1	999	2000	source2	gene1	known1	0	1
		chr3	1500	2500	novel	transcript1-2	gene1	-			1	1500	2500	source2	gene1	transcript1-2	0	1.5
		chr4	999	4000	novel	transcript2-1	gene2	+			2	999,2999	1500,4000	source2	gene2	transcript2-1	0,1	2
		chr4	999	5000	novel	transcript2	gene2	+			3	999,2999,4499	1500,4000,5000	source2	gene2	transcript2-2	0,1,2	3
	}]]\n
	file_write tmp/t2.tsv [string trim [deindent {
		chromosome	begin	end	category	name	gene	strand	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	source	gene_id	transcript_id	exon_number	count-t2
		chr3	999	1990	novel	transcript1	gene1	-			1	999	2000	source2	gene1	transcript1	0	4
		chr3	2000	2500	novel	transcript1-2	gene1	-			1	1500	2500	source2	gene1	transcript1-2	0	4.5
		chr4	999	4000	novel	transcript2-1	gene2	+			3	999,1999,3499	1500,3000,4000	source2	gene2	transcript2-2	0,1,2	5
		chr4	999	4000	novel	transcript2	gene2	+			2	999,2999	1500,4000	source2	gene2	transcript2-1	0,1	6
		chr4	1999	4000	known	known2-2	gene2	+			2	1999,2999	2500,4000	source2	gene2	known2-2	0,1	7
		chr4	3990	5000	novel	transcript2-3	gene2	+			2	3990	5000	source2	gene2	transcript2-3	0,1	8
		chr4	5500	6000	novel	transcript3	gene3	+			2	5500	6000	source2	gene3	transcript3	0,1	8
	}]]\n
	file_write tmp/expected.tsv [string trim [deindent {
		chromosome	begin	end	strand	exonStarts	exonEnds	cdsStart	cdsEnd	transcript	gene	geneid	category	name	exonCount	source	exon_number	count-t1	count-t2
		chr3	999	2000	-	999	2000			known1	gene1	gene1	known	known1	1	source2	0	1	4
		chr3	1500	2500	-	1500	2500			novelt_chr3_1500-e1000	gene1	gene1	novel	transcript1-2	1	source2	0	1.5	4.5
		chr4	999	4000	+	999,1999,3499	1500,3000,4000			novelt_chr4_999+e501i499e1001i499e501	gene2	gene2	novel	transcript2-1	3	source2	0,1,2	0.0	5
		chr4	999	4000	+	999,2999	1500,4000			novelt_chr4_999+e501i1499e1001	gene2	gene2	novel	transcript2-1	2	source2	0,1	2	6
		chr4	999	5000	+	999,2999,4499	1500,4000,5000			novelt_chr4_999+e501i1499e1001i499e501	gene2	gene2	novel	transcript2	3	source2	0,1,2	3	0.0
		chr4	1999	4000	+	1999,2999	2500,4000			known2-2	gene2	gene2	known	known2-2	2	source2	0,1	0.0	7
		chr4	3990	5000	+	3990	5000			novelt_chr4_3990+e1010	gene2	gene2	novel	transcript2-3	2	source2	0,1	0.0	8
		chr4	5500	6000	+	5500	6000			novelt_chr4_5500+e500	gene3	gene3	novel	transcript3	2	source2	0,1	0.0	8
	}]]\n
	file delete tmp/test.tsv
	cg multitranscript tmp/test.tsv tmp/t1.tsv tmp/t2.tsv
	exec diff tmp/test.tsv tmp/expected.tsv
} {}

test multitranscript {approx match (for novel)} {
	file_write tmp/t1.tsv [string trim [deindent {
		chromosome	begin	end	name	gene	strand	exonCount	exonStarts	exonEnds	gene_id	transcript_id	count-t1
		chr3	999	2000	transcript1-1n	gene1	-	1	999,	2000,	gene1	transcript1-1	11
		chr4	999	4000	transcript1-2n	gene2	+	2	999,2999,	1500,4000,	gene2	transcript1-2	12
		chr4	999	5000	transcript1-3n	gene2	+	3	999,2999,4499,	1500,4000,5000,	gene2	transcript1-3	13
		chr4	4850	4900	transcript1-4n	gene3	-	2	4850,	4900,	gene3	transcript1-4	14
	}]]\n
	file_write tmp/t2.tsv [string trim [deindent {
		chromosome	begin	end	name	gene	strand	exonCount	exonStarts	exonEnds	source	gene_id	transcript_id	count-t2
		chr3	998	2001	transcript2-1n	gene1	-	1	998,	2001,	source2	gene1	transcript2-1	21
		chr4	990	4000	transcript2-2n	gene2	+	2	990,2999,	1500,4000,	source2	gene2	transcript2-2	22
		chr4	999	4000	transcript2-3n	gene2	-	3	999,1999,3499,	1500,3000,4000,	source2	gene2	transcript2-3	23
		chr4	1999	4000	transcript2-4n	gene2	-	2	1999,2999,	2500,4000,	source2	gene2	transcript2-4	24
		chr4	4500	4800	transcript2-5n	gene3	-	2	4500,	4800,	source3	gene3	transcript2-5	25
		chr4	5000	5200	transcript2-6n	gene4	-	2	5000,	5200,	source4	gene4	transcript2-6	26
	}]]\n
	file_write tmp/expected.tsv [string trim [deindent {
		chromosome	begin	end	strand	exonStarts	exonEnds	cdsStart	cdsEnd	transcript	gene	geneid	category	name	exonCount	count-t1	count-t2
		chr3	998	2001	-	998	2001			transcript2-1	gene1	gene1		transcript2-1n	1	11	21
		chr4	990	4000	+	990,2999	1500,4000			transcript1-2	gene2	gene2		transcript1-2n	2	12	22
		chr4	999	4000	-	999,1999,3499	1500,3000,4000			transcript2-3	gene2	gene2		transcript2-3n	3	0.0	23
		chr4	999	5000	+	999,2999,4499	1500,4000,5000			transcript1-3	gene2	gene2		transcript1-3n	3	13	0.0
		chr4	1999	4000	-	1999,2999	2500,4000			transcript2-4	gene2	gene2		transcript2-4n	2	0.0	24
		chr4	4500	4800	-	4500	4800			transcript2-5	gene3	gene3		transcript2-5n	2	0.0	25
		chr4	4850	4900	-	4850	4900			transcript1-4	gene3	gene3		transcript1-4n	2	14	0.0
		chr4	5000	5200	-	5000	5200			transcript2-6	gene4	gene4		transcript2-6n	2	0.0	26
	}]]\n
	file delete tmp/test.tsv
	cg multitranscript -match transcript tmp/test.tsv tmp/t1.tsv tmp/t2.tsv
	exec diff tmp/test.tsv tmp/expected.tsv
} {}

test multicount {basic} {
	file_write tmp/gene_count-t1.tsv [string trim [deindent {
		gene	count-t1
		g1	11
		g2	12
		g3	13
		g4	14
	}]]\n
	file_write tmp/gene_count-t2.tsv [string trim [deindent {
		gene	count-t2
		g1	21
		g2	22
		g3	23
		g4	24
	}]]\n
	file_write tmp/gene_count-expected.tsv [string trim [deindent {
		gene	count-t1	count-t2
		g1	11	21
		g2	12	22
		g3	13	23
		g4	14	24
	}]]\n
	file delete tmp/gene_count-test.tsv
	cg multicount tmp/gene_count-test.tsv tmp/gene_count-t1.tsv tmp/gene_count-t2.tsv
	exec diff tmp/gene_count-test.tsv tmp/gene_count-expected.tsv
} {}

testsummarize
