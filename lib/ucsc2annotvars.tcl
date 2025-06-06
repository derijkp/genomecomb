#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc ucsc2annotvar {file outfile} {

	set f [gzopen $file]
	set header [gets $f]
	set header [split [string range $header 1 end] \t]
	if {$header ne "bin chrom chromStart chromEnd name alleleCount alleleFreq alleleScores"} {
		error "error in header"
	}
	set o [open $outfile w]
	puts $o [join {locus chromosome begin end type alleleSeq1 alleleSeq2 totalScore1 totalScore2 xRef geneId mrnaAcc proteinAcc orientation exonCategory exon codingRegionKnown aaCategory nucleotidePos proteinPos aaAnnot aaCall aaRef} \t]
	while {![eof $f]} {
		set line [getnotempty $f]
		foreach {bin chrom chromStart chromEnd name alleleCount alleleFreq alleleScores} $line break
		set alleles [split $name /]
		foreach {a1 a2} $alleles break
		if {[llength $alleles] == 1} {
			set a2 $a1
		}
		foreach {s1 s2} $alleleScores break
		if {[llength $alleleScores] == 1} {
			set s2 $s1
		}
		incr chromEnd -1
		set chrom [chr_clip $chrom]
		puts $o [join [list $bin $chrom $chromStart $chromEnd snp $a1 $a2] \t]
	}
	close $o
	gzclose $f

}


if 0 {

lappend auto_path $::env(HOME)/dev/comletegenomics/lib
package require Extral
set file /complgen/GS00028-DNA-C01/GS00028-ucsc2006-vars.tsv
set outfile /complgen/GS00028-DNA-C01/annot-GS00028-ucsc2006-vars.tsv
ucsc2annotvar $file $outfile

}
