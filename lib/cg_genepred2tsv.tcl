#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_genepred2tsv {args} {
	set filename -
	set outfile -
	set separate 0
	cg_options genepred2tsv args {
		-separate {
			set separate $value
		}
	} {filename outfile} 0 2
	set f [gzopen $filename]
	set o [wgzopen $outfile]
	puts $o [deindent {
		#filetype	tsv/transcriptsfile
		#fileversion	0.99
		#fields	table
		#fields	field	number	type	description
		#fields	name	1	String	Name of transcript (usually transcript_id from GTF)
		#fields	chromosome	1	String	Chromosome name
		#fields	strand	1	String	+ or - for strand
		#fields	begin	1	Integer	Transcription start position
		#fields	end	1	Integer	Transcription end position
		#fields	cdsStart	1	Integer	Coding region start
		#fields	cdsEnd	1	Integer	Coding region end
		#fields	exonCount	1	Integer	Number of exons
		#fields	exonStarts	E	Integer	Exon start positions
		#fields	exonEnds	E	Integer	Exon end positions
		#fields	score	1	Float	Score
		#fields	gene	1	String	Alternate name / name of gene (e.g. gene_id from GTF)
		#fields	cdsStartStat	1	String	Status of CDS start annotation (none, unknown, incomplete, or complete)
		#fields	cdsEndStat	1	String	Status of CDS end annotation (none, unknown, incomplete, or complete)
		#fields	exonFrames	E	Integer	Exon frame offsets {0,1,2}
	}]\n[join {
		name chromosome strand begin end cdsStart cdsEnd 
		exonCount exonStarts exonEnds score gene cdsStartStat cdsEndStat exonFrames
	} \t]
	while {[gets $f line] != -1} {
		set line [split $line \t]
		lset line 8 [string trimright [lindex $line 8] ,]
		lset line 9 [string trimright [lindex $line 9] ,]
		puts $o [join $line \t]
	}
	if {$f ne "stdin"} {catch {gzclose $f}}
	if {$o ne "stdout"} {catch {close $o}}
}
