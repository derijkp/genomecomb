#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

# todo:
# genotype in haploid calls (Y chromosome)

proc cg_gtf2sft {args} {
	cg_gtf2tsv {*}$args
}

proc cg_gtf2tsv {args} {
	set filename -
	set outfile -
	set separate 0
	set ignorecodon 1
	cg_options gtf2tsv args {
		-separate {
			set separate $value
		}
		-ignorecodon {
			set ignorecodon $value
		}
	} {filename outfile} 0 2

	catch {close $f} ;	catch {close $fb} ; catch {close $fa}
	set f [gzopen $filename]
	set comment {# -- tsv converted from gtf, original comments follow --}
	while {![eof $f]} {
		set line [gets $f]
		if {[string index $line 0] ne "#"} break
		append comment \n$line
	}
	append comment "\n# ----"
	set next 100000; set num 0
	set tempbase [tempfile]
	set tempattr [tempfile]
	set fb [open $tempbase w]
	set fa [open $tempattr w]
	set attrheader {}
	set attrtemplate {}
	unset -nocomplain attra
	if {!$separate} {
		set cheader [deindent {
			#filetype	tsv/transcriptsfile
			#fileversion	0.99
			#fields	table
			#fields	field	number	type	description
			#fields	chromosome	1	String	Chromosome name
			#fields	begin	1	Integer	Transcription start position
			#fields	end	1	Integer	Transcription end position
			#fields	transcript	1	String	Name of transcript (usually transcript_id from GTF)
			#fields	gene	1	String	Alternate name / name of gene (e.g. gene_id from GTF)
			#fields	strand	1	String	+ or - for strand
			#fields	cdsStart	1	Integer	Coding region start
			#fields	cdsEnd	1	Integer	Coding region end
			#fields	exonCount	1	Integer	Number of exons
			#fields	exonStarts	E	Integer	Exon start positions
			#fields	exonEnds	E	Integer	Exon end positions
			#fields	source	1	String	Source of data
		}]
		set nheader {chromosome begin end transcript gene strand cdsStart cdsEnd exonCount exonStarts exonEnds source}
		unset -nocomplain curchromosome
		set curtranscript {}
		set curlines {}
		set num 0
		set line [split $line \t]
		while 1 {
			while 1 {
				if {[eof $f] && $line eq ""} break
				foreach {chrom source type start end score strand phase attributes comments} $line break
				if {$type in "exon CDS start_codon stop_codon"} break
				if {[gets $f line] == -1} break
				set line [split $line \t]
			}
			if {[eof $f] && $line eq ""} break
			set a [dict create {*}[string_change $attributes {; " "}]]
			set curtranscript [dict get $a transcript_id]
			set curlines [list $line]
			while {[gets $f line] != -1} {
				if {[string index $line 0] eq "#"} continue
				set line [split $line \t]
				if {![llength $line]} continue
				foreach {chrom source type start end score strand phase attributes comments} $line break
				if {![inlist {exon CDS start_codon stop_codon} $type]} {
					continue
				}
				set a [dict create {*}[string_change $attributes {; " "}]]
				set transcript [dict get $a transcript_id]
				if {$transcript ne $curtranscript} break
				lappend curlines $line
			}
			if {![llength $curlines]} break
			# join $curlines \n

			if {$num >= $next} {putsprogress $curchromosome:$curbegin-$curend; incr next 100000}
			incr num

			# 
			set curlines [lsort -dict -index 3 $curlines]
			# start next
			foreach {curchromosome source type curbegin curend score curstrand phase attributes comments} [lindex $curlines 0] break
			unset -nocomplain curattra
			set cursource {}
			lappend cursource $source
			set curcdsStart {}
			set curcdsEnd {}
			set curexonCount 0
			set curexonStarts {}
			set curexonEnds {}
			set curscore {}
			foreach iline $curlines {
				foreach {chrom source type start end score strand phase attributes comments} $iline break
				foreach {key value}  [string_change $attributes {; " "}] {
					lappend curattra($key) $value
				}
				incr start -1
				switch $type {
					exon {
						lappend curexonStarts $start
						lappend curexonEnds $end
						lappend curscore $score
						if {$start < $curbegin} {set curbegin $start}
						if {$end > $curend} {set curend $end}
						incr curexonCount
					}
					CDS {
						if {$strand eq "+"} {
							if {$curcdsStart eq "" || $start < $curcdsStart} {set curcdsStart $start}
							set temp [expr {$end+3}]
							if {$curcdsEnd eq "" || $temp > $curcdsEnd} {set curcdsEnd $temp}
						} else {
							set temp [expr {$start-3}]
							if {$curcdsStart eq "" || $temp < $curcdsStart} {set curcdsStart $temp}
							if {$curcdsEnd eq "" || $end > $curcdsEnd} {set curcdsEnd $end}
						}
						if {$start < $curbegin} {set curbegin $start}
						if {$end > $curend} {set curend $end}
					}
					start_codon {
						if {!$ignorecodon} {
							if {$strand eq "+"} {
								set curcdsStart $start
							} else {
								set curcdsEnd $end
							}
						}
					}
					stop_codon {
						if {!$ignorecodon} {
							if {$strand eq "+"} {
								set curcdsEnd $end
							} else {
								set curcdsStart $start
							}
						}
					}
				}
			}
			if {![info exists curattra(gene_id)]} {
				error "field gene_id not found in attributes"
			}
			set curgene [list_remdup $curattra(gene_id)]
			if {![info exists curattra(transcript_id)]} {
				error "field transcript_id not found in attributes"
			}
			set curtranscript [list_remdup $curattra(transcript_id)]
			puts $fb [join [list $chrom $curbegin $curend $curtranscript $curgene $curstrand \
				$curcdsStart $curcdsEnd $curexonCount \
				[join $curexonStarts ,], [join $curexonEnds ,], [join [list_remdup $cursource] ,]] \t]
			set attrlist $attrtemplate
			foreach {key value} [array get curattra] {
				if {![info exists attra($key)]} {
					set attra($key) [llength $attrtemplate]
					lappend attrheader $key
					lappend attrtemplate {}
					lappend attrlist [join [list_remdup $value] ,]
				} else {
					lset attrlist $attra($key) [join [list_remdup $value] ,]
				}
			}
			puts $fa $attrlist
		}
	} else {
		set cheader [deindent {
			#filetype	tsv/transcriptsfile
			#fileversion	0.99
			#fields	table
			#fields	field	number	type	description
			#fields	chromosome	1	String	Chromosome name
			#fields	begin	1	Integer	Transcription start position
			#fields	end	1	Integer	Transcription end position
			#fields	type	1	String	type of element (transcript/exon)
			#fields	transcript	1	String	Name of transcript (usually transcript_id from GTF)
			#fields	gene	1	String	Alternate name / name of gene (e.g. gene_id from GTF)
			#fields	strand	1	String	+ or - for strand
			#fields	source	1	String	Source of data
			#fields	comments	1	String	extra info on element
		}]
		set nheader {chromosome begin end type transcript gene strand source comments}
		set num 0
		set curattr [dict create]
		while 1 {
			if {[string index $line 0] eq "#"} continue
			set line [split $line \t]
			if {![llength $line]} continue
			foreach {chrom source type start end score strand phase attributes comments} $line break
			incr start -1
			unset -nocomplain curattra
			foreach {key value}  [string_change $attributes {; " "}] {
				lappend curattra($key) $value
			}
			if {![info exists curattra(gene_id)]} {
				error "field gene_id not found in attributes"
			}
			set gene [list_remdup $curattra(gene_id)]
			if {![info exists curattra(transcript_id)]} {
				error "field transcript_id not found in attributes"
			}
			set transcript [list_remdup $curattra(transcript_id)]
			puts $fb [join [list $chrom $start $end $type $transcript $gene $strand $source $comments] \t]
			set attrlist $attrtemplate
			foreach {key value} [array get curattra] {
				if {![info exists attra($key)]} {
					set attra($key) [llength $attrtemplate]
					lappend attrheader $key
					lappend attrtemplate {}
					lappend attrlist [join [list_remdup $value] ,]
				} else {
					lset attrlist $attra($key) [join [list_remdup $value] ,]
				}
			}
			puts $fa $attrlist
			if {[gets $f line] == -1} break
		}
	}
	catch {close $fb} ; catch {close $fa}
	if {$f ne "stdin"} {catch {gzclose $f}}

	putsprogress "Assembling file"
	set o [wgzopen $outfile]
	puts $o $cheader
	puts $o $comment
	puts $o [join [list_concat $nheader $attrheader] \t]
	if {[info exists attrlist]} {
		set len [llength $attrlist]
		set fb [open $tempbase]
		set fa [open $tempattr]
		set len [llength $attrlist]
		while {![eof $fa]} {
			set linea [gets $fb]
			if {![string length $linea]} break
			set lineb [gets $fa]
			if {[llength $lineb] < $len} {
				lappend lineb {*}[list_fill [expr {$len-[llength $lineb]}] {}]
			}
			puts $o $linea\t[join $lineb \t]
		}
	}
	catch {close $fb} ; catch {close $fa}
	file delete $tempbase ; file delete $tempattr
	if {$o ne "stdout"} {catch {close $o}}
}
