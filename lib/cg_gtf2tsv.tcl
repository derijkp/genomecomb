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
	set transcripts 1
	set ignorecodon 1
	set sort 1
	cg_options gtf2tsv args {
		-separate {
			if {[true $value]} {set transcripts 0} else {set transcripts 1}
		}
		-transcripts {
			set transcripts $value
		}
		-ignorecodon {
			set ignorecodon $value
		}
		-sort {
			set sort $value
		}
	} {filename outfile} 0 2
	proc gtf2tsv_parse_attr {attributes} {
		set a [dict create {*}[string_change $attributes {; " "}]]
	}
	gtf2tsv $filename $outfile $transcripts $ignorecodon $sort
}

proc gtf2tsv {filename outfile {transcripts 1} {ignorecodon 1} {sort 1}} {
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
	unset -nocomplain genea
	if {$transcripts} {
		set cheader [deindent {
			#filetype	tsv/transcriptsfile
			#fileversion	0.99
			#fields	table
			#fields	field	number	type	description
			#fields	chromosome	1	String	Chromosome name
			#fields	begin	1	Integer	Transcription start position
			#fields	end	1	Integer	Transcription end position
			#fields	type	1	Integer	Type of element (transcript typically)
			#fields	transcript	1	String	Name of transcript (usually transcript_id from GTF)
			#fields	gene	1	String	Alternate name / name of gene (e.g. gene_name or gene_id from GTF)
			#fields	strand	1	String	+ or - for strand
			#fields	cdsStart	1	Integer	Coding region start
			#fields	cdsEnd	1	Integer	Coding region end
			#fields	exonCount	1	Integer	Number of exons
			#fields	exonStarts	E	Integer	Exon start positions
			#fields	exonEnds	E	Integer	Exon end positions
			#fields	source	1	String	Source of data
		}]
		set nheader {chromosome begin end type transcript gene strand cdsStart cdsEnd exonCount exonStarts exonEnds source}
		unset -nocomplain curchromosome
		set curtranscript {}
		set curlines {}
		set num 0
		set line [split $line \t]
		while 1 {
			while 1 {
				if {[eof $f] && $line eq ""} break
				foreach {chrom source type start end score strand phase attributes comments} $line break
				if {$type eq "transcript"} {
					set a [gtf2tsv_parse_attr $attributes]
					if {[dict exists $a Parent] && [dict exists $a transcript_id]} {
						set gene [dict get $a Parent]
						regsub ^gene: $gene {} gene
						set genea([dict get $a transcript_id]) $gene
					}
				} elseif {$type in "exon CDS start_codon stop_codon"} break
				if {[gets $f line] == -1} break
				set line [split $line \t]
			}
			if {[eof $f] && $line eq ""} break
			set a [gtf2tsv_parse_attr $attributes]
			if {[dict exists $a transcript_id]} {
				set curtranscript [dict get $a transcript_id]
			} else {
				set curtranscript [dict get $a Parent]
			}
			set curlines [list $line]
			while {[gets $f line] != -1} {
				if {[string index $line 0] eq "#"} continue
				set line [split $line \t]
				if {![llength $line]} continue
				foreach {chrom source type start end score strand phase attributes comments} $line break
				if {$type eq "transcript"} {
					set a [gtf2tsv_parse_attr $attributes]
					if {[dict exists $a Parent] && [dict exists $a transcript_id]} {
						set gene [dict get $a Parent]
						regsub ^gene: $gene {} gene
						set genea([dict get $a transcript_id]) $gene
					}
				} elseif {![inlist {exon CDS start_codon stop_codon} $type]} {
					continue
				}
				set a [gtf2tsv_parse_attr $attributes]
				if {[dict exists $a transcript_id]} {
					set transcript [dict get $a transcript_id]
				} else {
					set transcript [dict get $a Parent]
				}
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
				foreach {key value} [gtf2tsv_parse_attr $attributes] {
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
			if {[info exists curattra(transcript_id)]} {
				set curtranscript [list_remdup $curattra(transcript_id)]
			} elseif {[info exists curattra(Parent)]} {
				set curtranscript [list_remdup $curattra(Parent)]
			} else {
				error "field transcript_id not found in attributes"
			}
			regsub ^transcript: $curtranscript {} curtranscript
			if {[info exists curattra(gene_name)]} {
				set curgene [list_remdup $curattra(gene_name)]
			} elseif {[info exists curattra(gene_id)]} {
				set curgene [list_remdup $curattra(gene_id)]
			} elseif {[info exists genea($curtranscript)]} {
				set curgene $genea($curtranscript)
			} else {
				error "field gene_name or gene_id not found in attributes"
			}
			puts $fb [join [list $chrom $curbegin $curend transcript $curtranscript $curgene $curstrand \
				$curcdsStart $curcdsEnd $curexonCount \
				[join $curexonStarts ,] [join $curexonEnds ,] [join [list_remdup $cursource] ,]] \t]
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
		set nheader {chromosome begin end type strand transcript gene score source phase comments}
		set num 0
		set curattr [dict create]
		foreach key {chromosome type begin end score strand source phase transcript gene} {
			set transa($key) attr_$key
		}
		while 1 {
			if {[string index $line 0] eq "#"} continue
			set line [split $line \t]
			if {![llength $line]} continue
			foreach {chrom source type start end score strand phase attributes comments} $line break
			if {$type eq "transcript"} {
				set a [gtf2tsv_parse_attr $attributes]
				if {[dict exists $a Parent] && [dict exists $a transcript_id]} {
					set gene [dict get $a Parent]
					regsub ^gene: $gene {} gene
					set genea([dict get $a transcript_id]) $gene
				}
			}
			incr start -1
			unset -nocomplain curattra
			foreach {key value} [gtf2tsv_parse_attr $attributes] {
				lappend curattra($key) $value
			}
			if {[info exists curattra(transcript_id)]} {
				set transcript [list_remdup $curattra(transcript_id)]
				regsub ^transcript: $transcript {} transcript
			} else {
				set transcript .
			}
			if {[info exists curattra(gene_id)]} {
				set gene [list_remdup $curattra(gene_id)]
				regsub ^gene: $gene {} gene
			} elseif {[info exists genea($transcript)]} {
				set gene $genea($transcript)
			} else {
				set gene .
			}
			puts $fb [join [list $chrom $start $end $type $strand $transcript $gene $score $source $phase $comments] \t]
			set attrlist $attrtemplate
			foreach {key value} [array get curattra] {
				if {[info exists transa($key)]} {
					set key $transa($key)
				}
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
	if {$sort} {
		set tempout [tempfile].tsv.zst
		set o [wgzopen $tempout]
	} else {
		set o [wgzopen $outfile]
	}
	puts $o $cheader
	puts $o $comment
	puts $o [join [list_concat $nheader $attrheader] \t]
	set atrrsize [llength $attrheader]
	set fb [open $tempbase]
	set fa [open $tempattr]
	unset -nocomplain atemplate
	while {![eof $fa]} {
		if {[gets $fb lineb] == -1} break
		set linea [gets $fa]
		set len [llength $linea]
		if {$len < $atrrsize} {
			if {![info exists atemplate($len)]} {
				set atemplate($len) [string_fill \t [expr {$atrrsize-$len}]]
			}
			set linea [join $linea \t]$atemplate($len)
		} else {
			set linea [join $linea \t]
		}
		if {[eof $fa]} break
		puts $o $lineb\t$linea
	}
	catch {close $fb} ; catch {close $fa}
	catch {gzclose $o}
	file delete $tempbase ; file delete $tempattr
	if {$sort} {
		set f [open "|cg select -s - $tempout"]
		set o [wgzopen $outfile]
		fcopy $f $o
		gzclose $f ; gzclose $o
		file delete $tempout
	}
}
