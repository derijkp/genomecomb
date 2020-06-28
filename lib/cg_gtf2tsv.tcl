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
	cg_options gtf2tsv args {
		-separate {
			set separate $value
		}
	} {filename outfile} 0 2
	set f [gzopen $filename]

	set comment {# -- sft converted from gtf, original comments follow --}
	while {![eof $f]} {
		set line [gets $f]
		if {[string index $line 0] ne "#"} break
		append comment \n$line
	}
	append comment "\n# ----"
	set next 100000; set num 0
	set tempbase [tempfile]
	set tempattr [tempfile]
	catch {close $fb} ; catch {close $fa}
	set fb [open $tempbase w]
	set fa [open $tempattr w]
	set attrheader {}
	set attrtemplate {}
 	unset -nocomplain attra
	if {!$separate} {
		set nheader {chromosome begin end name gene strand cdsStart cdsEnd exonCount exonStarts exonEnds source}
		unset -nocomplain curchromosome
		set num 0
		while 1 {
			if {[string index $line 0] eq "#"} continue
			set line [split $line \t]
			set eof [eof $f]
			if {!$eof} {
				if {![llength $line]} continue
				foreach {chrom source type start end score strand phase attributes comments} $line break
				if {![inlist {exon CDS start_codon stop_codon 5UTR 3UTR inter inter_CNS intron_CNS} $type]} {
					set line [gets $f]
					continue
				}
				incr start -1
				set a [dict create {*}[string_change $attributes {; " "}]]
				set transcript [dict get $a transcript_id]
				set gene [dict get $a gene_id]
			} else {
				set transcript {}
				set gene {}
			}
			if {![info exists curchromosome] || $transcript ne $curtranscript} {
				# write previous to output
				if {[info exists curchromosome]} {
					if {$curstrand eq "-"} {
						set curexonStarts [list_reverse $curexonStarts]
						set curexonEnds [list_reverse $curexonEnds]
					}
					puts $fb [join [list $curchromosome $curbegin $curend $curtranscript $curgene $curstrand \
						$curcdsStart $curcdsEnd $curexonCount \
						[join $curexonStarts ,] [join $curexonEnds ,] $cursource] \t]
					set attrlist $attrtemplate
					dict for {key value} $curattr {
						if {![info exists attra($key)]} {
							set attra($key) [llength $attrtemplate]
							lappend attrheader $key
							lappend attrtemplate {}
							lappend attrlist $value
						} else {
							lset attrlist $attra($key) $value
						}
					}
					puts $fa [join $attrlist \t]
				}
				if {$eof} break
				# start next
				set curtranscript $transcript
				set curgene $gene
				set curchromosome $chrom
				set curbegin $start
				set curend $end
				set curstrand $strand
				set curcdsStart {}
				set curcdsEnd {}
				set curexonCount 0
				set curexonStarts {}
				set curexonEnds {}
				set cursource $source
				set curscore {}
				set curattr $a
			}
			if {$num >= $next} {putsprogress $curchromosome:$curbegin-$curend; incr next 100000}
			incr num
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
						set curcdsStart $start
						set curcdsEnd [expr {$end+3}]
					} else {
						set curcdsStart [expr {$start-3}]
						set curcdsEnd $end
					}
					if {$start < $curbegin} {set curbegin $start}
					if {$end > $curend} {set curend $end}
				}
				start_codon {
					if {$strand eq "+"} {
						set curcdsStart $start
					} else {
						set curcdsEnd $end
					}
				}
				stop_codon {
					if {$strand eq "+"} {
						set curcdsEnd $end
					} else {
						set curcdsStart $start
					}
				}
			}
			set curattr [dict merge $curattr $a]
			set line [gets $f]
		}
	} else {
		set nheader {chromosome begin end type transcript gene strand source comments}
		set num 0
		set curattr [dict create]
		while 1 {
			if {[string index $line 0] eq "#"} continue
			set line [split $line \t]
			if {![llength $line]} continue
			foreach {chrom source type start end score strand phase attributes comments} $line break
			incr start -1
			set attributes [dict create {*}[string_change $attributes {; " "}]]
			set transcript [dict get $attributes transcript_id]
			set gene [dict get $attributes gene_id]
			puts $fb [join [list $chrom $start $end $type $transcript $gene $strand $source $comments] \t]
			set attrlist $attrtemplate
			dict for {key value} $attributes {
				if {![info exists attra($key)]} {
					set attra($key) [llength $attrtemplate]
					lappend attrheader $key
					lappend attrtemplate {}
					lappend attrlist $value
				} else {
					lset attrlist $attra($key) $value
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
	puts $o $comment
	puts $o [join [list_concat $nheader $attrheader] \t]
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
	catch {close $fb} ; catch {close $fa}
	file delete $tempbase ; file delete $tempattr
	if {$o ne "stdout"} {catch {close $o}}
}
