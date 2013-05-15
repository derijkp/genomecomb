#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

# todo:
# genotype in haploid calls (Y chromosome)

proc cg_gtf2sft {args} {
	if {([llength $args] < 0) || ([llength $args] > 2)} {
		errorformat gtf2sft
		exit 1
	}
	if {[llength $args] > 0} {
		set filename [lindex $args 0]
		set f [gzopen $filename]
	} else {
		set f stdin
	}
	if {[llength $args] > 1} {
		set outfile [lindex $args 1]
		set o [open $outfile w]
	} else {
		set o stdout
	}

	set comment {# -- sft converted from gtf, original comments follow --}
	while {![eof $f]} {
		set line [gets $f]
		if {[string index $line 0] ne "#"} break
		append comment \n$line
	}
	append comment "\n# ----"
	set nheader {chromosome begin end name strand cdsStart cdsEnd exonCount exonStarts exonEnds source}
	set next 100000; set num 0
	set tempbase [tempfile]
	set tempattr [tempfile]
	catch {close $fb} ; catch {close $fa}
	set fb [open $tempbase w]
	set fa [open $tempattr w]
	set attrheader {}
	set attrtemplate {}
 	unset -nocomplain attra
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
		} else {
			set transcript {}
		}
		if {![info exists curchromosome] || $transcript ne $curtranscript} {
			# write previous to output
			if {[info exists curchromosome]} {
				if {$curstrand eq "-"} {
					set curexonStarts [list_reverse $curexonStarts]
					set curexonEnds [list_reverse $curexonEnds]
				}
				puts $fb [join [list $curchromosome $curbegin $curend $curtranscript $curstrand \
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
	catch {close $fb} ; catch {close $fa}

	putsprogress "Assembling file"
	puts $o $comment
	puts $o [join [list_concat $nheader $attrheader] \t]
	set fb [open $tempbase]
	set fa [open $tempattr]
	while {![eof $fa]} {
		set linea [gets $fb]
		set lineb [gets $fa]
		puts $o $linea\t$lineb
	}
	catch {close $fb} ; catch {close $fa}
	file delete $tempbase ; file delete $tempattr
	if {$o ne "stdout"} {catch {close $o}}
	if {$f ne "stdin"} {catch {close $f}}
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base $scriptname
	cg_gtf2sft {*}$argv
}
