#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

# todo:
# genotype in haploid calls (Y chromosome)

proc cg_gff2sft {args} {
	cg_gff2tsv {*}$args
}

proc cg_gff2tsv {args} {
	if {([llength $args] < 0) || ([llength $args] > 2)} {
		errorformat gff2tsv
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

	set comment {# -- sft converted from gff, original comments follow --}
	while {![eof $f]} {
		set line [gets $f]
		if {[string index $line 0] ne "#"} break
		append comment \n$line
	}
	append comment "\n# ----"
	set nheader {chromosome type begin end strand source phase}
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
	while {![eof $f]} {
		if {[string index $line 0] eq "#"} continue
		set line [gets $f]
		set line [split $line \t]
		if {![llength $line]} continue
		foreach {chrom source type start end score strand phase attributes comments} $line break
		incr start -1
		set a [dict create {*}[string_change $attributes {; " " = " "}]]
		puts $fb [join [list $chrom $type $start $end $strand $source $phase] \t]
		set attrlist $attrtemplate
		dict for {key value} $a {
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
		if {$num >= $next} {putsprogress $curchromosome:$curbegin-$curend; incr next 100000}
		incr num
	}
	catch {close $fb} ; catch {close $fa}

	putsprogress "Assembling file"
	puts -nonewline $o $comment
	puts $o [join [list_concat $nheader $attrheader] \t]
	set fb [open $tempbase]
	set fa [open $tempattr]
	while 1 {
		set linea [gets $fb]
		set lineb [gets $fa]
		if {[eof $fa]} break
		puts $o $linea\t$lineb
	}
	catch {close $fb} ; catch {close $fa}
	file delete $tempbase ; file delete $tempattr
	if {$o ne "stdout"} {catch {close $o}}
	if {$f ne "stdin"} {catch {gzclose $f}}
}