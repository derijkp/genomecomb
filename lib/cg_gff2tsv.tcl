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

	set comment {# -- tsv converted from gff, original comments follow --}
	while {![eof $f]} {
		set line [gets $f]
		if {[string index $line 0] ne "#"} break
		append comment \n$line
	}
	append comment "\n# ----\n"
	set nheader {chromosome type begin end score strand source phase}
	set next 100000; set num 0
	set tempbase [tempfile]
	set tempattr [tempfile]
	catch {close $fb} ; catch {close $fa}
	set fb [open $tempbase w]
	set fa [open $tempattr w]
	set attrheader {}
	set attrtemplate {}
 	unset -nocomplain attra
	set num 0
	foreach key {chromosome type begin end score strand source phase} {
		set transa($key) attr_$key
	}
	while 1 {
		while {$line eq "" || [string index $line 0] eq "#"} {
			if {[gets $f line] == -1} break
		}
		if {[eof $f]} break
		set line [split $line \t]
		if {![llength $line]} continue
		foreach {chrom source type start end score strand phase attributes comments} $line break
		incr start -1
		set data [split [string_change $attributes {"; " ";"}] {;=}]
		set a [dict create {*}$data]
		puts $fb [join [list $chrom $type $start $end $score $strand $source $phase] \t]
		set attrlist $attrtemplate
		dict for {key value} $a {
			if {[info exists transa($key)]} {
				set key $transa($key)
			}
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
		if {$num >= $next} {putsprogress $chrom:$start-$end; incr next 100000}
		incr num
		set line [gets $f]
	}
	catch {close $fb} ; catch {close $fa}

	putsprogress "Assembling file"
	set atrrsize [llength $attrheader]
	puts -nonewline $o $comment
	puts $o [join [list_concat $nheader $attrheader] \t]
	set fb [open $tempbase]
	set fa [open $tempattr]
	unset -nocomplain atemplate
	while 1 {
		set lineb [gets $fb]
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
	file delete $tempbase ; file delete $tempattr
	if {$o ne "stdout"} {catch {close $o}}
	if {$f ne "stdin"} {catch {gzclose $f}}
}
