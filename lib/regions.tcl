#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require Extral

proc open_region {f {headerVar {}}} {
	if {[string length $headerVar]} {
		upvar $headerVar header
	}
	set header [tsv_open $f]
	if {[string index $header 0] eq "#"} {
		set header [string range $header 1 end]
	}
	set poss2 [tsv_basicfields $header 3]
	return $poss2
}

proc get_region {f poss} {
	while 1 {
		set line [split [gets $f] \t]
		if {[llength $line]} break
		if {[eof $f]} break
	}
	set result [list_sub $line $poss]
}

proc refconsregions {varfile} {
	putslog "getting ref-(in)consistent regions from $varfile"
	cg select -q {$varType == "ref-consistent" || $varType == "ref-inconsistent" || $varType == "no-call-rc" || $varType == "no-call-ri"} $varfile rctemp.tsv
	cg regjoin rctemp.tsv >@stdout
	file delete rctemp.tsv
}

proc cg_refconsregions {args} {
	global scriptname action
	if {[llength $args] != 1} {
		error "format is: $scriptname $action variation_file\n - outputs regions annotated as ref-(in)consistent from variation file"
	}
	foreach {varfile} $args break
	refconsregions $varfile
}

proc nocallregions {varfile outfile} {
	putslog "getting partial no-call regions from $varfile"
	set h [cg select -h $varfile]
	if {[inlist $h allele]} {
		cg select -q {$varType == "no-call" && $allele != "all"} $varfile nctemp.tsv
	} else {
		cg select -q {$varType == "no-call" && $haplotype != "all"} $varfile nctemp.tsv
	}
	cg regjoin nctemp.tsv > $outfile
	file delete nctemp.tsv
}

proc findregionfile {file} {
	set tail [file tail $file]
	regsub {^[^-]*-} $tail {} tail
	return [file dir $file]/sreg-$tail
}
