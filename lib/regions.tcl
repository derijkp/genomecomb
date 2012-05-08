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
#	set chr [lindex $result 0]
#	if {![isint $chr]} {
#		lset result 0 [chr2num $chr]
#	}
#	return $result
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
		puts stderr "format is: $scriptname $action variation_file"
		puts stderr " - outputs regions annotated as ref-(in)consistent from variation file"
		exit 1
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

proc regsubtract {regfile1 regfile2} {
	global cache
	# catch {close $f1}
	# catch {close $f2}
	set f1 [open $regfile1]
	set poss1 [open_region $f1]
	set f2 [open $regfile2]
	set poss2 [open_region $f2]
	close $f1; close $f2
	# puts [list $regfile1 {*}$poss1 $regfile2 {*}$poss2]
	exec reg_subtract $regfile1 {*}$poss1 $regfile2 {*}$poss2 >@ stdout 2>@ stderr
}

proc reghisto {regfile} {
	global cache
	unset -nocomplain a
	set v [opencgifile $regfile header]
	set num 0
	while 1 {
		incr num
		if {![expr $num%100000]} {putsprogress $num}
		set vline [cggets $v]
		if {![llength $vline]} break
		foreach {vchr vstart vend} $vline break
		set len [expr {$vend-$vstart}]
		if {![info exists a($vchr,$len)]} {
			set a($vchr,$len) 1
		} else {
			incr a($vchr,$len)
		}
	}
	close $v
	unset -nocomplain tota
	set o stdout
	set names [lsort -dictionary [array names a]]
	foreach name $names {
		foreach {chr len} [split $name ,] break
		puts $o $chr\t$len\t$a($name)
		if {![info exists tota($len)]} {
			set tota($len) 1
		} else {
			incr tota($len) $a($name)
		}
	}
	set names [lsort -dictionary [array names tota]]
	set tot 0
	foreach len $names {
		puts $o total\t$len\t$tota($len)
		incr tot $tota($len)
	}
	puts $o all\tall\t$tot
	puts $o all\tallnonull\t[expr {$tot-[get tota(0) 0]}]
}

proc cg_reghisto {args} {
	global scriptname action
	if {[llength $args] != 1} {
		puts stderr "format is: $scriptname $action region_file"
		puts stderr ""
		exit 1
	}
	foreach {region_file} $args break
	reghisto $region_file
}

proc regjoin {regfile1 regfile2} {
	global cache
	# catch {close $f1}
	# catch {close $f2}
	set f1 [open $regfile1]
	set poss1 [open_region $f1]
	close $f1
	if {$regfile2 ne ""} {
		set f2 [open $regfile2]
		set poss2 [open_region $f2]
		close $f2
	} else {
		set poss2 {0 0 0}
	}
	# puts [list reg_join $regfile1 {*}$poss1 $regfile2 {*}$poss2]
	exec reg_join $regfile1 {*}$poss1 $regfile2 {*}$poss2 >@ stdout 2>@ stderr
}

proc cg_regsubtract {args} {
	if {[llength $args] != 2} {
		errorformat regsubtract
		exit 1
	}
	foreach {region_file1 region_file2} $args break
	regsubtract $region_file1 $region_file2
}

proc cg_regjoin {args} {
	if {([llength $args] != 1) && ([llength $args] != 2)} {
		errorformat regjoin
		exit 1
	}
	foreach {region_file1 region_file2} $args break
	regjoin $region_file1 $region_file2
}

proc cg_regcommon {args} {
	if {([llength $args] != 1) && ([llength $args] != 2)} {
		errorformat regcommon
		exit 1
	}
	foreach {region_file1 region_file2} $args break
	set tempfile [tempfile]
	cg regsubtract $region_file1 $region_file2 > $tempfile
	cg regsubtract $region_file1 $tempfile >@ stdout
}
