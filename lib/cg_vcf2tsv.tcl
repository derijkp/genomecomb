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

proc cg_vcf2tsv {args} {
	set splitalt 0
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-s {
				set splitalt [true $value]
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	set len [llength $args]
	set tempfile [scratchfile get]
	if {$len == 0} {
		exec vcf2tsv $splitalt | cg select -s - <@ stdin >@ stdout
	} elseif {$len == 1} {
		set infile [lindex $args 0]
		exec vcf2tsv $splitalt [gztemp $infile] | cg select -s - >@ stdout
	} elseif {$len == 2} {
		set infile [lindex $args 0]
		exec vcf2tsv $splitalt [gztemp $infile] | cg select -s - > [lindex $args 1]
	} else {
		errorformat vcf2tsv
		exit 1
	}
}

proc cg_vcf2sft {args} {
	cg_vcf2tsv {*}$args
}

proc cg_vcf2sft.old {args} {
	if {([llength $args] < 0) || ([llength $args] > 2)} {
		errorformat vcf2tsv
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
	set line [gets $f]
	if {![string match "##fileformat=VCF*" $line]} {
		error "input is not a vcf file"
	}
	array set conv_formata {
		AD alleledepth
		GT genotype
		DP coverage
		FT filter
		GL loglikelihood
		GQ genoqual
		HQ haploqual
		AN totalallelecount

		AC allelecount
		AF frequency
		AA Ancestralallele
		DB dbsnp
		H2 Hapmap2
	}
	unset -nocomplain a
	set comment {# -- sft converted from vcf, original comments follow --}
	while {![eof $f]} {
		set line [gets $f]
		append comment \n$line
		if {[string range $line 0 1] eq "##"} {
			regexp {##([^=]+)=(.*)$} $line temp key value
			lappend a($key) $value
		} elseif {[string index $line 0] eq "#"} {
			set header [string range $line 1 end]
			break
		}
	}
	append comment "\n# ----"
	puts $o $comment
	set samples [lrange $header 9 end]
	set nheader {chromosome begin end type ref alt name quality filter}
	set formatfields {GT}
	set headerfields {alleleSeq1 alleleSeq2 zyg phased}
	foreach temp [get a(FORMAT) ""] {
		regexp {ID=([^,]+)} $temp temp id
		if {[inlist {GT} $id]} continue
		lappend formatfields $id
		lappend headerfields [get conv_formata($id) $id]
	}
	if {[llength $samples] == 1} {
		lappend nheader {*}$headerfields
	} else {
		foreach sample $samples {
			foreach field $headerfields {
				lappend nheader $field-$sample
			}
		}
	}
	set infofields {}
	set num 0
	foreach temp $a(INFO) {
		regexp {ID=([^,]+)} $temp temp id
		lappend infofields $id
		if {$id eq "DP"} {
			lappend nheader totalcoverage
		} else {
			lappend nheader [get conv_formata($id) $id]
		}
		incr num
	}
	set extract [list_fill $num 0 3]
	puts $o [join $nheader \t]
	set next 100000; set num 0
	while {![eof $f]} {
		if {$num >= $next} {putsprogress $num; incr next 100000}
		set line [gets $f]
		if {[string index $line 0] eq "#"} continue
		set line [split $line \t]
		if {![llength $line]} continue
		foreach {chrom pos id ref alt qual filter info format} $line break
		set format [split $format :]
		set genos [lrange $line 9 end]
		set l1 [string length $ref]
		set l2 0
		set alts [split $alt ,]
		foreach calt $alts {
			if {[string length $calt] > $l2} {set l2 [string length $calt]}
		}
		if {$l1 == 1 && $l2 == 1} {
			set begin [expr {$pos - 1}]
#			if {[catch {set begin [expr {$pos - 1}]}]} {
#				regsub -all {[^0-9]} $pos {} pos
#				set begin [expr {$pos - 1}]
#			}
			set end $pos
			set type snp
		} else {
			set begin $pos
			set end [expr {$pos + $l1 - 1}]
			if {$l1 > 20} {
				set ref [expr {$l1 - 1}]
			} else {
				set ref [string range $ref 1 end]
			}
			set temp {}
			foreach calt $alts {
				lappend temp [string range $calt 1 end]
			}
			set alts $temp
			set alt [join $alts ,]
			if {$l1 == 1} {
				set type ins
			} elseif {$l2 == 1} {
				set type del
			} else {
				set type sub
			}
		}
		set genotypes [list $ref {*}$alts]
		set result [list $chrom $begin $end $type $ref $alt $id $qual $filter]
		set order [list_cor $format $formatfields]
		foreach sample $samples geno $genos {
			set temp [list_sub [split $geno :] $order]
			if {[lindex $order 0] eq "-1"} {
				# no genotype given -> ref call
				lappend result $ref $ref r 0 {*}[lrange $temp 1 end]
				continue
			}
			if {[string first / $temp] != -1} {set phased 0} else {set phased 1}
			set genotype [split [lindex $temp 0] {|/}]
			set a1 [lindex $genotype 0]
			if {$a1 eq "."} {
				set a1 ?
			} elseif {$a1 eq ""} {
				set a1 .
			} else {
				set a1 [lindex $genotypes $a1]
			}
			set a2 [lindex $genotype 1]
			if {$a2 eq "."} {
				set a2 ?
			} elseif {$a2 eq ""} {
				set a2 .
			} else {
				set a2 [lindex $genotypes $a2]
			}
			if {$a1 ne $ref} {
				if {$a2 eq $a1} {
					set zyg m
				} elseif {$a2 ne $ref} {
					set zyg c
				} else {
					set zyg t
				}
			} elseif {$a2 ne $ref} {
				set zyg t
			} else {
				set zyg r
			}
			lappend result $a1 $a2 $zyg $phased {*}[lrange $temp 1 end]
		}
		set dinfo [dict create]
		foreach {temp key value} [regexp -all -inline {([^;=]+)=?([^;=]*)} $info] {
			dict set dinfo $key $value
		}
		foreach field $infofields {
			if {[dict exists $dinfo $field]} {
				set v [dict get $dinfo $field]
				if {$v eq ""} {set v 1}
				lappend result $v
			} else {
				lappend result {}
			}
		}
		puts $o [join $result \t]
	}	
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
	cg_vcf2sft {*}$argv
}

if 0 {
	set filename data/test.vcf
	set filename data/test1000glow.vcf
	set args $filename
	set f [gzopen $filename]
	set o stdout

	cg vcf2sft $filename temp.tsv
}
