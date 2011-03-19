#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

# todo:
# genotype in haploid calls (Y chromosome)

proc cg_vcf2sft_help {} {
set help [file_read $::appdir/lib/cg_vcf2sft.help]
puts [string_change $help [list @BASE@ [get ::base {[info source]}]]]
}

proc cg_vcf2sft {args} {
	if {([llength $args] < 1) || ([llength $args] > 2)} {
		puts "Wrong number of arguments"
		cg_vcf2sft_help
		exit 1
	}
	if {[llength $args] > 0} {
		set filename [lindex $args 0]
		set f [rzopen $filename]
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
		DP coverage
		GL loglikelihood
		GQ genoqual
		GT genotype
		HQ haploqual
	}
	unset -nocomplain a
	while {![eof $f]} {
		set line [gets $f]
		if {[string range $line 0 1] eq "##"} {
			regexp {##([^=]+)=(.*)$} $line temp key value
			lappend a($key) $value
		} elseif {[string index $line 0] eq "#"} {
			set header [string range $line 1 end]
			break
		}
	}
	set samples [lrange $header 9 end]
	set nheader {chromsome begin end type ref alt}
	lappend nheader quality filter
	set formatfields {GT}
	set headerfields {alleleSeq1 alleleSeq2 fased}
	foreach temp $a(FORMAT) {
		set d [split [string range $temp 1 end-1] ,=]
		set id [dict get $d ID]
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
	puts $o [join $nheader \t]
	set next 100000; set num 0
	while {![eof $f]} {
		if {$num >= $next} {putslog $num; incr next 100000}
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		foreach {chrom pos id ref alt qual filter info format} $line break
		set format [split $format :]
		set genos [lrange $line 9 end]
		set l1 [string length $ref]
		set l2 0
		set alts [split $alt ,]
		set genotypes [list_concat $ref $alts]
		foreach calt $alts {
			if {[string length $calt] > $l2} {set l2 [string length $calt]}
		}
		if {$l1 == 1 && $l2 == 1} {
			set begin [expr {$pos - 1}]
			set end $pos
			set type snp
		} else {
			set begin $pos
			set end [expr {$pos + $l1 - 1}]
			set ref [string range $ref 1 end]
			set temp {}
			foreach calt $alts {
				lappend temp [string range $calt 1 end]
			}
			set alt [join $temp ,]
			if {$l1 == 1} {
				set type ins
			} elseif {$l2 == 1} {
				set type del
			} else {
				set type sub
			}
		}
		set result [list $chrom $begin $end $type $ref $alt $qual $filter]
		foreach sample $samples geno $genos {
			set order [list_cor $format $formatfields]
			set temp [list_sub [split $geno :] $order]
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
			lappend result $a1 $a2 $phased {*}[lrange $temp 1 end]
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
	set filename tests/data/test.vcf
	set args tests/data/test.vcf
	cg vcf2sft tests/data/test.vcf tests/data/test.sft
}
