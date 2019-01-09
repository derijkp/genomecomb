#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

# todo:
# genotype in haploid calls (Y chromosome)

proc cg_vcf2sft {args} {
	cg_vcf2tsv {*}$args
}

proc cg_vcf2tsv {args} {
	set split 1
	set sort 1
	set refout 0
	set typelist ". AD R RPA R AC A AF A"
	set collapse 0
	set removefields {}
	set keepfields *
	set locerror error
	set skiprefindels 0
	set meta {}
	cg_options vcf2tsv args {
		-s - -split {
			if {$value eq "ori"} {
				set split 0
				set collapse 0
			} elseif {[true $value]} {
				set split 1
				set collapse 0
			} else {
				set split 1
				set collapse 1
			}			
		}
		-sort {
			set sort [true $value]
		}
		-refout {
			# deletions and insertians start with a reference base
			# output a ref allele for this (if not overlapping a snp)
			set refout $value
		}
		-r - -removefields {
			set removefields $value
		}
		-t - -typelist {
			set typelist $value
		}
		-locerror {
			if {$locerror ni "error keep correct"} {error "wrong value $value for -locerror, should be one of: error keep correct"}
			set locerror $value
		}
		-keepfields {
			set keepfields $value
		}
		-skiprefindels {
			# sam varall sometimes contains long ref (alt =.) indicated as INDEL
			# they overlap with "snp" refs, causing less good results when used as a varall
			# this option to remove these
			set skiprefindels $value
		}
		-meta {
			set meta $value
		}
	} {infile outfile} 0 2
	set tempfile [tempfile]
	if {[info exists infile]} {
		# puts [list vcf2tsv $split $typelist $infile - $removefields $refout $keepfields $locerror $skiprefindels $tempfile $meta > tmp/out.tsv]
		set pipe [list exec {*}[gzcat $infile] $infile | vcf2tsv $split $typelist - - $removefields $refout $keepfields $locerror $skiprefindels $tempfile $meta]
	} else {
		set pipe [list exec vcf2tsv $split $typelist - - $removefields $refout $keepfields $locerror $skiprefindels $tempfile $meta]
	}
	if {$sort} {
		lappend pipe | cg select -s -
	}
	if {$collapse} {
		lappend pipe | cg collapsealleles
	}
	if {[info exists outfile]} {
		set compress [compresspipe $outfile]
		if {$compress ne ""} {
			lappend pipe {*}$compress
		}
		lappend pipe > $outfile 2>@ stderr
	} else {
		lappend pipe >@ stdout 2>@ stderr
	}
	# putsvars pipe
	set error [catch $pipe msg opt]
	if {$error} {
		if {$::errorCode eq "NONE"} return
		if {[string match {CHILDKILLED * SIGPIPE *} $::errorCode]} return
		# puts stderr [list set ::errorCode $::errorCode \; set msg $msg \; set opt $opt]
		dict unset opt -level
		return -options $opt "error converting vcf file: $msg"
	}
}

proc vcf2tsv_header {f {samplesVar {}} {commentVar {}} {formatfieldsVar {}} {infofieldsVar {}}} {
	if {$samplesVar ne ""} {upvar $samplesVar samples}
	if {$commentVar ne ""} {upvar $commentVar comment}
	if {$formatfieldsVar ne ""} {upvar $formatfieldsVar formatfields}
	if {$infofieldsVar ne ""} {upvar $infofieldsVar infofields}
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
	return $nheader
}

proc cg_vcf2tsv.old {args} {
	if {([llength $args] < 0) || ([llength $args] > 2)} {
		errorformat vcf2tsv
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
	set nheader [vcf2tsv_header $f samples comment formatfields infofields]
	puts -nonewline $o $comment
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
	if {$f ne "stdin"} {catch {gzclose $f}}
}

if 0 {
	set filename data/test.vcf
	set filename data/test1000glow.vcf
	set args $filename
	set f [gzopen $filename]
	set o stdout

	cg vcf2tsv $filename temp.tsv
}
