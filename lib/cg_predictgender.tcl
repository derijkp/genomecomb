#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_predictgender {args} {
	set varfile {}
	set targetfile {}
	set xreg ""
	set yreg ""
	set refreg ""
	set dbdir ""
	cg_options predictgender args {
		-dbdir {
			set dbdir [file_absolute $value]
		}
		-varfile {
			set varfile $value
		}
		-targetfile {
			set targetfile $value
		}
		-refreg {
			set refreg $value
		}
		-xreg {
			set xreg $value
		}
		-yreg {
			set yreg $value
		}
	} {bamfile outfile} 1 2 {
		predicts gender based on bam and var file, also gives some metrics indicative of gender
	}
	set dbdir [dbdir $dbdir]
	if {$dbdir eq ""} {
		if {$xreg eq ""} {error "no xreg, give either -dbdir or -xreg option to specify non-pseudoautosomal region on X"}
		if {$yreg eq ""} {error "no yreg, give either -dbdir or -yreg option to specify non-pseudoautosomal region on Y"}
	} else {
		if {$xreg eq ""} {
			set c [split [cg select -q {$chromosome in "chrX X"} -f {begin end} [gzfile $dbdir/extra/reg_*_pseudoautosomal.tsv]] \n]
			set xreg chrX:[lindex $c 1 1]-[lindex $c 2 0]
		}
		if {$yreg eq ""} {
			set c [split [cg select -q {$chromosome in "chrY Y"} -f {begin end} [gzfile $dbdir/extra/reg_*_pseudoautosomal.tsv]] \n]
			set yreg chrY:[lindex $c 1 1]-[lindex $c 2 0]
		}
		if {$refreg eq ""} {
			set c [split [cg select -q {$chromosome in "chr22 22"} -f {begin end size} [gzfile $dbdir/extra/reg_*_sequencedgenome.tsv]] \n]
			set c [lindex [lsort -index 2 -dict $c] end-1]
			set refreg chr22:[lindex $c 0]-[lindex $c 1]
		}
	}
	if {[file isdir $bamfile]} {
		set dir $bamfile
		set bamfile [lindex [bsort [glob -nocomplain $dir/map-*.bam $dir/map-*.cram]] 0]
		set sample [file tail $dir]
	} else {
		set dir [file dir $bamfile]
		set sample [lindex [split [file root [file tail $bamfile]] -] end]
	}
	if {$varfile eq ""} {
		set varfile [gzfile $dir/var-gatk-*.tsv]
		if {![file exists $varfile]} {set varfile [gzfile $dir/var-*.tsv]}
		if {![file exists $varfile]} {set varfile ""}
	}
	if {$targetfile eq ""} {
		set targetfile [gzfile $dir/reg_*targets*.tsv]
		if {![file exists $targetfile]} {set targetfile ""}
	}
	if {$targetfile ne ""} {
		set refsize [lindex [exec cg select -q "region(\"$refreg\") == 1" $targetfile | cg covered] end]
		set xsize [lindex [exec cg select -q "region(\"$xreg\") == 1" $targetfile | cg covered] end]
		set ysize [lindex [exec cg select -q "region(\"$yreg\") == 1" $targetfile | cg covered] end]
		set tempxreg [tempfile]
		set tempyreg [tempfile]
		exec cg select -f {chromosome begin end} -sh /dev/null -q "region(\"$xreg\") == 1" $targetfile | head -500 > $tempxreg
		exec cg select -f {chromosome begin end} -sh /dev/null -q "region(\"$yreg\") == 1" $targetfile | head -500 > $tempyreg
		set end [lindex [exec tail -1 $tempxreg] end]
		if {![isint $end]} {
			set xcov ?
		} else {
			set xcov [median [exec samtools depth -a -r [lindex [split $xreg -] 0]-$end -b $tempxreg $bamfile | cut -d \t -f 3]]
		}
		set end [lindex [exec tail -1 $tempyreg] end]
		if {![isint $end]} {
			set ycov ?
		} else {
			set ycov [median [exec samtools depth -a -r [lindex [split $yreg -] 0]-$end -b $tempyreg $bamfile | cut -d \t -f 3]]
		}
	} else {
		if {![file exists $dbdir]} {
			error "no targetfile was given or found in sampledir, -dbdir must be given"
		}
		if {![regexp {([0-9]+)-([0-9]+)$} $refreg temp x1 x2]} {
			error "refreg has wrong format: $refreg"
		}
		set refsize [expr {$x2-$x1}]
		#
		if {![regexp {(.+):([0-9]+)-([0-9]+)$} $xreg temp chr x1 x2]} {
			error "xreg has wrong format: $xreg"
		}
		set xsize [expr {$x2-$x1}]
		set tempxreg [tempfile]
		set rfile [gzfile $dbdir/gene_*_refGene*.tsv]
		exec cg select -q "region(\"$xreg\") == 1" $rfile | cg gene2reg | cg select -q {$type eq "CDS"} -f {chromosome begin end} -sh /dev/null | head -500 > $tempxreg
		set end [lindex [exec tail -1 $tempxreg] end]
		if {![isint $end]} {
			set xcov ?
		} else {
			set xcov [median [exec samtools depth -a -r [lindex [split $xreg -] 0]-$end -b $tempxreg $bamfile | cut -d \t -f 3]]
		}
		#
		if {![regexp {(.+):([0-9]+)-([0-9]+)$} $yreg temp chr y1 y2]} {
			error "xreg has wrong format: $yreg"
		}
		set ysize [expr {$y2-$y1}]
		set tempyreg [tempfile]
		exec cg select -q "region(\"$yreg\") == 1" $rfile | cg gene2reg | cg select -q {$type eq "CDS"} -f {chromosome begin end} -sh /dev/null | head -500 > $tempyreg
		set end [lindex [exec tail -1 $tempyreg] end]
		if {![isint $end]} {
			set ycov ?
		} else {
			set ycov [median [exec samtools depth -a -r [lindex [split $yreg -] 0]-$end -b $tempyreg $bamfile | cut -d \t -f 3]]
		}
	}
	if {![isdouble $xcov] || $xcov < 4} {
		set yxcovratio ?
	} else {
		set yxcovratio [expr {double($ycov)/$xcov}]
	}
	set refcount [exec samtools view -q 20 -c $bamfile $refreg]
	set xcount [exec samtools view -q 20 -c $bamfile $xreg]
	set ycount [exec samtools view -q 20 -c $bamfile $yreg]
	set refncount [expr {$refcount/double($refsize)}]
	if {$xsize != 0} {
		set xncount [expr {$xcount/double($xsize)}]
	} else {
		set xncount ?
	}
	if {$ysize != 0} {
		set yncount [expr {$ycount/double($ysize)}]
	} else {
		set yncount ?
	}
	if {![isdouble $xcount] || $xcount == 0 || ![isdouble $ycount]} {
		set yxratio ?
	} else {
		set yxratio [expr {double($ycount)/$xcount}]
	}
	if {![isdouble $xncount] || $xncount == 0 || ![isdouble $yncount]} {
		set yxnratio ?
	} else {
		set yxnratio [expr {double($yncount)/$xncount}]
	}
	set pctheterozygous ?
	set pcthqheterozygous ?
	set hthqvars ?
	set totalvars 0
	set htvars 0
	if {$varfile ne ""} {
		catch {
			# pctheterozygous
			if {"zyg" ni [cg select -h $varfile]} {
				set f {zyg=zyg($alleleSeq1,$alleleSeq2,$ref,$alt)}
			} else {
				set f {}
			}
			set c [exec cg select -f $f -q "region(\"$xreg\") and \$coverage >= 20" -g zyg $varfile]
			array set a $c
			set htvars [expr {[get a(t) 0] + [get a(c) 0]}]
			set totalvars [expr {$htvars + [get a(m) 0]}]
			set pctheterozygous [format %.4f [expr {100.0*$htvars/$totalvars}]]
			# pctheterozygous
			set c [exec cg select -f $f -q "region(\"$xreg\") and \$coverage >= 20 and \$quality >= 50" -g zyg $varfile]
			array set a $c
			set hthqvars [expr {$a(t) + $a(c)}]
			set totalvars [expr {$hthqvars + $a(m)}]
			set pcthqheterozygous [format %.4f [expr {100.0*$hthqvars/$totalvars}]]
		}
	}
	if {[isdouble $yxcovratio] && $yxcovratio < 0.1} {
		set predgender f
	} elseif {[isdouble $yxcovratio] && $yxcovratio > 0.3} {
		set predgender m
	} elseif {[isdouble $yxnratio] && $yxnratio > 0.9} {
		if {[isdouble $pcthqheterozygous] && $pcthqheterozygous > 60} {
			set predgender u
		} else {
			set predgender m
		}
	} elseif {[isdouble $yxnratio] && $yxnratio < 0.5} {
		if {[isdouble $pcthqheterozygous] && $pcthqheterozygous < 40} {
			set predgender u
		} else {
			set predgender f
		}
	} elseif {[isdouble $yxnratio] && $yxnratio > 0.7 && [isdouble $pcthqheterozygous] && $pcthqheterozygous < 50} {
		set predgender m
	} elseif {[isdouble $yxnratio] && $yxnratio < 0.6 && [isdouble $pcthqheterozygous] && $pcthqheterozygous > 50} {
		set predgender f
	} else {
		set predgender u
	}
	# write results
	if {![info exists outfile]} {set o stdout} else {set o [open $outfile w]}
	puts $o [join {sample source parameter value} \t]
	puts $o $sample\tgenomecomb\tpg_xtotalvars\t$totalvars
	puts $o $sample\tgenomecomb\tpg_xheterozygous\t$htvars
	puts $o $sample\tgenomecomb\tpg_xhqheterozygous\t$hthqvars
	puts $o $sample\tgenomecomb\tpg_pctheterozygous\t$pctheterozygous
	puts $o $sample\tgenomecomb\tpg_refsize\t$refsize
	puts $o $sample\tgenomecomb\tpg_xsize\t$xsize
	puts $o $sample\tgenomecomb\tpg_ysize\t$ysize
	puts $o $sample\tgenomecomb\tpg_refcount\t$refcount
	puts $o $sample\tgenomecomb\tpg_xcount\t$xcount
	puts $o $sample\tgenomecomb\tpg_ycount\t$ycount
	puts $o $sample\tgenomecomb\tpg_refncount\t[oformat $refncount 4]
	puts $o $sample\tgenomecomb\tpg_xncount\t[oformat $xncount 4]
	puts $o $sample\tgenomecomb\tpg_yncount\t[oformat $yncount 4]
	puts $o $sample\tgenomecomb\tpg_yxratio\t[oformat $yxratio 4]
	puts $o $sample\tgenomecomb\tpg_pcthqheterozygous\t[oformat $pcthqheterozygous 4]
	puts $o $sample\tgenomecomb\tpg_yxnratio\t[oformat $yxnratio 4]
	puts $o $sample\tgenomecomb\tpg_xcov\t$xcov
	puts $o $sample\tgenomecomb\tpg_ycov\t$ycov
	puts $o $sample\tgenomecomb\tpg_yxcovratio\t[oformat $yxcovratio 3]
	puts $o $sample\tgenomecomb\tpredgender\t$predgender
	if {$o ne "stdout"} {close $o}
}

# for hg19
# set xreg chrX:2699520-154931044
# set yreg chrY:2649520-59034050
# set refreg chr22:20609431-50364777
