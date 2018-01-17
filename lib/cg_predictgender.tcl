#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_predictgender {args} {
	set bamfile {}
	set varfile {}
	set xreg ""
	set yreg ""
	set refreg ""
	set dbdir ""
	cg_options predictgender args {
		-dbdir {
			set dbdir [file_absolute $value]
		}
		-bamfile {
			set bamfile $value
		}
		-varfile {
			set varfile $value
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
	} {file outfile} 1 2 {
		predicts gender based on bam and var file, also gives some metrics indicative of gender
	}
	set dbdir [dbdir $dbdir]
	if {$dbdir eq ""} {
		if {$xreg eq ""} {set xreg chrX:2699520-154931044}
		if {$yreg eq ""} {set yreg chrY:2649520-59034050}
		if {$refreg eq ""} {set refreg chr22:20609431-50364777}
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
	set dir $file
	set sample [file tail $dir]
	if {$bamfile eq ""} {
		set bamfile [lindex [glob -nocomplain $dir/map-*.bam] 0]
	}
	if {$varfile eq ""} {
		set varfile [gzfile $dir/var-gatk-*.tsv]
		if {$varfile eq ""} {set varfile [gzfile $dir/var-*.tsv]}
	}
	#
	set targetfile [gzfile $dir/reg_*targets*.tsv]
	if {$targetfile ne ""} {
		set refsize [lindex [exec cg select -q "region(\"$refreg\") == 1" $targetfile | cg covered] end]
		set xsize [lindex [exec cg select -q "region(\"$xreg\") == 1" $targetfile | cg covered] end]
		set ysize [lindex [exec cg select -q "region(\"$yreg\") == 1" $targetfile | cg covered] end]
	} else {
		if {![regexp {([0-9]+)-([0-9]+)$} $refreg temp x1 x2]} {
			error "refreg has wrong format: $refreg"
		}
		set refsize [expr {$x2-$x1}]
		if {![regexp {([0-9]+)-([0-9]+)$} $xreg temp x1 x2]} {
			error "xreg has wrong format: $xreg"
		}
		set xsize [expr {$x2-$x1}]
		if {![regexp {([0-9]+)-([0-9]+)$} $yreg temp y1 y2]} {
			error "xreg has wrong format: $yreg"
		}
		set ysize [expr {$y2-$y1}]
	}
	set refcount [exec samtools view -q 5 -c $bamfile $refreg]
	set xcount [exec samtools view -q 5 -c $bamfile $xreg]
	set ycount [exec samtools view -q 5 -c $bamfile $yreg]
	set refncount [expr {$refcount/double($refsize)}]
	set xncount [expr {$xcount/double($xsize)}]
	set yncount [expr {$ycount/double($ysize)}]
	set xyratio [expr {double($ycount)/$xcount}]
	if {$xyratio <= 0.03} {
		set xygender f
	} elseif {$xyratio >= 0.05} {
		set xygender m
	} else {
		set xygender u
	}
	# pctheterozygous
	set c [exec cg select -q "region(\"$xreg\") and \$coverage >= 20" -g zyg $varfile]
	array set a $c
	set htvars [expr {$a(t) + $a(c)}]
	set totalvars [expr {$htvars + $a(m)}]
	set pctheterozygous [expr {100.0*$htvars/$totalvars}]
	# pctheterozygous
	set c [exec cg select -q "region(\"$xreg\") and \$coverage >= 20 and \$quality >= 50" -g zyg $varfile]
	array set a $c
	set hthqvars [expr {$a(t) + $a(c)}]
	set totalvars [expr {$hthqvars + $a(m)}]
	set pcthqheterozygous [expr {100.0*$hthqvars/$totalvars}]
	# write results
	if {![info exists outfile]} {set o stdout} else {set o [open $outfile w]}
	puts $o [join {sample source parameter value} \t]
	puts $o $sample\tgenomecomb\tpg_pctheterozygous\t[format %.4f $pctheterozygous]
	puts $o $sample\tgenomecomb\tpg_pcthqheterozygous\t[format %.4f $pctheterozygous]
	puts $o $sample\tgenomecomb\tpg_refsize\t$refsize
	puts $o $sample\tgenomecomb\tpg_xsize\t$xsize
	puts $o $sample\tgenomecomb\tpg_ysize\t$ysize
	puts $o $sample\tgenomecomb\tpg_refcount\t$refcount
	puts $o $sample\tgenomecomb\tpg_xcount\t$xcount
	puts $o $sample\tgenomecomb\tpg_ycount\t$ycount
	puts $o $sample\tgenomecomb\tpg_refncount\t[format %.4f $refncount]
	puts $o $sample\tgenomecomb\tpg_xncount\t[format %.4f $xncount]
	puts $o $sample\tgenomecomb\tpg_yncount\t[format %.4f $yncount]
	puts $o $sample\tgenomecomb\tpg_xyratio\t[format %.4f $xyratio]
	puts $o $sample\tgenomecomb\tpg_xheterozygous\t$htvars
	puts $o $sample\tgenomecomb\tpg_xtotal\t$totalvars
	puts $o $sample\tgenomecomb\txygender\t$xygender
	if {$o ne "stdout"} {close $o}
}
