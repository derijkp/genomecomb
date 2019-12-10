proc cg_report_vars {args} {
	set sample {}
	set targetfile {}
	set refcodingfile {}
	cg_options reports_vars args {
		-sample {
			set sample $value
		}
		-targetfile {
			set targetfile $value
		}
		-refcodingfile {
			set refcodingfile $value
		}
	} {varfile resultfile}
	set tempfile [tempfile]
	set fields {}
	set header [cg select -h $varfile]
	if {"coverage" in $header && "quality" in $header} {
		set usehq 1
		lappend fields {hq=if(def($coverage,0) >= 20 and def($quality,0) >= 50,1,0)}
	} elseif {"coverage-$sample" in $header && "quality-$sample" in $header} {
		set usehq 1
		lappend fields "hq=if(def(\$coverage-$sample,0) >= 20 and def(\$quality-$sample,0) >= 50,1,0)"
	} else {
		set usehq 0
		lappend fields {hq=0}
	}
	set annotfiles {}
	set targetfile [gzfile $targetfile]
	if {[file exists $targetfile]} {
		set tempdir [tempdir]
		mklink $targetfile $tempdir/reg_targets.tsv[gzext $targetfile]
		lappend annotfiles $tempdir/reg_targets.tsv
		lappend fields {target=if($targets ne "",1,0)}
	} else {
		lappend fields {target=0}
	}
	set refcodingfile [gzfile $refcodingfile]
	if {[file exists $refcodingfile]} {
		lappend annotfiles $refcodingfile
		lappend fields {refcoding=if($refcoding ne "",1,0)}
	} else {
		lappend fields {refcoding=0}
	}
	if {"ref" in $header && "alt" in $header} {
		# titv ratio WGS ~ 2.0-2.1, WES ~ 3.0-3.3
		set usetitv 1
		# transition A<->G and C<->T, transversion: A<->C A<->T G<->T G<->C
		lappend fields {titv=if($type ne "snp","","$ref$alt" in "AG GA CT TC","i","v")}
	} else {
		set usetitv 0
		lappend fields {titv=""}
	}
	if {[inlist $header type]} {
		lappend fields type
	} else {
		lappend fields {type=""}
	}
	if {[inlist $header zyg]} {
		lappend fields zyg
	} elseif {[inlist $header zyg-$sample]} {
		lappend fields "zyg=\$zyg-$sample"
	} else {
		lappend fields {zyg=""}
	}
	if {[llength $annotfiles]} {
		cg annotate $varfile $tempfile {*}$annotfiles
	} else {
		mklink $varfile $tempfile
	}
	if {[inlist $header sequenced]} {
		set query {$sequenced == "v"}
	} elseif {[inlist $header sequenced-$sample]} {
		set query "\$sequenced-$sample == \"v\""
	} elseif {[inlist $header zyg]} {
		set query {$zyg in "m t c"}
	} elseif {[inlist $header zyg-$sample]} {
		set query "\$zyg-$sample in \"m t c\""
	} else {
		set query {}
	}
	set temp [exec cg select -q $query \
		-f $fields -g {hq * target * refcoding * titv * type * zyg *} $tempfile]

	set vars 0 ; set qvars 0 ; set qvars_target 0; set qvars_refcoding 0
	set vars_snp 0 ; set qvars_snp 0 ; set qvars_target_snp 0; set qvars_refcoding_snp 0
	set titva(i) 0 ; set titva(v) 0
	set qtitva(i) 0 ; set qtitva(v) 0
	unset -nocomplain typea ; unset -nocomplain typea
	unset -nocomplain zyga ; unset -nocomplain qzyga
	foreach line [lrange [split [string trim $temp] \n] 1 end] {
		foreach {hq ontarget refcoding titv type zyg count} [split $line \t] break
		if {$hq} {
			incr qvars $count
			if {$type eq "snp"} {incr qvars_snp $count}
			if {$ontarget} {
				incr qvars_target $count
				if {$type eq "snp"} {incr qvars_target_snp $count}
			}
			if {$refcoding} {
				incr qvars_refcoding $count
				if {$type eq "snp"} {incr qvars_refcoding_snp $count}
			}
			incr qtitva($titv) $count
			if {![info exists qzyga($zyg)]} {set qzyga($zyg) $count} else {incr qzyga($zyg) $count}
			if {![info exists qtypea($type)]} {set qtypea($type) $count} else {incr qtypea($type) $count}
		}
		incr vars $count
		if {$type eq "snp"} {incr vars_snp $count}
		incr titva($titv) $count
		if {![info exists zyga($zyg)]} {set zyga($zyg) $count} else {incr zyga($zyg) $count}
		if {![info exists typea($type)]} {set typea($type) $count} else {incr typea($type) $count}
	}
	set f [open $resultfile.temp w]
	puts $f [join {sample source parameter value} \t]
	puts $f $sample\tgenomecomb\tvars\t$vars
	puts $f $sample\tgenomecomb\tvars_snp\t$vars_snp
	puts $f $sample\tgenomecomb\tqvars\t$qvars
	puts $f $sample\tgenomecomb\tqvars_snp\t$qvars_snp
	if {[file exists $targetfile]} {
		puts $f $sample\tgenomecomb\tqvars_target\t$qvars_target
		puts $f $sample\tgenomecomb\tqvars_target_snp\t$qvars_target_snp
	}
	if {[file exists $refcodingfile]} {
		puts $f $sample\tgenomecomb\tqvars_refcoding\t$qvars_refcoding
		puts $f $sample\tgenomecomb\tqvars_refcoding_snp\t$qvars_refcoding_snp
	}
	proc printvarinfo {f sample name varVar} {
		upvar $varVar var
		foreach temp [lsort -dict [list_remove [array names var] {}]] {
			puts $f $sample\tgenomecomb\t$name$temp\t$var($temp)
		}
	}
	printvarinfo $f $sample vars_zyg_ zyga
	printvarinfo $f $sample vars_type_ typea
	printvarinfo $f $sample vars_titv_ titva
	if {$titva(v) != 0} {
		puts $f $sample\tgenomecomb\tvars_titv\t[format %.1f [expr double ($titva(i))/$titva(v)]]
	}
	printvarinfo $f $sample qvars_zyg_ qzyga
	printvarinfo $f $sample qvars_type_ qtypea
	printvarinfo $f $sample qvars_titv_ qtitva
	if {$qtitva(v) != 0} {
		puts $f $sample\tgenomecomb\tqvars_titv\t[format %.1f [expr double ($qtitva(i))/$qtitva(v)]]
	}
	close $f
	file rename -force -- $resultfile.temp $resultfile

}