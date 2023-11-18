proc sc_pseudobulk_job {args} {
	upvar job_logdir job_logdir
	cg_options sc_pseudobulk_job args {
	} {scgenefile scisoformfile groupfile} 3 3
	#
	set scgenefile [file_absolute $scgenefile]
	set scisoformfile [file_absolute $scisoformfile]
	set groupfile [file_absolute $groupfile]
	set rootname [file_rootname $groupfile]
	set dir [file dir $scgenefile]
	set pb_genefile $dir/pb_gene_counts-$rootname.tsv.zst
	set pb_isoformfile $dir/pb_isoform_counts-$rootname.tsv.zst
	job sc_pseudobulk-$rootname -deps {
		$scgenefile $scisoformfile $groupfile
	} -targets {
		$pb_genefile $pb_isoformfile
	} -vars {
		scgenefile scisoformfile groupfile pb_genefile pb_isoformfile rootname
	} -code {
		# read group file
		set f [gzopen $groupfile]
		set header [tsv_open $f]
		set cellpos [lsearch $header cell]
		if {$cellpos == -1} {error "group file $groupfile does not have a cell field"}
		set grouppos [lsearch $header group]
		if {$grouppos == -1} {error "group file $groupfile does not have a group field"}
		set poss [list $cellpos $grouppos]
		set infoposs [list_cor $header {group_filtered score ncells}]
		set infoposs [list_sub $infoposs -exclude [list_find $infoposs -1]]
		set colinfofields [list_sub $header $infoposs]
		unset -nocomplain a
		unset -nocomplain ga
		while {[gets $f line] != -1} {
			set line [split $line \t]
			foreach {cell group} [list_sub $line $poss] break
			set group [regsub -all " " $group _]
			set a($cell) $group
			if {![info exists ga($group)]} {
				set ga($group) [list_sub $line $infoposs]
			}
		}
		gzclose $f

		# write pseudobulk gene file
		catch {close $f} ; catch {close $o}
		set tempfile [tempfile]
		cg select -s geneid $scgenefile $tempfile
		set temppb_genefile $pb_genefile.temp[gzext $pb_genefile]
		set o [wgzopen $temppb_genefile]
		set f [open $tempfile]
		set header [tsv_open $f]
		set poss [list_sub [tsv_basicfields $header 7 0] {0 1 2 6}]
		set common {chromosome begin end strand}
		set temp [list_find $poss -1]
		set newheader [list_sub $common -exclude $temp]
		set coltypes [list_fill [llength $newheader] info]
		set poss [list_sub $poss -exclude $temp]
		lappend newheader gene geneid
		lappend coltypes info id
		lappend poss {*}[list_cor $header {gene geneid}]
		set idposs [list_cor $header {geneid cell}]
		set datafields [list_remove $header {*}$newheader cell]
		set dataposs [list_cor $header $datafields]
		set types [bsort [array names ga]]
		foreach datafield $datafields {
			foreach type $types {
				lappend newheader $datafield-$type-$rootname
			}
		}
		puts $o [join $newheader \t]
		set curgeneid {}
		foreach datafield $datafields {
			foreach type $types {
				set totala($datafield-$type) 0
			}
		}
		while 1 {
			if {[gets $f line] == -1} break
			set line [split $line \t]
			foreach {geneid cell} [list_sub $line $idposs] break
			if {$geneid ne $curgeneid} {
				if {$curgeneid ne ""} {
					set resultline $curgeneinfo
					foreach datafield $datafields {
						foreach type $types {
							if {[string is int $da($datafield-$type)]} {
								lappend resultline $da($datafield-$type)
							} else {
								lappend resultline [format %.1f $da($datafield-$type)]
							}
							set totala($datafield-$type) [expr {$totala($datafield-$type)+$da($datafield-$type)}]
						}
					}
					puts $o [join $resultline \t]
				}
				set curgeneid $geneid
				set curgeneinfo [list_sub $line $poss]
				foreach datafield $datafields {
					foreach type $types {
						set da($datafield-$type) 0
					}
				}
			}
			set type $a($cell)
			foreach datafield $datafields value [list_sub $line $dataposs] {
				set da($datafield-$type) [expr {$da($datafield-$type)+$value}]
			}
		}
		gzclose $o
		gzclose $f
		set temppb_genefile2 $pb_genefile.temp2[gzext $pb_genefile]
		cg select -s - $temppb_genefile $temppb_genefile2
		file delete $temppb_genefile
		#
		set o [wgzopen [file root [gzroot $pb_genefile]].colinfo.tsv]
		set sample [lindex [split $rootname -] end]
		puts $o [join [list field fieldtype group {*}$colinfofields total sample rootname] \t]
		set empty [list_fill [expr {2+[llength $ga([lindex [array names ga] 0])]}] {}]
		foreach field $newheader coltype $coltypes {
			if {$coltype eq ""} break
			puts $o [join [list $field $coltype {*}$empty] \t]
		}
		foreach datafield $datafields {
			foreach type $types {
				set field $datafield-$type-$rootname
				puts $o [join [list $field data $type {*}$ga($type) [format %.1f $totala($datafield-$type)] $sample $rootname] \t]
			}
		}
		gzclose $o
		file rename -force $temppb_genefile2 $pb_genefile

		# write pseudobulk isoform file
		catch {close $f} ; catch {close $o}
		set tempfile [tempfile]
		cg select -overwrite 1 -s transcript $scisoformfile $tempfile
		set temppb_isoformfile $pb_isoformfile.temp[gzext $pb_isoformfile]
		set o [wgzopen $temppb_isoformfile]
		set f [open $tempfile]
		set header [tsv_open $f]
		set dataposs [list_find -regexp $header count]
		set datafields [list_sub $header $dataposs]
		set infofields [list_remove [list_sub $header -exclude $dataposs] ROW cell]
		set idposs [list_cor $header {transcript cell}]
		set infoposs [list_cor $header $infofields]
		set newheader $infofields
		set coltypes [list_fill [llength $newheader] info]
		lset coltypes [lsearch $newheader transcript] id
		set types [bsort [array names ga]]
		foreach datafield $datafields {
			foreach type $types {
				lappend newheader $datafield-$type-$rootname
			}
		}
		puts $o [join $newheader \t]
		set curtranscript {}
		foreach datafield $datafields {
			foreach type $types {
				set totala($datafield-$type) 0
			}
		}
		while 1 {
			if {[gets $f line] == -1} break
			set line [split $line \t]
			foreach {transcript cell} [list_sub $line $idposs] break
			if {$transcript ne $curtranscript} {
				if {$curtranscript ne ""} {
					set resultline $curisoforminfo
					foreach datafield $datafields {
						foreach type $types {
							if {[string is int $da($datafield-$type)]} {
								lappend resultline $da($datafield-$type)
							} else {
								lappend resultline [format %.1f $da($datafield-$type)]
							}
							set totala($datafield-$type) [expr {$totala($datafield-$type)+$da($datafield-$type)}]
						}
					}
					puts $o [join $resultline \t]
				}
				set curtranscript $transcript
				set curisoforminfo [list_sub $line $infoposs]
				foreach datafield $datafields {
					foreach type $types {
						set da($datafield-$type) 0
					}
				}
			}
			set type $a($cell)
			foreach datafield $datafields value [list_sub $line $dataposs] {
				set da($datafield-$type) [expr {$da($datafield-$type)+$value}]
			}
		}
		gzclose $o
		gzclose $f
		set temppb_isoformfile2 $pb_isoformfile.temp2[gzext $pb_isoformfile]
		cg select -s - $temppb_isoformfile $temppb_isoformfile2
		file delete $temppb_isoformfile
		#
		set o [wgzopen [file root [gzroot $pb_isoformfile]].colinfo.tsv]
		set sample [lindex [split $rootname -] end]
		puts $o [join [list field fieldtype group {*}$colinfofields total sample rootname] \t]
		set empty [list_fill [expr {2+[llength $ga([lindex [array names ga] 0])]}] {}]
		foreach field $newheader coltype $coltypes {
			if {$coltype eq ""} break
			puts $o [join [list $field $coltype {*}$empty] \t]
		}
		foreach datafield $datafields {
			foreach type $types {
				set field $datafield-$type-$rootname
				puts $o [join [list $field data $type {*}$ga($type) [format %.1f $totala($datafield-$type)] $sample $rootname] \t]
			}
		}
		gzclose $o
		file rename -force $temppb_isoformfile2 $pb_isoformfile
	}
}

proc cg_sc_pseudobulk {args} {
	set args [job_init {*}$args]
	sc_pseudobulk_job {*}$args
	job_wait
}

