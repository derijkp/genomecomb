#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc multicompar_annot_join {cur1 cur2} {
	global joinposs1 mergeposs1 joinposs2 mergeposs2 dummy1 dummy2 restposs1 restposs2 refpos1 refpos2 altpos1 altpos2 alleleposs1 alleleposs2 listfields1 listfields2 sequenced2pos
	if {[inlist {{} -} $cur1]} {
		set region [list_sub $cur2 $joinposs2]
		set merge [list_sub $cur2 $mergeposs2]
		if {$cur1 eq "-"} {
			set cur1 $dummy1
		} else {
			set cur1 [list_change $dummy1 {- {}}]
		}
		if {$sequenced2pos != -1} {
			set sequenced [lindex $cur2 $sequenced2pos]
		} else {
			set sequenced v
		}
	} elseif {[inlist {{} -} $cur2]} {
		set region [list_sub $cur1 $joinposs1]
		set merge [list_sub $cur1 $mergeposs1]
		if {$cur2 eq "-"} {
			set cur2 $dummy2
		} else {
			set cur2 [list_change $dummy2 {- {}}]
		}
		set sequenced ?
	} else {
		set region [list_sub $cur1 $joinposs1]
		set merge {}
		foreach el1 [list_sub $cur1 $mergeposs1] el2 [list_sub $cur2 $mergeposs2] {
			lappend merge [list_union $el1 $el2]
		}
		if {$sequenced2pos != -1} {
			set sequenced [lindex $cur2 $sequenced2pos]
		} else {
			set sequenced v
		}
	}
	if {$refpos1 == -1} {
		if {$refpos2 == -1} {
			set ref ?
		} else {
			set ref [lindex $cur2 $refpos2]
		}
	} else {
		set ref [lindex $cur1 $refpos1]
		if {$ref eq "?"} {
			if {$refpos2 != -1} {set ref [lindex $cur2 $refpos2]}
		}
	}
	# merge alt
	if {$altpos1 == -1} {
		set alt1 [list_remove [list_remdup [list_sub $cur1 $alleleposs1]] $ref]
	} else {
		set alt1 [split [lindex $cur1 $altpos1] ,]
		if {[inlist ? $alt1]} {set alt1 ""}
	}
	if {$altpos2 == -1} {
		set alt2 [list_remove [list_remdup [list_sub $cur2 $alleleposs2]] $ref]
	} else {
		set alt2 [split [lindex $cur2 $altpos2] ,]
		if {[inlist ? $alt2]} {set alt2 ""}
	}
	# putsvars cur1 cur2 alt1 alt2
	if {$alt1 eq $alt2} {
		set alt $alt1
	} elseif {![llength $alt1]} {
		set alt $alt2
	} elseif {![llength $alt2]} {
		set alt $alt1
	} else {
		set alt [list_union $alt1 $alt2]
		set a1poss [list_cor $alt1 $alt]
		set a2poss [list_cor $alt2 $alt]
		foreach p $listfields1 {
			set v [lindex $cur1 $p]
			set v [split $v ,]
			set v [list_sub $v $a1poss]
			set v [list_change $v {{} -}]
			lset cur1 $p [join $v ,]
		}
		foreach p $listfields2 {
			set v [lindex $cur2 $p]
			set v [split $v ,]
			set v [list_sub $v $a2poss]
			set v [list_change $v {{} -}]
			lset cur2 $p [join $v ,]
		}
	}
	set remove [list $ref ? - N @]
	set type [lindex $region 3]
	if {$type eq "snp"} {
		lappend remove {}
	}
	set alt [list_lremove $alt $remove]
	set result $region
	lappend result $ref
	lappend result [join $alt ,]
	lappend result {*}[list_sub $cur1 $restposs1]
	lappend result $sequenced
	lappend result {*}[list_sub $cur2 $restposs2]
	lappend result {*}$merge
	return [join $result \t]
}

proc multicompar_getcomp {line poss split file prevcomp1} {
	set comp1 [list_sub $line $poss]
	if {![llength $line]} {return ~}
	if {$split && [regexp , [lindex $comp1 end]]} {
		error "split mode does not allow multiallelic variants: $comp1 in file $file"
	}
	set comp1 [join $comp1 " "]
	set prevcheck [loc_compare $prevcomp1 $comp1]
	if {$prevcheck >= 0} {
		if {$prevcheck > 0} {
			error "sorting error in \"$file\": \"$prevcomp1\" comes before \"$comp1\""
		} elseif {$prevcheck == 0 && !$split} {
			error "error in \"$file\": file uses split alleles (\"$prevcomp1\" occurs more than once and you are not running multicompar with the -split option)"
		}
	}
	return $comp1
}

proc multicompar {compar_file dir {split 0} {listfields {}}} {
	global cache joinposs1 joinposs2 comparposs1 mergeposs1 comparposs2 mergeposs2 dummy1 dummy2 restposs1 restposs2 refpos1 refpos2 altpos1 altpos2 alleleposs1 alleleposs2 listfields1 listfields2 sequenced2pos
	catch {close $f1}; catch {close $f2}; catch {close $o}
	set mergefields {xRef geneId mrnaAcc proteinAcc symbol orientation component componentIndex hasCodingRegion impact nucleotidePos proteinPos annotationRefSequence sampleSequence genomeRefSequence pfam}
	set allelefields {alleleSeq1 alleleSeq2}
	#
	if {[file isdir $dir]} {
		set name [file tail $dir]
		set file2 [gzfile $dir/fannotvar-$name.tsv]
	} else {
		set base [file root [file tail [gzroot $dir]]]
		set name [sourcename $base]
		set file2 [gzfile $dir]
	}
	if {![file exists $compar_file]} {
		file_write $compar_file [join {chromosome begin end type ref alt} \t]
	}
	set f1 [open $compar_file]
	set header1 [tsv_open $f1]
	if {[inlist $header1 alleleSeq1-$name]} {
		error "$name already present in $compar_file"
	}
	set comparposs1 [tsv_basicfields $header1 6 0]
	if {([lsearch $comparposs1 -1] != -1)} {
		puts stderr "header error in comparfile $compar_file"
		exit 1
	}
	set refpos1 [lindex $comparposs1 4]
	set altpos1 [lindex $comparposs1 5]
	if {$altpos1 == -1} {
		set alleleposs1 [list_find -glob $header1 alleleSeq*]
	}
	set joinposs1 [lrange $comparposs1 0 3]
	if {$split} {
		set comparposs1 [list_sub $comparposs1 {0 1 2 3 5}]
	} else {
		set comparposs1 $joinposs1
	}
	set tp1 [lindex $comparposs1 0]
	set dummy1 [list_fill [llength $header1] ?]
	set f2 [gzopen $file2]
	set header2 [tsv_open $f2]
	set comparposs2 [tsv_basicfields $header2 6 0]
	if {([lsearch [lrange $comparposs2 0 3] -1] != -1)} {
		puts stderr "header error in fannot_varfile2 $compar_file"
		exit 1
	}
	set refpos2 [lindex $comparposs2 4]
	set altpos2 [lindex $comparposs2 5]
	if {$altpos2 == -1} {
		set alleleposs2 [list_find -glob $header2 alleleSeq*]
	}
	set joinposs2 [lrange $comparposs2 0 3]
	if {$split} {
		set comparposs2 [list_sub $comparposs2 {0 1 2 3 5}]
	} else {
		set comparposs2 $joinposs2
	}
	set tp2 [lindex $comparposs2 0]
	set mergefields [list_common $mergefields $header2]
	set nonmergefields [list_lremove $header2 $mergefields]
	set mergeposs1 [list_cor $header1 $mergefields]
	set mergeposs2 [list_cor $header2 $mergefields]
	set dummy2 [list_fill [llength $header2] ?]
	set listfields1 [list_find -glob $header1 l:*]
	set listfields2 [list_find -glob $header2 l:*]
	foreach p $listfields {
		set poss [list_find -glob $header1 $p]
		if {[llength $poss]} {lappend listfields1 {*}$poss}
		set poss [list_find -glob $header2 $p]
		if {[llength $poss]} {lappend listfields2 {*}$poss}
	}
	set listfields1 [lsort -integer [list_remdup $listfields1]]
	set listfields2 [lsort -integer [list_remdup $listfields2]]
	# make output header
	set restfields1 [list_lremove [list_sub $header1 -exclude $comparposs1] $mergefields]
	set restfields1 [list_remove $restfields1 ref reference alt]
	set restposs1 [list_cor $header1 $restfields1]
	set oheader [list_concat {chromosome begin end type ref alt} $restfields1]
	set restfields2 [list_lremove [list_sub $header2 -exclude $comparposs2] $mergefields]
	set sequenced2pos [lsearch $header2 sequenced]
	set restfields2 [list_remove $restfields2 ref reference alt sequenced]
	set restposs2 [list_cor $header2 $restfields2]
	lappend oheader sequenced-$name
	foreach field $restfields2 {
		lappend oheader ${field}-$name
	}
	set oheader [list_concat $oheader $mergefields]
	# start
	set prevcomp1 {}
	set prevcomp2 {}
	set o [open $compar_file.temp w]
	puts $o [join $oheader \t]
	set cur1 [split [gets $f1] \t]
	if {[llength $cur1]} {lset cur1 $tp1 [chr_clip [lindex $cur1 $tp1]]}
	set comp1 [multicompar_getcomp $cur1 $comparposs1 $split $compar_file $prevcomp1]
	set cur2 [split [gets $f2] \t]
	if {[llength $cur2]} {lset cur2 $tp2 [chr_clip [lindex $cur2 $tp2]]}
	set comp2 [multicompar_getcomp $cur2 $comparposs2 $split $file2 $prevcomp2]
	set num 1; set next 100000
	while {![eof $f1] || ![eof $f2]} {
		incr num
		if {$num > $next} {putslog $num ; incr next 100000}
		set d [loc_compare $comp1 $comp2]
		if {$d == 0} {
			puts $o [multicompar_annot_join $cur1 $cur2]
			set cur1 [compare_annot_getline $f1]
			if {[llength $cur1]} {lset cur1 $tp1 [chr_clip [lindex $cur1 $tp1]]}
			set comp1 [multicompar_getcomp $cur1 $comparposs1 $split $compar_file $prevcomp1]
			set prevcomp1 $comp1
			set cur2 [compare_annot_getline $f2]
			if {[llength $cur2]} {lset cur2 $tp2 [chr_clip [lindex $cur2 $tp2]]}
			set comp2 [multicompar_getcomp $cur2 $comparposs2 $split $file2 $prevcomp2]
			set prevcomp2 $comp2
		} elseif {$d < 0} {
			while {$d < 0} {
				puts $o [multicompar_annot_join $cur1 -]
				if {[eof $f1]} break
				set cur1 [compare_annot_getline $f1]
				if {[llength $cur1]} {lset cur1 $tp1 [chr_clip [lindex $cur1 $tp1]]}
				set comp1 [multicompar_getcomp $cur1 $comparposs1 $split $compar_file $prevcomp1]
				set prevcomp1 $comp1
				if {![llength $cur1]} break
				incr num
				if {![expr {$num % 100000}]} {putslog $num}
				set d [loc_compare $comp1 $comp2]
			}
		} else {
			while {$d > 0} {
				puts $o [multicompar_annot_join - $cur2]
				if {[eof $f2]} break
				set cur2 [compare_annot_getline $f2]
				if {[llength $cur2]} {lset cur2 $tp2 [chr_clip [lindex $cur2 $tp2]]}
				set comp2 [multicompar_getcomp $cur2 $comparposs2 $split $file2 $prevcomp2]
				set prevcomp2 $comp2
				incr num
				if {![expr {$num % 100000}]} {putslog $num}
				set d [loc_compare $comp1 $comp2]
			}
		}
	}

	close $f1; close $f2; close $o
	catch {file rename -force $compar_file $compar_file.old}
	file rename -force $compar_file.temp $compar_file
}

proc multicompar_reannot_find {basedir sample args} {
	if {![llength $args]} {set args [list {}]}
	set sampledir [lindex [split $sample -] end]
	foreach pattern $args {
		set test [gzfile [file join $basedir $sample $pattern]]
		if {[file exists $test]} {
			return $test
		}
	}
	foreach pattern $args {
		set test [gzfile [file join [file dir $basedir] $sample $pattern]]
		if {[file exists $test]} {
			return $test
		}
	}
	foreach pattern $args {
		set test [gzfile [file join $basedir $sampledir $pattern]]
		if {[file exists $test]} {
			return $test
		}
	}
	foreach pattern $args {
		set test [gzfile [file join [file dir $basedir] $sampledir $pattern]]
		if {[file exists $test]} {
			return $test
		}
	}
	foreach pattern $args {
		set test [gzfile [file join $basedir $pattern]]
		if {[file exists $test]} {
			return $test
		}
	}
	return {}
}

proc multicompar_reannot {compar_file {force 0} {regonly 0} {skipincomplete 0}} {

	set compar_file [file_absolute $compar_file]
	set basedir [file dir $compar_file]
	catch {close $f}; catch {close $o}
	set f [gzopen $compar_file]
	set header [tsv_open $f]
	set pos -1
	unset -nocomplain samplea
	set samples {}
	foreach field $header {
		incr pos
		set temp [split $field -]
		if {[llength $temp] > 1} {
			set sample [join [lrange $temp 1 end] -]
			lappend samplea(poss,$sample) $pos
			lappend samplea(fields,$sample) [lindex $temp 0]
			list_addnew samples $sample
		}
	}
	set referencepos [lsearch $header reference]
	if {$referencepos == -1} {set referencepos [lsearch $header ref]}
	putslog "file contains samples: $samples"
	foreach sample $samples {
		set samplea(a1,$sample) [lsearch $header alleleSeq1-$sample]
		set samplea(a2,$sample) [lsearch $header alleleSeq2-$sample]
		set samplea(rpos,$sample) [lsearch $header refscore-$sample]
		set samplea(cpos,$sample) [lsearch $header coverage-$sample]
		set samplea(seq,$sample) [lsearch $header sequenced-$sample]
		set samplea(dir,$sample) [multicompar_reannot_find $basedir $sample]
		set samplea(regionfile,$sample) [multicompar_reannot_find $basedir $sample sreg-$sample.tsv]
		if {[file exists $samplea(dir,$sample)/allpos]} {
			set samplea(type,$sample) rtg
			annot_rtg_init $samplea(dir,$sample)
			set samplea(rtgposs,$sample) {}
			foreach field {
				alleleSeq1 alleleSeq2 posterior coverage correction
				numA numC numG numT percA percC percG percT nonidentityposterior
			} {
				lappend samplea(rtgposs,$sample) [lsearch $header ${field}-$sample]
			}
		} elseif {[file exists $samplea(regionfile,$sample)]} {
			set samplea(type,$sample) cg
			if {!$regonly} {
				putslog "Using $samplea(dir,$sample) to reannot"
				annot_coverage_init $samplea(dir,$sample) $sample
			}
			putslog "Using $samplea(regionfile,$sample) to reannot"
			annot_region_init $samplea(regionfile,$sample)
		} else {
			if {!$skipincomplete} {
				error "no sorted region file (sreg-$sample.tsv) or allpos dir (for rtg) found: not properly processed sample"
			}
			set samples [list_remove $samples $sample]
		}
		set samplea(varall,$sample) [multicompar_reannot_find $basedir $sample varall-$sample.tsv]
		if {$samplea(varall,$sample) ne ""} {
			if {!$regonly} {
				putslog "Using $samplea(varall,$sample) to reannot"
				annot_varall_init $samplea(varall,$sample) $sample $header
			}
		}
		set samplea(todo,$sample) {}
		foreach temp {refcons nocall cluster} {
			set regfile [multicompar_reannot_find $basedir $sample reg_${temp}-$sample.tsv]
			if {$regfile ne "" && [inlist $samplea(fields,$sample) $temp]} {
				lappend samplea(todo,$sample) [list [lsearch $header ${temp}-$sample] 1 $regfile]
			}
		}
		list_foreach {field value regfile} $samplea(todo,$sample) {
			annot_region_init $regfile
		}
	}
	# start processing
	set o [open $compar_file.temp w]
	puts $o [join $header \t]
	set poss [list_cor $header {chromosome begin end}]
	set num 0
	while {![eof $f]} {
		incr num
		if {![expr $num%10000]} {putslog $num}
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set reference [lindex $line $referencepos]
		set loc [list_sub $line $poss]
		foreach {chr begin end} $loc break
		foreach sample $samples {
			#if {!$force && ([lindex $line $samplea(a1,$sample)] ne "?")} continue
			list_foreach {field value regfile} $samplea(todo,$sample) {
				if {[lindex $line $field] == "-"} continue
				if {!$force && [lindex $line $field] != "?"} continue
				set r [annot_region_get $regfile $chr $begin $end]
				if {$r} {lset line $field $value} else {lset line $field {}}
			}
			if {$samplea(varall,$sample) ne ""} {
				if {!$regonly} {
					annot_varall_annot $samplea(varall,$sample) $sample $loc $force line
				}
			}
			if {$samplea(type,$sample) eq "cg"} {
				if {!$regonly && ($force || ([lindex $line $samplea(rpos,$sample)] eq "?") || ([lindex $line $samplea(cpos,$sample)] eq "?"))} {
					foreach {r c} [annot_coverage_get $samplea(dir,$sample) $sample $chr $begin] break
					if {$samplea(rpos,$sample) != -1} {lset line $samplea(rpos,$sample) $r}
					lset line $samplea(cpos,$sample) $c
				}
				set seq [lindex $line $samplea(seq,$sample)]
				if {$seq eq "?"} {
					set r [annot_region_in $samplea(regionfile,$sample) $chr $begin $end]
					set a1 [inlist [list - ?] [lindex $line $samplea(a1,$sample)]]
					set a2 [inlist [list - ?] [lindex $line $samplea(a2,$sample)]]
					if {$r} {
						lset line $samplea(seq,$sample) r
						if {$a1} {lset line $samplea(a1,$sample) $reference}
						if {$a2} {lset line $samplea(a2,$sample) $reference}
					} else {
						lset line $samplea(seq,$sample) u
						if {$a1} {lset line $samplea(a1,$sample) -}
						if {$a2} {lset line $samplea(a2,$sample) -}
					}
				}
			} else {
				set sub [list_sub $line $samplea(rtgposs,$sample)]
				if {!$force && ![inlist [list_sub $line $samplea(rtgposs,$sample)] ?]} continue
				set temp [annot_rtg_get $samplea(dir,$sample) $chr $begin]
				set rtgdata [lrange $temp 4 end-1]
				if {[llength $rtgdata] < [llength $samplea(rtgposs,$sample)]} {
					lset line $samplea(seq,$sample) u
					foreach pos $samplea(rtgposs,$sample) v $rtgdata {
						if {$pos == -1} continue
						lset line $pos -
					}
				} elseif {[llength $rtgdata] > [llength $samplea(rtgposs,$sample)]} {
					error "error in rtg data\nsample=$sample\nline=$line"
				} else {
					set coverage [lindex $rtgdata 3]
					if {$coverage < 10} {
						lset line $samplea(seq,$sample) u
					} else {
						lset line $samplea(seq,$sample) r
					}
					foreach pos $samplea(rtgposs,$sample) v $rtgdata {
						if {$pos == -1} continue
						lset line $pos $v
					}
				}
			}
		}
		puts $o [join $line \t]
	}

	close $o
	close $f
	foreach sample $samples {
		list_foreach {field value regfile} $samplea(todo,$sample) {
			annot_region_close $regfile
		}
		if {$samplea(varall,$sample) ne ""} {
			if {!$regonly} {
				annot_varall_close $samplea(varall,$sample) $sample
			}
		}
		if {$samplea(type,$sample) eq "cg"} {
			if {!$regonly} {
				annot_coverage_close $samplea(dir,$sample) $sample
			}
			annot_region_close $samplea(regionfile,$sample)
		} else {
			annot_rtg_close $samplea(dir,$sample)
		}
	}
	file rename -force $compar_file $compar_file.old
	file rename -force $compar_file.temp $compar_file
}

proc cg_multicompar.old {args} {
	set reannot 0
	set regonly 0
	set split 0
	set listfields {}
	set pos 0
	while 1 {
		set key [lindex $args $pos]
		switch -- $key {
			-reannot {
				putslog "Also reannot"
				set reannot 1
			}
			-reannotregonly {
				putslog "Also reannot"
				set reannot 1
				set regonly 1
			}
			-split {
				incr pos
				set split [true [lindex $args $pos]]
			}
			-listfields {
				incr pos
				set listfields [lindex $args $pos]
			}
			-- break
			default {
				break
			}
		}
		incr pos
	}
	set args [lrange $args $pos end]
	if {([llength $args] < 1)} {
		puts "Wrong number of arguments"
		errorformat multicompar
		exit 1
	}
	foreach {compar_file} $args break
	set dirs [lrange $args 1 end]
	foreach dir $dirs {
		set dir [file_absolute $dir]
		if {![file isdir $dir] && [llength [cg select -n $dir]]} {
			set header [cg select -h $dir]
			set samples [samples $header]
			set basicfields [list_sub $header -exclude [list_find -regexp $header -]]
			putslog "$dir is already a multicompar, merging"
			foreach sample $samples {
				putslog "Adding $dir sample $sample"
				catch {file delete $dir.temp}
				file mkdir $dir.temp
				set fields $basicfields
				foreach field [list_sub $header [list_find -regexp $header -$sample]] {
					set nfield [join [lrange [split $field -] 0 end-1] -]
					lappend fields $nfield=\$$field
				}
				cg select -f $fields $dir $dir.temp/vars-$sample.tsv
				multicompar $compar_file $dir.temp/vars-$sample.tsv $split $listfields
				file delete $dir.temp/vars-$sample.tsv
			}
			file delete -force $dir.temp
		} else {
			putslog "Adding $dir"
			multicompar $compar_file $dir $split $listfields
		}
	}
	if {$reannot} {
		putslog "Reannotating $compar_file"
		multicompar_reannot $compar_file 0 $regonly
	}
}

proc cg_multicompar_reannot {args} {
	if {[llength $args] < 1} {
		errorformat multicompar_reannot
		exit 1
	}
	set compar_file [list_shift args]
	set force 0
	set regonly 0
	set skipincomplete 0
	foreach option $args {
		switch $option {
			force {set force 1}
			regonly {set regonly 1}
			skipincomplete {set skipincomplete 1}
			default {error "unrecognized option $option"}
		}
	}
	putslog "Reannotating $compar_file"
	multicompar_reannot $compar_file $force $regonly $skipincomplete
}
