#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

target cg_chromosomes {cg_chromosomes\((.*)\)} -deps {$target1/oricg/ASM/REF/coverage*} -code {
	set result {}
	set files [gzfiles $dep]
	foreach file $files {
		lappend result [chr_clip [lindex [split $file -] 1]]
	}	
	set ::cg_chromosomes($target1) $result
}

target cg_svar {(.*)/svar-([^/-]*).tsv} -deps {$target1/oricg/ASM/var-*-ASM*.tsv} -code {
	set varfile [gzfiles $dep]
	if {[llength $varfile] != 1} {error "could not identify varfile"}
	set varfile [lindex $varfile 0]
	set f [gzopen $varfile]
	set info {}
	while {![eof $f]} {
		set line [gets $f]
		if {[string index $line 0] ne "#"} break
		lappend info $line
	}
	catch {close $f}
	if {[file exists $target1/info.txt]} {
		set test [split [file_read $target1/info.txt] \n]
		if {$info ne $test} {
			error "$target1 already has info.txt that differs from data in the source $dir"
		}
	}
	file_write $target1/info.txt [join $info \n]
	putslog "Sort var file ($varfile)"
	cg select -s "chromosome begin end varType" $varfile $target.temp
	file rename -force $target.temp $target
}

target cg_sgene {(.*)/sgene-([^/-]*).tsv} -deps {$target1/oricg/ASM/gene-*-ASM*.tsv*} -code {
	set genefile [gzfile $dep]
	set genefile [list_sub $genefile -exclude [list_find -regexp $genefile gene-var-summary]]
	if {[llength $genefile] != 1} {error "could not identify genefile"}
	putslog "Sort gene file ($genefile)"
	cg select -s "chromosome begin end" $genefile $target.temp
	file rename -force $target.temp $target
}

# annotated vars file
target cg_annotvar {(.*)/annotvar-([^/-]*).tsv} -deps {$target1/svar-$target2.tsv $target1/sgene-$target2.tsv} -code {
	foreach {svarfile sgenefile} $deps break
	putslog "Create annotated varfile $target"
	cg var2annot $svarfile $sgenefile $target.temp
	file rename -force $target.temp $target
}

# if not from cg, we do not have svar and sgene, take from first var_* file found
target cg_annotvar_other {(.*)/annotvar-([^/-]*).tsv} -deps {$target1/var_*.tsv} -code {
	set file [gzfile $dep]
	set ext [file extension $file]
	if {$ext ne ".tsv"} {
		mklink $file $target$ext
	} else {
		mklink $file $target
	}
}

# if we also do not find a var_* file, take from first variant* file found
target cg_annotvar_other {(.*)/annotvar-([^/-]*).tsv} -deps {$target1/variant*.tsv} -code {
	set file [gzfile $dep]
	set ext [file extension $file]
	if {$ext ne ".tsv"} {
		mklink $file $target$ext
	} else {
		mklink $file $target
	}
}

# if not from cg, we do not have svar and sgene, take from first var_file found
target cg_annotvar_other {(.*)/annotvar-([^/-]*).tsv} -deps {$target1/var_*.tsv} -code {
	mklink $dep $target
}

target cg_sreg {(.*)/sreg-([^/-]*).tsv} -deps {$target1/oricg/ASM/reg-*-ASM*.tsv*} -code {
	set regfile [gzfile $dep]
	putslog "Sort region file ($regfile)"
	cg select -s "chromosome begin end" $regfile $target.temp
	file rename -force $target.temp $target
}

target cg_regfromsvar {(.*)/sreg-([^/-]*).tsv} -deps {$target1/svar-$target2.tsv} -code {
	foreach svarfile $deps break
	putslog "Extract $target from $svarfile"
	cg select -q {$varType != "no-call" && $varType != "no-ref"} -f "chromosome begin end" $svarfile $target.temp
	cg regjoin $target.temp > $target.temp2
	file rename -force $target.temp2 $target
	file delete $target.temp
}

target cg_reg_refcons {(.*)/reg_refcons-([^/-]*).tsv} -deps {$target1/svar-$target2.tsv} -code {
	putslog "Find refcons regions for $dep"
	cg refconsregions $dep > $target.temp
	file rename -force $target.temp $target
}

target cg_reg_nocall {(.*)/reg_nocall-([^/-]*).tsv} -deps {$target1/svar-$target2.tsv} -code {
	putslog "Find partial no-call regions for dep"
	if {[catch {
		nocallregions $dep $target.temp
	}]} {
		puts stderr "Could not make $target (old version files ?)"
	} else {
		file rename -force $target.temp $target
	}
}

target cg_coverage {(.*)/coverage/coverage-([^/-]*)-([^/-]*).bcol} -deps {$target1/oricg/ASM/REF/coverage*-chr$target3-*} -code {
	# make coverage files
	set file [gzfile $dep]
	file mkdir $target1/coverage
	set sample $target2
	set chr $target3
	set chr [chr_clip $chr]
	set header [cg select -h $file]
	foreach posfield {offset pos} {
		if {[lsearch $header $posfield] != -1}  break
	}
	if {$posfield == -1} {
		exiterror "No position/offset field found in $file"
	}
	foreach covfield {uniqueSequenceCoverage coverage} {
		if {[lsearch $header $covfield] != -1}  break
	}
	if {$covfield == -1} {
		exiterror "No coverage/uniqueSequenceCoverage field found in $file"
	}
	set other [list_remove $header $posfield $covfield]
	foreach field $other {
		set base $target1/coverage/$field-$sample-$chr
		if {![file exists $base.bcol]} {
			putslog "Making $base.bcol"
			if {[catch {
				exec [catprog $file] $file | cg bcol make -p $posfield -t s -n -1 $base $field <  $file
			} e]} {
				exec [catprog $file] $file | cg bcol make -p $posfield -t i -n -1 $base $field <  $file
			}
		}
	}
	set base [file root $target]
	if {![file exists $base.bcol]} {
		putslog "Making $base.bcol"
		exec [catprog $file] $file | cg bcol make -p $posfield -t su $base $covfield
	}
}

target cg_cpCNV {(.*)/CNV} -deps {$target1/oricg/ASM/CNV} -code {
	putslog "Copying CNV"
	file delete -force $target.temp
	catch {
		file copy $dep $target.temp
		file rename $target.temp $target
	} result
	putslog $result
}

target cg_cpSV {(.*)/SV} -deps {$target1/oricg/ASM/SV} -code {
	putslog "Copying SV"
	file delete -force $target.temp
	catch {
		file copy $dep $target.temp
		file rename $target.temp $target
	} result
	putslog $result
}

target cg_cgsv {(.*)/cgsv-(.*).tsv} -deps {$target1/SV/allJunctionsBeta-*.tsv*} -code {
	cg convcgsv [gzfile $dep] $target
}

target cg_cgsv {(.*)/cgsv-(.*).tsv} -deps {$target1/SV/annotatedJunctionsAlpha-*.tsv*} -code {
	cg convcgsv [gzfile $dep] $target
}

target cg_cgcnv {(.*)/cgcnv-(.*).tsv} -deps {$target1/CNV/cnvSegmentsBeta-*.tsv*} -code {
	cg convcgcnv [gzfile $dep] $target
}

target cg_cgcnv {(.*)/cgcnv-(.*).tsv} -deps {$target1/CNV/cnvSegmentsAlpha-*.tsv*} -code {
	cg convcgcnv [gzfile $dep] $target
}

target reg_cluster {(.*)/reg_cluster-(.*).tsv} -deps {$target1/annotvar-$target2.tsv} -code {
	cg clusterregions < $dep > $target.temp
	file rename -force $target.temp $target
}

target reg_ns {(.*)/reg_ns-(.*).tsv} -deps {$target1/annotvar-$target2.tsv} -code {
	putslog "Find regions with N's for $dep"
	cg select -f {chromosome begin end} -q {$alleleSeq1 ~ /[N?]/ || $alleleSeq2 ~ /[N?]/} < $dep > $target.temp
	file rename -force $target.temp $target
}

target reg_lowscore {(.*)/reg_lowscore-(.*).tsv} -deps {$target1/annotvar-$target2.tsv} -code {
	set header [cg select -h $dep]
	if {[llength [list_common $header {totalScore1 totalScore2}]] == 2} {
		putslog "Find regions with lowscores for $dep"
		cg select -f {chromosome begin end} -q {$totalScore1 < 60 || $totalScore2 < 60} < $dep > $target.temp
		file rename -force $target.temp $target
	}
}

target cg_fannotvar {(.*)/fannotvar-(.*).tsv} -deps {$target1/annotvar-$target2.tsv ($target1/reg_refcons-$target2.tsv) ($target1/reg_cluster-$target2.tsv)} -code {
	# add filterdata to annotvar
	set todo {}
	if {[file exists $target1/reg_refcons-$target2.tsv]} {
		lappend todo [list refcons rc $target1/reg_refcons-$target2.tsv]
	}
	# lappend todo [list nocall nc reg_nocall-$target2.tsv]
	if {[file exists $target1/reg_refcons-$target2.tsv]} {
		lappend todo [list cluster cl $target1/reg_cluster-$target2.tsv]
	}
	# lappend todo [list trf trf $dbdir/regdb-simple_repeats.tsv]
	# lappend todo [list str str $dbdir/regdb-microsatelite.tsv]
	# lappend todo [list segdup sd $dbdir/regdb-segdups.tsv]
	# lappend todo [list selfchain sc $dbdir/regdb-selfchain.tsv]
	# lappend todo [list repeat rp $dbdir/regdb-repeatmasker.tsv]
	# lappend todo [list rna rna $dbdir/regdb-rnagenes.tsv]
	# lappend todo [list more5pct m5 $dbdir/regdb-1000genomesmore5pct.tsv]
	# lappend todo [list more1pct m1 $dbdir/regdb-1000genomesmore1pct.tsv]
	# foreach file [lsort -dictionary [glob -nocomplain $dbdir/checked*.tsv]] {
	# 	set value [lindex [split [file root [file tail $file]] _] 1]
	# 	if {$value eq ""} {set value checked}
	# 	lappend todo [list checked $value $file]
	# }
	annot_annotvar $dep $target $todo $target1
}

target reg_covered {(.*)/reg-(.*).covered} -deps {$target1/sreg-$target2.tsv} -code {
	putslog "Genomic coverage of sequenced regions"
	cg covered $dep > $target.temp
	file rename -force $target.temp $target
}

target cg_filteredrefcons {(.*)/filteredrefcons-(.*).covered} -deps {$target1/sreg-$target2.tsv $target1/reg_refcons-$target2.tsv} -code {
	putslog "Coverage of refcons region"
	cg regsubtract $target1/sreg-$target2.tsv $target1/reg_refcons-$target2.tsv > $target1/filteredrefcons-$target2.tsv.temp
	file rename -force $target1/filteredrefcons-$target2.tsv.temp $target1/filteredrefcons-$target2.tsv
	cg covered $target1/filteredrefcons-$target2.tsv > $target.temp
	file rename -force $target.temp $target
}

target cg_filteredns {(.*)/filteredns-(.*).tsv} -deps {$target1/sreg-$target2.tsv $target1/reg_ns-$target2.tsv} -code {
	putslog "Coverage of ns region"
	cg regsubtract $target1/sreg-$target2.tsv $target1/reg_ns-$target2.tsv > $target.temp
	file rename -force $target.temp $target
}

target cg_filteredns_covered {(.*)/filteredns-(.*).covered} -deps {$target1/filteredns-$target2.tsv} -code {
	putslog "Making filteredns-$target2.covered"
	cg covered $dep > $target.temp
	file rename -force $target.temp $target
}

target cg_filteredlowscore {(.*)/filteredlowscore-(.*).tsv} -deps {$target1/sreg-$target2.tsv $target1/reg_lowscore-$target2.tsv} -code {
	cg regsubtract $target1/sreg-$target2.tsv $target1/reg_lowscore-$target2.tsv > $target.temp
	file rename -force $target.temp $target
}

target cg_filteredlowscore_covered {(.*)/filteredlowscore-(.*).covered} -deps {$target1/filteredlowscore-$target2.tsv} -code {
	cg covered $dep > $target.temp
	file rename -force $target.temp $target
}

target cg_refcons_histo {(.*)/histo-refcons-(.*).tsv} -deps {$target1/reg_refcons-$target2.tsv} -code {
	putslog "Making $target"
	cg reghisto $dep > $target.temp
	file rename -force $target.temp $target
}

target cg_filteredcluster {(.*)/filteredcluster-(.*).tsv} -deps {$target1/sreg-$target2.tsv $target1/reg_cluster-$target2.tsv} -code {
	putslog "Coverage of clusters region"
	cg regsubtract $target1/sreg-$target2.tsv $target1/reg_cluster-$target2.tsv > $target.temp
	file rename -force $target.temp $target
}

target cg_filteredcluster_covered {(.*)/filteredcluster-(.*).covered} -deps {$target1/filteredcluster-$target2.tsv} -code {
	putslog "Making $target"
	cg covered $dep > $target.temp
	file rename -force $target.temp $target
}

target cg_process_sample {(.*)/cg_process_sample-(.*).finished} -deps {
	$target1/annotvar-$target2.tsv
	$target1/sreg-$target2.tsv
	$target1/reg_refcons-$target2.tsv
	$target1/coverage/coverage-$target2-${_cg_chromosomes($target1)}.bcol
	$target1/fannotvar-$target2.tsv
	($target1/reg_nocall-$target2.tsv)
	($target1/SV) ($target1/CNV)
	($target1/cgsv-$target2.tsv)
	($target1/cgcnv-$target2.tsv)
	($target1/reg_cluster-$target2.tsv)
	($target1/reg_ns-$target2.tsv)
	($target1/reg_lowscore-$target2.tsv)
	($target1/reg-$target2.covered)
	($target1/filteredrefcons-$target2.covered)
	($target1/filteredns-$target2.covered)
	($target1/filteredlowscore-$target2.covered)
	($target1/filteredcluster-$target2.covered)
	($target1/histo-refcons-$target2.tsv)
} -code {
	file_write $target [timestamp]
}

# Added for autoloading lib with cgmake_lib
proc cgmakelib_process_sample {args} {
	return targetdir
}

