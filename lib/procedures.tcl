
proc cg {args} {
	# puts "cg $args"
	eval exec cg $args 2>@ stderr
}

proc process_sample {dir dbdir {force 0}} {
	set keepdir [pwd]
	cd $dir
	puts stderr "Processing sample $dir"
	set name [file tail $dir]
	set varfile [lindex [glob ASM/var-*-ASM.tsv] 0]
	set genefile [lindex [glob ASM/gene-*-ASM.tsv] 0]
	set regfile [lindex [glob ASM/reg-*-ASM.tsv] 0]
	# sort files
	if {$force || ![file exists svar-$name.tsv]} {
		puts stderr "Sort var file ($varfile)"
		cg select -s "chromosome begin end varType" < $varfile > temp.tsv
		file rename temp.tsv svar-$name.tsv
	}
	if {$force || ![file exists sgene-$name.tsv]} {
		puts stderr "Sort gene file ($genefile)"
		cg select -s "chromosome begin end varType" < $genefile > temp.tsv
		file rename temp.tsv sgene-$name.tsv
	}
	if {$force || ![file exists sreg-$name.tsv]} {
		puts stderr "Sort region file ($regfile)"
		cg select -s "chromosome begin end" < $regfile > temp.tsv
		file rename temp.tsv sreg-$name.tsv
	}
	# annotated vars file
	if {$force || ![file exists annotvar-$name.tsv]} {
		puts stderr "Create annotated varfile annotvar-$name.tsv"
		cg var2annot svar-$name.tsv sgene-$name.tsv temp.tsv
		file rename temp.tsv annotvar-$name.tsv
	}
	# sample specific filters
	if {$force || ![file exists reg_refcons-$name.tsv]} {
		puts stderr "Find refcons regions for var-$name.tsv"
		cg refconsregions svar-$name.tsv > temp.tsv
		file rename temp.tsv reg_refcons-$name.tsv
	}
	if {$force || ![file exists reg_cluster-$name.tsv]} {
		puts stderr "Find cluster regions for svar-$name.tsv"
		cg clusterregions < annotvar-$name.tsv > temp.tsv
		file rename temp.tsv reg_cluster-$name.tsv
	}
	if {$force || ![file exists reg_ns-$name.tsv]} {
		puts stderr "Find regions with N's for svar-$name.tsv"
		cg select -f {chromosome begin end} -q {$alleleSeq1 ~ /[N?]/ || $alleleSeq2 ~ /[N?]/} < annotvar-$name.tsv > temp.tsv
		file rename temp.tsv reg_ns-$name.tsv
	}
	if {$force || ![file exists reg_lowscore-$name.tsv]} {
		puts stderr "Find regions with lowscores for svar-$name.tsv"
		cg select -f {chromosome begin end} -q {$totalScore1 < 60 || $totalScore2 < 60} < annotvar-$name.tsv > temp.tsv
		file rename temp.tsv reg_lowscore-$name.tsv
	}
	if {$force || ![file exists fannotvar-$name.tsv]} {
		# add filterdata to annotvar
		set todo {}
		lappend todo [list refcons rc $dir/reg_refcons-$name.tsv]
		lappend todo [list cluster cl $dir/reg_cluster-$name.tsv]
		lappend todo [list trf trf $dbdir/regdb-simple_repeats.tsv]
		lappend todo [list str str $dbdir/regdb-microsatelite.tsv]
		lappend todo [list segdup sd $dbdir/regdb-segdups.tsv]
		lappend todo [list selfchain sc $dbdir/regdb-selfchain.tsv]
		lappend todo [list repeat rp $dbdir/regdb-repeatmasker.tsv]
		lappend todo [list rna rna $dbdir/regdb-rnagenes.tsv]
		foreach file [lsort -dictionary [glob $dbdir/checked*.tsv]] {
			set value [lindex [split [file root [file tail $file]] _] 1]
			if {$value eq ""} {set value checked}
			lappend todo [list checked $value $file]
		}
		annot_annotvar annotvar-$name.tsv fannotvar-$name.tsv $todo $dir
	}
	# coverage
	if {$force || ![file exists reg-$name.covered]} {
		puts stderr "Coverage of sequenced regions"
		cg covered sreg-$name.tsv > temp.tsv
		file rename temp.tsv reg-$name.covered
	}
	if {$force || ![file exists filteredrefcons-$name.covered]} {
		puts stderr "Coverage of refcons region"
		cg regsubtract sreg-$name.tsv reg_refcons-$name.tsv > temp.tsv
		file rename temp.tsv filteredrefcons-$name.tsv
		cg covered filteredrefcons-$name.tsv > temp.tsv
		file rename temp.tsv filteredrefcons-$name.covered
	}
	if {$force || ![file exists filteredns-$name.covered]} {
		puts stderr "Coverage of ns region"
		cg regsubtract sreg-$name.tsv reg_ns-$name.tsv > temp.tsv
		file rename temp.tsv filteredns-$name.tsv
		cg covered filteredns-$name.tsv > temp.tsv
		file rename temp.tsv filteredns-$name.covered
	}
	if {$force || ![file exists filteredlowscore-$name.covered]} {
		puts stderr "Coverage of lowscore region"
		cg regsubtract sreg-$name.tsv reg_lowscore-$name.tsv > temp.tsv
		file rename temp.tsv filteredlowscore-$name.tsv
		cg covered filteredlowscore-$name.tsv > temp.tsv
		file rename temp.tsv filteredlowscore-$name.covered
	}
	if {$force || ![file exists histo-refcons-$name.tsv]} {
		cg reghisto reg_refcons-$name.tsv > temp.tsv
		file rename temp.tsv histo-refcons-$name.tsv
	}
	if {$force || ![file exists filteredcluster-$name.covered]} {
		puts stderr "Coverage of clusters region"
		cg regsubtract sreg-$name.tsv reg_cluster-$name.tsv > temp.tsv
		file rename temp.tsv filteredcluster-$name.tsv
		cg covered filteredcluster-$name.tsv > temp.tsv
		file rename temp.tsv filteredcluster-$name.covered
	}
	cd $keepdir
}

proc process_compare_checkfield {file field} {
	set f [open $file]
	set header [tsv_open $f]
	close $f
	inlist $header $field
}

proc process_compare {dir1 dir2 dbdir resultsdir {force 0}} {
	set dir1 [file normalize $dir1]
	set dir2 [file normalize $dir2]
	set keepdir [pwd]
	file mkdir $resultsdir
	cd $resultsdir
	set name1 [file tail $dir1]
	set name2 [file tail $dir2]
	# compare
	if {$force || (![file exists compar_${name1}_${name2}.tsv] && ![file exists fcompar_${name1}_${name2}.tsv])} {
		puts stderr "comparing $name1 and $name2 (-> $resultsdir)"
		cg compare_annot $name1 $dir1/fannotvar-$name1.tsv $dir1/sreg-$name1.tsv $name2 $dir2/fannotvar-$name2.tsv $dir2/sreg-$name2.tsv compar_${name1}_${name2}.tsv
	}
	if {$force || ![file exists fcompar_${name1}_${name2}.tsv]} {
		reannot_compare compar_${name1}_${name2}.tsv $dir1 $dir2 fcompar_${name1}_${name2}.tsv
	}
	if {$force || ![file exists pvtcompar_${name1}_${name2}.tsv]} {
		cg compare_pvt < fcompar_${name1}_${name2}.tsv > temp.tsv
		file rename temp.tsv pvtcompar_${name1}_${name2}.tsv
		cg compare_pvtsummary < pvtcompar_${name1}_${name2}.tsv > temp.tsv
		file rename temp.tsv summarycompar_${name1}_${name2}.tsv
	}
	cd $keepdir
}

proc process_compare_coverage {dir1 dir2 dbdir resultsdir {force 0}} {
	cd $resultsdir
	# region files
	set dir1 [file normalize $dir1]
	set dir2 [file normalize $dir2]
	set name1 [file tail $dir1]
	set name2 [file tail $dir2]
	set target ${name1}_${name2}_regcommon.tsv.covered
	if {$force || ![file exists $target]} {
		cg regsubtract $dir1/sreg-${name1}.tsv $dir2/sreg-${name2}.tsv > temp.tsv
		cg regsubtract $dir1/sreg-${name1}.tsv temp.tsv > ${name1}_${name2}_regcommon.tsv
		cg covered ${name1}_${name2}_regcommon.tsv > $target
	}
	# region from each sample
	set todo {}
	lappend todo [list refcons $dir1/reg_refcons-${name1}.tsv $dir2/reg_refcons-${name2}.tsv]
	lappend todo [list ns $dir1/reg_ns-${name1}.tsv $dir2/reg_ns-${name2}.tsv]
	lappend todo [list lowscore $dir1/reg_lowscore-${name1}.tsv $dir2/reg_lowscore-${name2}.tsv]
	lappend todo [list cluster $dir1/reg_cluster-${name1}.tsv $dir2/reg_cluster-${name2}.tsv]
	lappend todo [list trf $dbdir/regdb-simple_repeats.tsv {}]
	lappend todo [list str $dbdir/regdb-microsatelite.tsv {}]
	lappend todo [list b20 $dir1/ASM/REF/below20coverage.regions $dir2/ASM/REF/below20coverage.regions]
	lappend todo [list a100 $dir1/ASM/REF/above100coverage.regions $dir2/ASM/REF/above100coverage.regions]
	lappend todo [list segdup $dbdir/regdb-segdups.tsv {}]
	lappend todo [list selfchain $dbdir/regdb-selfchain.tsv {}]
	lappend todo [list repeat $dbdir/regdb-repeatmasker.tsv {}]
	list_foreach {dbname rfile1 rfile2} $todo {
		set target ${name1}_${name2}_regcommon-$dbname.tsv.covered
		if {$force || ![file exists $target]} {
			puts "Making $target"
			cg regsubtract ${name1}_${name2}_regcommon.tsv $rfile1 > temp.tsv
			if {$rfile2 ne ""} {
				file rename -force temp.tsv temp-old.tsv
				cg regsubtract temp-old.tsv $rfile2 > temp.tsv
			}
			file rename -force temp.tsv ${name1}_${name2}_regcommon-$dbname.tsv
			cg covered ${name1}_${name2}_regcommon-$dbname.tsv > temp.covered
			file rename -force temp.covered $target
		} else {
			puts "$target already there"
		}
	}
	# subsequent subtractions (filters)
	set newname ${name1}_${name2}_regcommon-rc
	set todo [lrange $todo 1 end]
	list_foreach {dbname rfile1 rfile2} $todo {
		set oldname $newname
		set newname $oldname-$dbname
		set target $newname.tsv.covered
		if {$force || ![file exists $target]} {
			puts "Making $target"
			cg regsubtract $oldname.tsv  $rfile1 > temp.tsv
			if {$rfile2 ne ""} {
				file rename -force temp.tsv temp-old.tsv
				cg regsubtract temp-old.tsv $rfile2 > temp.tsv
			}
			file rename -force temp.tsv $newname.tsv
			cg covered $newname.tsv > temp.covered
			file rename -force temp.covered $target
		} else {
			puts "$target already there"
		}
	}
}

if 0 {

	lappend auto_path /complgen/bin/complgen/apps/cg/lib
	lappend auto_path ~/dev/completegenomics/lib
	package require Tclx
	signal -restart error SIGINT
	package require Extral

	set basedir /media/passport/complgen
	set basedir /complgen
	set dir1 $basedir/GS102
	set dir2 $basedir/GS103
	set dbdir /complgen/refseq
	set resultsdir $basedir/compar_GS102_GS103
	set force 0
	set name1 [file tail $dir1]
	set name2 [file tail $dir2]
	set file1 $dir1/fannotvar-$name1.tsv
	set file2 $dir2/fannotvar-$name2.tsv
	set regfile1 $dir1/sreg-$name1.tsv
	set regfile2 $dir2/sreg-$name2.tsv
	set outfile $resultsdir/compar_${name1}_${name2}.tsv
	file mkdir $resultsdir
	cd $resultsdir

	annot_compare_coverage compar_GS102_GS103.tsv $dir1 $dir2 ftemp.tsv

}
