
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
		cg select -s "chromosome begin end" < $varfile > svar-$name.tsv
	}
	if {$force || ![file exists sgene-$name.tsv]} {
		puts stderr "Sort gene file ($genefile)"
		cg select -s "chromosome begin end" < $genefile > sgene-$name.tsv
	}
	if {$force || ![file exists sreg-$name.tsv]} {
		puts stderr "Sort region file ($regfile)"
		cg select -s "chromosome begin end" < $regfile > sreg-$name.tsv
	}
	# annotated vars file
	if {$force || ![file exists annotvar-$name.tsv]} {
		puts stderr "Create annotated varfile annotvar-$name.tsv"
		cg var2annot svar-$name.tsv sgene-$name.tsv annotvar-$name.tsv
	}
	# sample specific filters
	if {$force || ![file exists reg_refcons-$name.tsv]} {
		puts stderr "Find refcons regions for var-$name.tsv"
		cg refconsregions svar-$name.tsv > reg_refcons-$name.tsv
	}
	if {$force || ![file exists reg_cluster-$name.tsv]} {
		puts stderr "Find cluster regions for svar-$name.tsv"
		cg clusterregions < annotvar-$name.tsv > reg_cluster-$name.tsv
	}
	if {$force || ![file exists reg_ns-$name.tsv]} {
		puts stderr "Find regions with N's for svar-$name.tsv"
		cg select -f {chromosome begin end} -q {$alleleSeq1 ~ /[N?]/ || $alleleSeq2 ~ /[N?]/} < annotvar-$name.tsv > reg_ns-$name.tsv
	}
	if {$force || ![file exists reg_lowscore-$name.tsv]} {
		puts stderr "Find regions with lowscores for svar-$name.tsv"
		cg select -f {chromosome begin end} -q {$totalScore1 < 60 || $totalScore2 < 60} < annotvar-$name.tsv > reg_lowscore-$name.tsv
	}
	if {$force || ![file exists fannotvar-$name.tsv]} {
		# add filterdata to annotvar
		file copy -force annotvar-$name.tsv temp.tsv
		annot_coverage $dir temp.tsv ftemp.tsv
		file copy -force ftemp.tsv temp.tsv
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
		list_foreach {field value regfile} $todo {
			if {![file exists $regfile]} {
				puts stderr "$field: $regfile does not exists"
				continue
			}
			puts "Annotating $field $value"
			cg annot_compare_region temp.tsv $regfile $field $value "" > tempfiltered.tsv
			if {[file exists tempfiltered.tsv]} {
				file rename -force tempfiltered.tsv temp.tsv
			}
		}
		file rename temp.tsv fannotvar-$name.tsv
	}
	# coverage
	if {$force || ![file exists reg-$name.covered]} {
		puts stderr "Coverage of sequenced regions"
		cg covered sreg-$name.tsv > reg-$name.covered
	}
	if {$force || ![file exists filteredrefcons-$name.covered]} {
		puts stderr "Coverage of refcons region"
		cg regsubtract sreg-$name.tsv reg_refcons-$name.tsv > filteredrefcons-$name.tsv
		cg covered filteredrefcons-$name.tsv > filteredrefcons-$name.covered
	}
	if {$force || ![file exists filteredns-$name.covered]} {
		puts stderr "Coverage of ns region"
		cg regsubtract sreg-$name.tsv reg_ns-$name.tsv > filteredns-$name.tsv
		cg covered filteredns-$name.tsv > filteredns-$name.covered
	}
	if {$force || ![file exists filteredlowscore-$name.covered]} {
		puts stderr "Coverage of lowscore region"
		cg regsubtract sreg-$name.tsv reg_lowscore-$name.tsv > filteredlowscore-$name.tsv
		cg covered filteredlowscore-$name.tsv > filteredlowscore-$name.covered
	}
	if {$force || ![file exists histo-refcons-$name.tsv]} {
		cg reghisto reg_refcons-$name.tsv > histo-refcons-$name.tsv
	}
	if {$force || ![file exists filteredcluster-$name.covered]} {
		puts stderr "Coverage of clusters region"
		cg regsubtract sreg-$name.tsv reg_cluster-$name.tsv > filteredcluster-$name.tsv
		cg covered filteredcluster-$name.tsv > filteredcluster-$name.covered
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
		cg compare_pvt < fcompar_${name1}_${name2}.tsv > pvtcompar_${name1}_${name2}.tsv
		cg compare_pvtsummary < pvtcompar_${name1}_${name2}.tsv > summarycompar_${name1}_${name2}.tsv
	}
	cd $keepdir
}

if 0 {

	set basedir /media/passport/complgen
	set basedir /complgen
	lappend auto_path ~/dev/completegenomics/lib
	package require Tclx
	signal -restart error SIGINT
	package require Extral
	set dir1 $basedir/GS102
	set dir2 $basedir/GS103
	set dbdir /data/db
	set resultsdir $basedir/testcompar_GS102_GS103
	set force 0
	set name1 [file tail $dir1]
	set name2 [file tail $dir2]
	set file1 $dir1/fannotvar-$name1.tsv
	set file2 $dir2/fannotvar-$name2.tsv
	set regfile1 $dir1/sreg-$name1.tsv
	set regfile2 $dir2/sreg-$name2.tsv
	set outfile compar_${name1}_${name2}.tsv
	file mkdir $resultsdir
	cd $resultsdir

	annot_compare_coverage compar_GS102_GS103.tsv $dir1 $dir2 ftemp.tsv

}
