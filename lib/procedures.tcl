
proc cg {args} {
	puts "cg $args"
	eval exec cg $args 2>@ stderr
}

proc process_sample {dir {force 0}} {
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
		cg select -s "chromosome begin" < $varfile > svar-$name.tsv
	}
	if {$force || ![file exists sgene-$name.tsv]} {
		puts stderr "Sort gene file ($genefile)"
		cg select -s "chromosome begin" < $genefile > sgene-$name.tsv
	}
	if {$force || ![file exists sreg-$name.tsv]} {
		puts stderr "Sort region file ($regfile)"
		cg select -s "chromosome begin" < $regfile > sreg-$name.tsv
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
		cg compare_annot $name1 $dir1/annotvar-$name1.tsv $dir1/sreg-$name1.tsv $name2 $dir2/annotvar-$name2.tsv $dir2/sreg-$name2.tsv compar_${name1}_${name2}.tsv
	}
	if {$force || ![file exists fcompar_${name1}_${name2}.tsv]} {
		file copy compar_${name1}_${name2}.tsv temp.tsv
	} else {
		file copy -force fcompar_${name1}_${name2}.tsv temp.tsv
	}
	set todo {}
	lappend todo [list refcons rc $dir1/reg_refcons-$name1.tsv]
	lappend todo [list refcons rc $dir2/reg_refcons-$name2.tsv]
	lappend todo [list cluster cl $dir1/reg_cluster-$name1.tsv]
	lappend todo [list cluster cl $dir2/reg_cluster-$name2.tsv]
	lappend todo [list trf trf $dbdir/regdb-simple_repeats.tsv]
	lappend todo [list str str $dbdir/regdb-microsatelite.tsv]
	lappend todo [list repeat rp $dbdir/regdb-repeatmasker.tsv]
	lappend todo [list segdup sd $dbdir/regdb-segdups.tsv]
	lappend todo [list rna rna $dbdir/regdb-rnagenes.tsv]
	lappend todo [list selfchain sc $dbdir/regdb-selfchain.tsv]
	lappend todo [list a100 a100 $dir1/ASM/REF/above100coverage.regions]
	lappend todo [list a100 a100 $dir2/ASM/REF/above100coverage.regions]
	lappend todo [list a70 a70 $dir1/ASM/REF/above70coverage.regions]
	lappend todo [list a70 a70 $dir2/ASM/REF/above70coverage.regions]
	lappend todo [list b30 b30 $dir1/ASM/REF/below30coverage.regions]
	lappend todo [list b30 b30 $dir2/ASM/REF/below30coverage.regions]
	lappend todo [list b20 b20 $dir1/ASM/REF/below20coverage.regions]
	lappend todo [list b20 b20 $dir2/ASM/REF/below20coverage.regions]
	lappend todo [list b15 b15 $dir1/ASM/REF/below15coverage.regions]
	lappend todo [list b15 b15 $dir2/ASM/REF/below15coverage.regions]
	lappend todo [list checked checked checked_${name1}_${name2}.tsv]
	list_foreach {field value regfile} $todo {
		if {![file exists $regfile]} {
			puts stderr "$field: $regfile does not exists"
			continue
		}
		if {!$force && [process_compare_checkfield temp.tsv $field]} {
			puts stderr "skipping $field: already present"
			continue
		}
		cg annot_compare_region temp.tsv $regfile $field $value "" > tempfiltered.tsv
		file delete temp.tsv
		file rename tempfiltered.tsv temp.tsv
	}
	file rename temp.tsv fcompar_${name1}_${name2}.tsv
	cg compare_pvt < fcompar_${name1}_${name2}.tsv > pvtcompar_${name1}_${name2}.tsv
	cg compare_pvtsummary < pvtcompar_${name1}_${name2}.tsv > summarycompar_${name1}_${name2}.tsv
	cd $keepdir
}

if 0 {
	lappend auto_path ~/dev/completegenomics/lib
	package require Tclx
	signal -restart error SIGINT
	package require Extral
	set dir /media/passport/complgen/GS102
	set dir1 /media/passport/complgen/GS102
	set dir2 /media/passport/complgen/GS103
	set dbdir /data/db
	set resultsdir /media/passport/complgen/compar_GS102_GS103
	set force 0
}
