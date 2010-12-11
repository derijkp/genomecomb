
proc cg {args} {
	# puts "cg $args"
	eval exec cg $args 2>@ stderr
}

proc process_sample {dir destdir dbdir {force 0}} {
	set keepdir [pwd]
	set dir [file normalize $dir]
	set destdir [file normalize $destdir]
	set dbdir [file normalize $dbdir]
	file mkdir $destdir
	cd $destdir
	set name [file tail $destdir]
	set chromosomes {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 M X Y}
	puts stderr "Processing sample $dir -> $destdir"
	# sort files
	if {$force || ![file exists svar-$name.tsv] || ![file exists $destdir/info.txt]} {
		set varfile [glob $dir/ASM/var-*-ASM*.tsv*]
		if {[llength $varfile] != 1} {error "could not identify varfile"}
		if {[llength $varfile] > 1} {error "could not identify varfile"}
		set f [rzopen $varfile]
		set info {}
		while {![eof $f]} {
			set line [gets $f]
			if {[string index $line 0] ne "#"} break
			lappend info $line
		}
		catch {close $f}
		if {[file exists $destdir/info.txt]} {
			set test [split [file_read $destdir/info.txt] \n]
			if {$info ne $test} {
				error "$destdir already has info.txt that differs from data in the source $dir"
			}
		}
		file_write $destdir/info.txt [join $info \n]
		puts stderr "Sort var file ($varfile)"
		cg select -s "chromosome begin end varType" $varfile temp.tsv
		file rename -force temp.tsv svar-$name.tsv
	}
	if {$force || ![file exists sgene-$name.tsv]} {
		set genefile [list_lremove [glob $dir/ASM/gene-*-ASM*.tsv*] [glob -nocomplain $dir/ASM/gene-var-summary-*-ASM*.tsv*]]
		if {[llength $genefile] != 1} {error "could not identify genefile"}
		puts stderr "Sort gene file ($genefile)"
		cg select -s "chromosome begin end" $genefile temp.tsv
		file rename -force temp.tsv sgene-$name.tsv
	}
	if {$force || ![file exists sreg-$name.tsv]} {
		puts stderr "Sort region file ($regfile)"
		set regfile [glob -nocomplain $dir/ASM/reg-*-ASM*.tsv*]
		if {[file exists $regfile]} {
			cg select -s "chromosome begin end" $regfile temp.tsv
			file rename -force temp.tsv sreg-$name.tsv
		} else {
			cg select -q {$varType != "no-call" && $varType != "no-ref"} -f "chromosome begin end" svar-$name.tsv temp.tsv
			cg regjoin temp.tsv > temp2.tsv
			file rename -force temp2.tsv sreg-$name.tsv
		}
	}
	# sort coverage files
	file mkdir coverage
	foreach {chr} $chromosomes  {
		set resultfile coverage/coverageRefScore-$chr-$name.tsv
		if {$force || (![file exists $resultfile] && ![file exists $resultfile.rz])} {
			set oricov [lindex [glob -nocomplain \
				$dir/ASM/REF/coverageRefScore-$chr-*-ASM*.tsv \
				$dir/ASM/REF/coverageRefScore-chr$chr-*-ASM*.tsv \
				$dir/ASM/REF/coverageRefScore-$chr-*-ASM*.tsv.* \
				$dir/ASM/REF/coverageRefScore-chr$chr-*-ASM*.tsv.*] 0]
			if {[file exists $oricov]} {
				putslog "Sorting to create $resultfile"
				cg select -s offset $oricov coverage/temp.tsv
				file rename -force coverage/temp.tsv $resultfile
			}
		}
		if {$force || ![file exists $resultfile.offset_index]} {
			if {[file exists $resultfile]} {
				putslog "Creating index $resultfile.offset_index"
				process_indexcompress $resultfile
			}
		}
	}
	# annotated vars file
	if {$force || ![file exists annotvar-$name.tsv]} {
		puts stderr "Create annotated varfile annotvar-$name.tsv"
		# set file svar-$name.tsv; set genefile sgene-$name.tsv; set outfile temp.tsv
		cg var2annot svar-$name.tsv sgene-$name.tsv temp.tsv
		file rename -force temp.tsv annotvar-$name.tsv
	}
	# sample specific filters
	if {$force || ![file exists reg_refcons-$name.tsv]} {
		puts stderr "Find refcons regions for var-$name.tsv"
		cg refconsregions svar-$name.tsv > temp.tsv
		file rename -force temp.tsv reg_refcons-$name.tsv
	}
	if {$force || ![file exists reg_nocall-$name.tsv]} {
		puts stderr "Find partial no-call regions for var-$name.tsv"
		if {[catch {
			nocallregions svar-$name.tsv temp.tsv
		}]} {
			puts stderr "Could not make reg_nocall-$name.tsv (old version files ?)"
		} else {
			file rename -force temp.tsv reg_nocall-$name.tsv
		}
	}
	if {$force || ![file exists reg_cluster-$name.tsv]} {
		puts stderr "Find cluster regions for svar-$name.tsv"
		cg clusterregions < annotvar-$name.tsv > temp.tsv
		file rename -force temp.tsv reg_cluster-$name.tsv
	}
	if {$force || ![file exists reg_ns-$name.tsv]} {
		puts stderr "Find regions with N's for svar-$name.tsv"
		cg select -f {chromosome begin end} -q {$alleleSeq1 ~ /[N?]/ || $alleleSeq2 ~ /[N?]/} < annotvar-$name.tsv > temp.tsv
		file rename -force temp.tsv reg_ns-$name.tsv
	}
	if {$force || ![file exists reg_lowscore-$name.tsv]} {
		puts stderr "Find regions with lowscores for svar-$name.tsv"
		cg select -f {chromosome begin end} -q {$totalScore1 < 60 || $totalScore2 < 60} < annotvar-$name.tsv > temp.tsv
		file rename -force temp.tsv reg_lowscore-$name.tsv
	}
	if {$force || ![file exists fannotvar-$name.tsv]} {
		# add filterdata to annotvar
		set todo {}
		lappend todo [list refcons rc reg_refcons-$name.tsv]
		lappend todo [list nocall nc reg_nocall-$name.tsv]
		lappend todo [list cluster cl reg_cluster-$name.tsv]
		lappend todo [list trf trf $dbdir/regdb-simple_repeats.tsv]
		lappend todo [list str str $dbdir/regdb-microsatelite.tsv]
		lappend todo [list segdup sd $dbdir/regdb-segdups.tsv]
		lappend todo [list selfchain sc $dbdir/regdb-selfchain.tsv]
		lappend todo [list repeat rp $dbdir/regdb-repeatmasker.tsv]
		lappend todo [list rna rna $dbdir/regdb-rnagenes.tsv]
		lappend todo [list more5pct m5 $dbdir/regdb-1000genomesmore5pct.tsv]
		lappend todo [list more1pct m1 $dbdir/regdb-1000genomesmore1pct.tsv]
		foreach file [lsort -dictionary [glob -nocomplain $dbdir/checked*.tsv]] {
			set value [lindex [split [file root [file tail $file]] _] 1]
			if {$value eq ""} {set value checked}
			lappend todo [list checked $value $file]
		}
		annot_annotvar annotvar-$name.tsv fannotvar-$name.tsv $todo $destdir
	}
	# coverage
	if {$force || ![file exists reg-$name.covered]} {
		puts stderr "Coverage of sequenced regions"
		cg covered sreg-$name.tsv > temp.tsv
		file rename -force temp.tsv reg-$name.covered
	}
	if {$force || ![file exists filteredrefcons-$name.covered]} {
		puts stderr "Coverage of refcons region"
		cg regsubtract sreg-$name.tsv reg_refcons-$name.tsv > temp.tsv
		file rename -force temp.tsv filteredrefcons-$name.tsv
		cg covered filteredrefcons-$name.tsv > temp.tsv
		file rename -force temp.tsv filteredrefcons-$name.covered
	}
	if {$force || ![file exists filteredns-$name.tsv]} {
		puts stderr "Coverage of ns region"
		cg regsubtract sreg-$name.tsv reg_ns-$name.tsv > temp.tsv
		file rename -force temp.tsv filteredns-$name.tsv
	}
	if {$force || ![file exists filteredns-$name.covered]} {
		cg covered filteredns-$name.tsv > temp.tsv
		file rename -force temp.tsv filteredns-$name.covered
	}
	if {$force || ![file exists filteredlowscore-$name.tsv]} {
		puts stderr "Coverage of lowscore region"
		cg regsubtract sreg-$name.tsv reg_lowscore-$name.tsv > temp.tsv
		file rename -force temp.tsv filteredlowscore-$name.tsv
	}
	if {$force || ![file exists filteredlowscore-$name.covered]} {
		cg covered filteredlowscore-$name.tsv > temp.tsv
		file rename -force temp.tsv filteredlowscore-$name.covered
	}
	if {$force || ![file exists histo-refcons-$name.tsv]} {
		cg reghisto reg_refcons-$name.tsv > temp.tsv
		file rename -force temp.tsv histo-refcons-$name.tsv
	}
	if {$force || ![file exists filteredcluster-$name.tsv]} {
		puts stderr "Coverage of clusters region"
		cg regsubtract sreg-$name.tsv reg_cluster-$name.tsv > temp.tsv
		file rename -force temp.tsv filteredcluster-$name.tsv
	}
	if {$force || ![file exists filteredcluster-$name.covered]} {
		cg covered filteredcluster-$name.tsv > temp.tsv
		file rename -force temp.tsv filteredcluster-$name.covered
	}
	if {$force || ![file exists reg_below20-$name.tsv]} {
		puts stderr "Make region file reg_below20-$name.tsv"
		file delete temp.tsv
		set f [open temp.tsv w]
		puts $f "chromosome\tbegin\tend"
		close $f
		foreach {chr} $chromosomes  {
			set file [lindex [glob -nocomplain coverage/coverageRefScore-$chr-$name.tsv coverage/coverageRefScore-$chr-$name.tsv.rz coverage/coverageRefScore-$chr-$name.tsv.gz] 0]
			if {($file eq "" ) && ($chr eq "Y")} continue
			puts stderr "Processing $file"
			set f [rzopen $file]
			while {![eof $f]} {
				set header [gets $f]
				if {[string index $header 0] ne "#"} break
			}
			catch {close $f}
			set poscol [lsearch $header offset]
			set coveragecol [lsearch $header coverage]
			if {[inlist {.rz .gz} [file extension $file]]} {set cat zcat} else {set cat cat}
			set error [catch {
				exec $cat $file | getregions $chr $poscol $coveragecol 20 0 0 >> temp.tsv
			} errmessage]
			if {$error && ![regexp {decompression OK, trailing garbage ignored} $errmessage]} {
				error $errmessage
			}
		}
		file rename -force temp.tsv reg_below20-$name.tsv
	}
	if {$force || ![file exists reg_below20-$name.covered]} {
		puts stderr "Make reg_below20-$name.covered"
		cg covered reg_below20-$name.tsv > temp.tsv
		file rename temp.tsv reg_below20-$name.covered
	}
	if {$force || ![file exists reg_below10-$name.tsv]} {
		puts stderr "Make region file reg_below10-$name.tsv"
		file delete temp.tsv
		set f [open temp.tsv w]
		puts $f "chromosome\tbegin\tend"
		close $f
		foreach {chr} $chromosomes  {
			set file [lindex [glob -nocomplain coverage/coverageRefScore-$chr-$name.tsv coverage/coverageRefScore-$chr-$name.tsv.rz coverage/coverageRefScore-$chr-$name.tsv.gz] 0]
			if {($file eq "" ) && ($chr eq "Y")} continue
			puts stderr "Processing $file"
			set f [rzopen $file]
			while {![eof $f]} {
				set header [gets $f]
				if {[string index $header 0] ne "#"} break
			}
			catch {close $f}
			set poscol [lsearch $header offset]
			set coveragecol [lsearch $header coverage]
			if {[inlist {.rz .gz} [file extension $file]]} {set cat zcat} else {set cat cat}
			set error [catch {
				exec $cat $file | getregions $chr $poscol $coveragecol 10 0 0 >> temp.tsv
			} errmessage]
			if {$error && ![regexp {decompression OK, trailing garbage ignored} $errmessage]} {
				error $errmessage
			}
		}
		file rename -force temp.tsv reg_below10-$name.tsv
	}
	if {$force || ![file exists reg_below10-$name.covered]} {
		puts stderr "Make reg_below10-$name.covered"
		cg covered reg_below10-$name.tsv > temp.tsv
		file rename temp.tsv reg_below10-$name.covered
	}
	if {$force || ![file exists reg_above100-$name.tsv]} {
		puts stderr "Make region file reg_above100-$name.tsv"
		file delete temp.tsv
		set f [open temp.tsv w]
		puts $f "chromosome\tbegin\tend"
		close $f
		foreach {chr} $chromosomes  {
			set file [lindex [glob -nocomplain coverage/coverageRefScore-$chr-$name.tsv coverage/coverageRefScore-$chr-$name.tsv.rz coverage/coverageRefScore-$chr-$name.tsv.gz] 0]
			if {($file eq "" ) && ($chr eq "Y")} continue
			puts stderr "Processing $file"
			set f [rzopen $file]
			while {![eof $f]} {
				set header [gets $f]
				if {[string index $header 0] ne "#"} break
			}
			catch {close $f}
			set poscol [lsearch $header offset]
			set coveragecol [lsearch $header coverage]
			if {[inlist {.rz .gz} [file extension $file]]} {set cat zcat} else {set cat cat}
			set error [catch {
				exec $cat $file | getregions $chr $poscol $coveragecol 100 1 0 >> temp.tsv
			} errmessage]
			if {$error && ![regexp {decompression OK, trailing garbage ignored} $errmessage]} {
				error $errmessage
			}
		}
		file rename -force temp.tsv reg_above100-$name.tsv
	}
	if {$force || ![file exists reg_above100-$name.covered]} {
		puts stderr "Make reg_above100-$name.covered"
		cg covered reg_above100-$name.tsv > temp.tsv
		file rename temp.tsv reg_above100-$name.covered
	}
	puts stderr "Finished $destdir"
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

proc process_indexcompress {file} {
	set ext [file extension $file]
	if {$ext eq ".gz"} {
		gunzip $file
		set file [file root $file]
	}
	set f [open $file]
	set header [tsv_open $f]
	foreach field {offset end1} {
		set fpos [lsearch $header $field]
		if {$fpos != -1} break
	}
	if {$fpos == -1} {error "no column offset, end1 in file $file"}
	if {([gets $f] eq "") && [eof $f]} return
	close $f
	if {![file exists $file.${field}_index]} {
		tsv_index $field $file
	}
	putslog "Compressing $file"
	exec razip -c $file > $file.rz.temp
	file rename -force $file.rz.temp $file.rz
	file delete $file
}

proc process_compare_coverage {dir1 dir2 dbdir resultsdir {force 0}} {
	set dir1 [file normalize $dir1]
	set dir2 [file normalize $dir2]
	set dbdir [file normalize $dbdir]
	file mkdir $resultsdir
	cd $resultsdir
	# region files
	set name1 [file tail $dir1]
	set name2 [file tail $dir2]
	set target ${name1}_${name2}_regcommon.tsv.covered
	if {$force || ![file exists $target]} {
		putslog "Making $target"
		cg regsubtract $dir1/sreg-${name1}.tsv $dir2/sreg-${name2}.tsv > temp.tsv
		cg regsubtract $dir1/sreg-${name1}.tsv temp.tsv > ${name1}_${name2}_regcommon.tsv
		cg covered ${name1}_${name2}_regcommon.tsv > $target
	}
	# region from each sample
	set filters {}
	lappend filters refcons [list $dir1/reg_refcons-${name1}.tsv $dir2/reg_refcons-${name2}.tsv]
	lappend filters ns [list $dir1/reg_ns-${name1}.tsv $dir2/reg_ns-${name2}.tsv]
	lappend filters str [list $dbdir/regdb-microsatelite.tsv {}]
	lappend filters trf [list $dbdir/regdb-simple_repeats.tsv {}]
	lappend filters cluster [list $dir1/reg_cluster-${name1}.tsv $dir2/reg_cluster-${name2}.tsv]
	lappend filters rtg [list $dir1/reg_rtg-${name1}.tsv $dir2/reg_rtg-${name2}.tsv]
	lappend filters segdup [list $dbdir/regdb-segdups.tsv {}]
	lappend filters b10 [list $dir1/reg_below10-${name1}.tsv $dir2/reg_below10-${name2}.tsv]
	lappend filters b20 [list $dir1/reg_below20-${name1}.tsv $dir2/reg_below20-${name2}.tsv]
	lappend filters a100 [list $dir1/reg_above100-${name1}.tsv $dir2/reg_above100-${name2}.tsv]
	lappend filters lowscore [list $dir1/reg_lowscore-${name1}.tsv $dir2/reg_lowscore-${name2}.tsv]
	lappend filters selfchain [list $dbdir/regdb-selfchain.tsv {}]
	lappend filters repeat [list $dbdir/regdb-repeatmasker.tsv {}]
	set todo [list_unmerge $filters]
	foreach dbname $todo {
		foreach {rfile1 rfile2} [dict get $filters $dbname] break
		set target ${name1}_${name2}_regcommon-$dbname.tsv.covered
		if {$force || ![file exists $target]} {
			putslog "Making $target"
			cg regsubtract ${name1}_${name2}_regcommon.tsv $rfile1 > temp.tsv
			if {$rfile2 ne ""} {
				file rename -force temp.tsv temp-old.tsv
				cg regsubtract temp-old.tsv $rfile2 > temp.tsv
			}
			file rename -force temp.tsv ${name1}_${name2}_regcommon-$dbname.tsv
			cg covered ${name1}_${name2}_regcommon-$dbname.tsv > temp.covered
			file rename -force temp.covered $target
		} else {
			putslog "$target already there"
		}
	}
	# subsequent subtractions (filters)
	set orders {
		{ns trf str cluster rtg segdup b20 a100 lowscore selfchain repeat}
		{ns trf str cluster rtg segdup b10 a100 lowscore selfchain repeat}
		{ns trf str segdup cluster lowscore b20 a100 rtg}
		{ns trf str segdup cluster lowscore b10 a100 rtg}
		{ns trf str segdup cluster lowscore rtg}
		{ns b20 a100 cluster lowscore trf str segdup rtg}
		{ns b10 a100 cluster lowscore trf str segdup rtg}
		{ns cluster lowscore trf str segdup rtg}
		{ns b20 a100 lowscore cluster trf str segdup rtg}
		{ns b10 a100 lowscore cluster trf str segdup rtg}
		{ns lowscore cluster trf str segdup rtg}
	}
	foreach order $orders {
		set newname ${name1}_${name2}_regcommon-refcons
		foreach dbname $order {
			foreach {rfile1 rfile2} [dict get $filters $dbname] break
			set oldname $newname
			set newname $oldname-$dbname
			set target $newname.tsv.covered
			if {$force || ![file exists $target]} {
				putslog "Making $target"
				cg regsubtract $oldname.tsv  $rfile1 > temp.tsv
				if {$rfile2 ne ""} {
					file rename -force temp.tsv temp-old.tsv
					cg regsubtract temp-old.tsv $rfile2 > temp.tsv
				}
				file rename -force temp.tsv $newname.tsv
				cg covered $newname.tsv > temp.covered
				file rename -force temp.covered $target
			} else {
				putslog "$target already there"
			}
		}
	}
	exec grep total *.covered > summary_genomecoverage_${name1}_${name2}.tsv
}

if 0 {

	lappend auto_path /complgen/bin/complgen/apps/cg/lib
	lappend auto_path ~/dev/completegenomics/lib
	package require Tclx
	signal -restart error SIGINT
	package require Extral

78638

	set dbdir /complgen/refseq
	set basedir /complgen/projects/dlb1
	set dir $basedir/NNN
	set destdir $basedir/dlb_d_d388
	set force 0

	set dir1 $basedir/GS102
	set dir2 $basedir/GS103
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

