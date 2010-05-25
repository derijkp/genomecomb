proc multicompar_annot_join {cur1 cur2} {
	global comparposs1 mergeposs1 comparposs2 mergeposs2 dummy1 dummy2 restposs1 restposs2
	global compare_annot_join_trans
	if {[inlist {{} -} $cur1]} {
		set region [list_sub $cur2 $comparposs2]
		set merge [list_sub $cur2 $mergeposs2]
		if {$cur1 eq "-"} {
			set cur1 $dummy1
		} else {
			set cur1 [list_change $dummy1 {- {}}]
		}
	} elseif {[inlist {{} -} $cur2]} {
		set region [list_sub $cur1 $comparposs1]
		set merge [list_sub $cur1 $mergeposs1]
		if {$cur2 eq "-"} {
			set cur2 $dummy2
		} else {
			set cur2 [list_change $dummy2 {- {}}]
		}
	} else {
		set region [list_sub $cur1 $comparposs1]
		set merge {}
		foreach el1 [list_sub $cur1 $mergeposs1] el2 [list_sub $cur2 $mergeposs2] {
			lappend merge [list_union $el1 $el2]
		}
	}
	set result $region
	lappend result {*}[list_sub $cur1 $restposs1]
	lappend result {*}[list_sub $cur2 $restposs2]
	lappend result {*}$merge
	return [join $result \t]
}

proc multicompar {file1 dir} {
	global cache comparposs1 mergeposs1 comparposs2 mergeposs2 dummy1 dummy2 restposs1 restposs2

	catch {close $f1}; catch {close $f2}; catch {close $o}
	set comparfields {chromosome begin end type}
	set mergefields {reference xRef trf str segdup selfchain repeat rna checked geneId mrnaAcc proteinAcc orientation exonCategory exon codingRegionKnown aaCategory nucleotidePos proteinPos aaAnnot aaCall aaRef effect neffect}
	set allelefields {alleleSeq1 alleleSeq2}
	#
	set name [file tail $dir]
	set file2 $dir/fannotvar-$name.tsv
	if {![file exists $file1]} {
		file_write $file1 [join {chromosome begin end type} \t]
	}
	set f1 [open $file1]
	set header1 [split [gets $f1] \t]
	if {[inlist $header1 alleleSeq1-$name]} {
		error "$name already present in $file1"
	}
	set comparposs1 [list_cor $header1 $comparfields]
	if {([lsearch $comparposs1 -1] != -1)} {
		puts stderr "header error in comparfile $file1"
		exit 1
	}
	set mergeposs1 [list_remove [list_cor $header1 $mergefields] -1]
	set dummy1 [list_fill [llength $header1] ?]
	set f2 [open $file2]
	set header2 [split [gets $f2] \t]
	set comparposs2 [list_cor $header2 $comparfields]
	if {([lsearch $comparposs2 -1] != -1)} {
		puts stderr "header error in fannot_varfile2 $file1"
		exit 1
	}
	set mergeposs2 [list_remove [list_cor $header2 $mergefields] -1]
	set dummy2 [list_fill [llength $header2] ?]
	# start
	set o [open $file1.temp w]
	# make output header
	set restfields1 [list_lremove [list_lremove $header1 $mergefields] $comparfields]
	set restposs1 [list_cor $header1 $restfields1]
	set oheader [list_concat $comparfields $restfields1]
	set restfields2 [list_lremove [list_lremove $header2 $mergefields] $comparfields]
	set restposs2 [list_cor $header2 $restfields2]
	foreach field $restfields2 {
		lappend oheader ${field}-$name
	}
	set oheader [list_concat $oheader $mergefields]
	puts $o [join $oheader \t]
	set cur1 [split [gets $f1] \t]
	set comp1 [list_sub $cur1 $comparposs1]
	set cur2 [split [gets $f2] \t]
	set comp2 [list_sub $cur2 $comparposs2]
	set num 1
	while {![eof $f1] || ![eof $f2]} {
		incr num
		if {![expr {$num % 100000}]} {putslog $num}
		set d [comparepos $comp1 $comp2]
		if {$d == 0} {
			puts $o [multicompar_annot_join $cur1 $cur2]
			set cur1 [compare_annot_getline $f1]
			set comp1 [list_sub $cur1 $comparposs1]
			set cur2 [compare_annot_getline $f2]
			set comp2 [list_sub $cur2 $comparposs2]
		} elseif {$d < 0} {
			while {[comparepos $comp1 $comp2] < 0} {
				puts $o [multicompar_annot_join $cur1 -]
				if {[eof $f1]} break
				set cur1 [compare_annot_getline $f1]
				set comp1 [list_sub $cur1 $comparposs1]
				if {![llength $cur1]} break
			}
		} else {
			while {[comparepos $comp1 $comp2] > 0} {
				puts $o [multicompar_annot_join - $cur2]
				if {[eof $f2]} break
				set cur2 [compare_annot_getline $f2]
				set comp2 [list_sub $cur2 $comparposs2]
				if {![llength $cur2]} break
			}
		}
	}
	close $f1; close $f2; close $o
	file rename -force $file1 $file1.old
	file rename $file1.temp $file1
}

proc multicompar_reannot {compar_file} {
	set compar_file [file normalize $compar_file]
	set basedir [file dir [file dir $compar_file]]
	catch {close $f}; catch {close $o}
	set f [open $compar_file]
	set header [split [gets $f] \t]
	set pos -1
	unset -nocomplain samplea
	set samples {}
	foreach field $header {
		incr pos
		set temp [split $field -]
		if {[llength $temp] == 2} {
			set sample [lindex $temp 1]
			lappend samplea(poss,$sample) $pos
			lappend samplea(fields,$sample) [lindex $temp 0]
			list_addnew samples $sample
		}
	}
	set referencepos [lsearch $header reference]
	foreach sample $samples {
		set samplea(a1,$sample) [lsearch $header alleleSeq1-$sample]
		set samplea(a2,$sample) [lsearch $header alleleSeq2-$sample]
		set samplea(rpos,$sample) [lsearch $header refscore-$sample]
		set samplea(cpos,$sample) [lsearch $header coverage-$sample]
		set samplea(regionfile,$sample) $basedir/$sample/sreg-$sample.tsv
		if {[inlist $samplea(fields,$sample) refcons]} {
			lappend samplea(todo,$sample) [list [lsearch $header refcons-$sample] rc $basedir/$sample/reg_refcons-$sample.tsv]
		}
		if {[inlist $samplea(fields,$sample) cluster]} {
			lappend samplea(todo,$sample) [list [lsearch $header cluster-$sample] cl $basedir/$sample/reg_cluster-$sample.tsv]
		}
		list_foreach {field value regfile} $samplea(todo,$sample) {
			annot_region_init $regfile
		}
		annot_coverage_init $basedir/$sample
		annot_region_init $samplea(regionfile,$sample)
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
		foreach {chr begin end} [list_sub $line $poss] break
		foreach sample $samples {
			if {[lindex $line $samplea(a1,$sample)] ne "?"} continue
			if {([lindex $line $samplea(rpos,$sample)] eq "?") || ([lindex $line $samplea(cpos,$sample)] eq "?")} {
				foreach {r c} [annot_coverage_get $basedir/$sample $chr $begin] break
				lset line $samplea(rpos,$sample) $r
				lset line $samplea(cpos,$sample) $c
			}
			list_foreach {field value regfile} $samplea(todo,$sample) {
				if {[lindex $line $field] == "-"} continue
				if {[lindex $line $field] != "?"} continue
				set r [annot_region_get $regfile $chr $begin $end]
				if {$r} {lset line $field $value} else {lset line $field {}}
			}
			set r [annot_region_get $samplea(regionfile,$sample) $chr $begin $end]
			if {$r} {
				lset line $samplea(a1,$sample) $reference
				lset line $samplea(a2,$sample) $reference
			} else {
				lset line $samplea(a1,$sample) -
				lset line $samplea(a2,$sample) -
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
		annot_coverage_close $basedir/$sample
		annot_region_close $samplea(regionfile,$sample)
	}
	file rename -force $compar_file $compar_file.old
	file rename $compar_file.temp $compar_file
}

if 0 {

	lappend auto_path ~/dev/completegenomics/lib
	lappend auto_path /complgen/bin/complgen/apps/cg/lib
	package require Extral
	package require Tclx
	signal -restart error SIGINT
set compar_file /complgen/multicompar/compar-X.tsv

	set basedir /media/passport/complgen
	set basedir /complgen
	set dbdir /complgen/refseq
	set file1 /complgen/multicompar/compar.tsv
	set dir /complgen/GS102
	set dir /complgen/GS103
	set dir /complgen/GS100
multicompar $file1 $dir
	set dir /complgen/GS101A01
multicompar $file1 $dir
	set dir /complgen/GS101A02
multicompar $file1 $dir

	set compar_file /complgen/multicompar/compar.tsv

}
