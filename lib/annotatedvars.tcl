#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc assert {check message} {
	if {![uplevel expr $check]} {
		error $message
	}
}

proc readlocus {f {pos 0}} {
	global cache
	set cor1 [get cache(cor1) ""]
	set line $cache($f)
	set locus [lindex $line $pos]
	set list [list $line]
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {[llength $cor1]} {set line [list_sub $line $cor1]}
		if {![llength $line]} continue
		if {[lindex $line $pos] eq $locus} {
			lappend list $line
		} else {
			while {![eof $f] && ![llength $line]} {
				set line [split [gets $f] \t]
				if {[llength $cor1]} {set line [list_sub $line $cor1]}
			}
			set cache($f) $line
			break
		}
	}
	return $list
}

proc joinvarlines list {
	set line [lindex $list 0] 
	if {[llength $list] <= 1} {
		return $line
	}
	set ref [join [list_subindex $list 6] {}]
	lset line 6 $ref
	set allele [join [list_subindex $list 7] {}]
	lset line 7 $allele
	lset line 4 [lindex $list end 4]
}

proc var2annotvar_readonevar_merge {list} {
	global cache
	set list [lsort -integer -index 3 [lsort -integer -index 4 $list]]

	set keeplist $list
	# make "ali"
	set firstpos [lindex $list 0 3]
	set extrapos $cache(extrapos)
	set extranum $cache(extranum)
	set seq(1) {}; set seq(2) {}; set seq(1,i) {}; set seq(2,i) {}; set seq(1,q) {}; set seq(2,q) {}; set seq(1,e) {}; set seq(2,e) {}
	set i 0
	list_foreach {bin hp chr start end type ref oalt score} $list {
		if {$start == $end} {
			lset list $i 5 ins
		}
		set len [expr {$end-$start}]
		set alt $oalt
		if {[regexp {\?} $alt]} {
			set alt ?
		} elseif {[string length $oalt] != [string length $ref]} {
			regsub -all . $alt \@ alt
		}
		if {$len > [string length $alt]} {
			if {$alt eq "?"} {set r ?} else {set r -}
			append alt [string_fill $r [expr {$len-[string length $alt]}]]
		}
		if {$len < [string length $alt]} {
			lappend seq($hp,i) [expr {$end-$firstpos}] [string range $alt $len end]
			set alt [string range $alt 0 [expr {$len-1}]]
		}
		if {[string length $alt]} {
			set s [expr {$start-$firstpos}]
			set e [expr {$end-1-$firstpos}]
			set seq($hp) [string_replace $seq($hp) $s $e $alt]
			if {[llength $seq($hp,q)] <= $e} {
				lappend seq($hp,q) {*}[list_fill [expr {$e-[llength $seq($hp,q)]+1}] {}]
			}
			set seq($hp,q) [lreplace $seq($hp,q) $s $e {*}[list_fill [expr {$e-$s+1}] $score]]
			if {$extranum} {
				set extra [lrange [lindex $list $i] $extrapos end]
				if {[llength $seq($hp,e)] <= $e} {
					lappend seq($hp,e) {*}[list_fill [expr {$e-[llength $seq($hp,e)]+1}] {}]
				}
				set seq($hp,e) [lreplace $seq($hp,e) $s $e {*}[list_fill [expr {$e-$s+1}] $extra]]
			}
		}
		incr i
	}

	set types [list_subindex $list 5]
	set poss [list_find -regexp $types {snp|ins|del|delins|sub}]
	set wlist [list_sub $list $poss]
	lappend wlist {}
	set line1 [list_shift wlist]
	foreach {bin1 hp1 chr1 begin1 end1 type1 ref1} $line1 break
	set line2 {}
	set rlist {}
	while {[llength $wlist]} {
		set line [list_shift wlist]
		foreach {bin hp chr begin end type ref} $line break
		if {$line ne "" && $chr eq $chr1 && $begin == $begin1 && $end == $end1} {
			lappend rlist $line1 $line
			set line1 [list_shift wlist]
			foreach {bin1 hp1 chr1 begin1 end1 type1 ref1} $line1 break
		} else {
			if {$hp1 == 1} {set uhp 2} else {set uhp 1}
			set line2 $line1
			set s [expr {$begin1-$firstpos}]
			set e [expr {$end1-$firstpos-1}]
			set alt [string range $seq($uhp) $s $e]
			foreach {iseq pos} [list_reverse $seq($uhp,i)] {
				if {$pos >= $s && $pos < $e} {
					set alt [string_replace $alt [expr {$pos-$s}] -1 $iseq]
				}
			}
			regsub -all -- - $alt {} alt
			if {[regexp {[?]} $alt]} {set alt ?}
			if {[regexp {@} $alt]} {set alt @}
			lset line2 8 [lindex $seq($uhp,q) $s]
			lset line2 9 {}
			lset line2 7 $alt
			if {$extranum} {
				set extra [lindex $seq($uhp,e) $s]
				set line2 [lreplace $line2 $extrapos end {*}$extra]
			}
			lappend rlist $line1 $line2
			set line1 $line
			foreach {bin1 hp1 chr1 begin1 end1 type1 ref1} $line1 break
		}
		if {![llength $line1]} break
	}
	return $rlist
}

proc cg_var2annot {args} {
	global scriptname action
	if {[llength $args] != 3} {
		puts stderr "format is: $scriptname $action variation_file sorted_gene_file output_file"
		puts stderr " - makes a new variation file, where the gene information is included"
		exit 1
	}
	foreach {file genefile outfile} $args break
	var2annotvar $file $genefile $outfile
}

proc var2annotvar_readonevar f {
	global cache list
	# join $cache($f,rov) \n
	set extrapos $cache(extrapos)
	set extranum $cache(extranum)
	if {[info exists cache($f,rov)]} {
		set line1 [list_shift cache($f,rov)]
		set line2 [list_shift cache($f,rov)]
		if {![llength $cache($f,rov)]} {unset cache($f,rov)}
	} else {
		set line1 {}
		set line2 {}
		while {![eof $f]} {
			while {![eof $f]} {
				set list [readlocus $f]
				set types [list_subindex $list 5]
				set type [list_remdup [list_remove $types = ref-consistent no-call-rc ref no-ref no-call no-ref PAR-called-in-X]]
				if {$type ne ""} break
			}
			if {[inlist $types sub]} {
				set temp {}
				foreach line $list {
					if {[lindex $line 5] eq "sub" && [string length [lindex $line 6]] <= 3} {
						set extra [lrange $line $extrapos end]
						foreach {bin al chr start end type ref alt s1 s2 s3} $line break
						set len [string length $ref]
						if {$len == [string length $alt]} {
							string_foreach r $ref e $alt {
								if {$r ne $e} {
									lappend temp [list $bin $al $chr $start [expr {$start+1}] snp $r $e $s1 $s2 $s3 {*}$extra]
								} else {
									lappend temp [list $bin $al $chr $start [expr {$start+1}] ref $r $r {} {} {} {*}$extra]
								}
								incr start
							}
						} else {
							lappend temp $line
						}
					} else {
						lappend temp $line
					}
				}
				set list $temp
			}
			if {$type eq ""} {return {}}
#putsvars list
			set rlist [var2annotvar_readonevar_merge $list]
			set line1 [list_shift rlist]
			set line2 [list_shift rlist]
			if {[llength $rlist]} {
				set cache($f,rov) $rlist
			}
			if {![llength $line1]} continue
			break
		}
	}
	if {![llength $line1]} {
		return {}
	}
	set extrapos $cache(extrapos)
	foreach {locus haplotype chromosome begin end varType reference alleleSeq totalScore hapLink xRef} $line1 break
	foreach {locus2 haplotype2 chromosome2 begin2 end2 varType2 reference2 alleleSeq2 totalScore2 hapLink2 xRef2} $line2 break
	set type [list [get varType unkown] [get varType2 unknown]]
	set type [list_remdup [list_remove $type = ref-consistent no-call-rc ref-inconsistent ref no-ref no-call no-ref PAR-called-in-X unknown]]
	if {$type ne ""} {
		set alleleSeq2 [get alleleSeq2 $reference]
		set totalScore2 [get totalScore2 ""]
#		if {$alleleSeq2 < $alleleSeq} {
#			set temp $alleleSeq ; set alleleSeq $alleleSeq2 ; set alleleSeq2 $temp
#			set temp $totalScore ; set totalScore $totalScore2 ; set totalScore2 $temp
#		}
		if {$alleleSeq eq "?"} {set alleleSeq -}
		if {$alleleSeq2 eq "?"} {set alleleSeq2 -}
		set result [list $locus $chromosome $begin $end [join $type _] $reference $alleleSeq $alleleSeq2 $totalScore $totalScore2 $xRef]
		if {$extranum} {
			foreach v1 [lrange $line1 $extrapos end] v2 [lrange $line2 $extrapos end] {
				lappend result $v1 $v2
			}
		}
		return $result
	} else {
		return {}
	}
}

if 0 {
	set source temp.tsv
	set dest tempfiltered.tsv

	set dir /complgen/GS102
	set source testannotvar.tsv
	set dest ftestannotvar.tsv
}

proc annot_coverage {dir source dest} {

	catch {close $o} ; catch {close $f} ; catch {close $fc}
	set f [open $source]
	set o [open $dest w]
	set header [split [gets $f] \t]
	set len [llength $header]
	set ipos [lsearch $header totalScore2]
	incr ipos
	puts $o [join [linsert $header $ipos refscore coverage] \t]
	set poss [list_cor $header {chromosome begin end}]
	annot_coverage_init $dir
	set num 0
	while {![eof $f]} {
		incr num
		if {![expr {$num%10000}]} {putsprogress $num}
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		if {[llength $line] < $len} {
			lappend line {*}[list_fill [expr {$len-[llength $line]}] {}]
		} elseif {[llength $line] > $len} {
			set line [lrange $line 0 [expr {$len-1}]]
		}
		foreach {chr begin end} [list_sub $line $poss] break
		foreach {refscore coverage} {{} {}} break
		foreach {refscore coverage} [annot_coverage_get $dir $chr $begin] break
		set line [linsert $line $ipos $refscore $coverage]
		puts $o [join $line \t]
	}
	annot_coverage_close $dir

	close $o
	close $f
}

if 0 {
	set name GS102
	set dbdir /data/db
	set source /complgen/$name/annotvar-$name.tsv
	set outfile /complgen/$name/fannotvar-$name.tsv
	set dir1 /complgen/$name
	set todo {}
	lappend todo [list refcons rc $dir/reg_refcons-$name.tsv]
	lappend todo [list cluster cl $dir/reg_cluster-$name.tsv]
	lappend todo [list trf trf $dbdir/regdb-simple_repeats.tsv]
	lappend todo [list str str $dbdir/regdb-microsatelite.tsv]
	lappend todo [list segdup sd $dbdir/regdb-segdups.tsv]
	lappend todo [list selfchain sc $dbdir/regdb-selfchain.tsv]
	lappend todo [list repeat rp $dbdir/regdb-repeatmasker.tsv]
	lappend todo [list rna rna $dbdir/regdb-rnagenes.tsv]
}

# todo list of lists: field value regfile
proc annot_annotvar {source outfile todo {dir {}}} {
	puts "annotating $source -> $outfile"

	catch {close $f} ; catch {close $o}
	set f [open $source]
	set header [tsv_open $f]
	if {$dir ne ""} {
		annot_coverage_init $dir
		set ipos [lsearch $header totalScore2]
		incr ipos
		set header [linsert $header $ipos refscore coverage]
		set addref 1
	} else {
		set addref 0
	}
	set extra {}
	set lpos 0
	set wtodo $todo
	list_foreach {field value regfile} $wtodo {
		puts "Init $regfile"
		set pos [lsearch $header $field]
		if {$pos == -1} {
			lappend header $field
			set pos [lsearch $header $field]
			lappend extra {}
		}
		lset wtodo $lpos 0 $pos 
		annot_region_init $regfile
		incr lpos
	}
	set o [open $outfile.temp w]
	set poss [list_cor $header {chromosome begin end}]
	puts $o [join $header \t]
	set num 0
	while {![eof $f]} {
		incr num
		if {![expr $num%10000]} {putsprogress $num}
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		foreach {chr begin end} [list_sub $line $poss] break
		if {$addref} {
			# refscore
			foreach {refscore coverage} {{} {}} break
			foreach {refscore coverage} [annot_coverage_get $dir $chr $begin] break
			set line [linsert $line $ipos $refscore $coverage]
		}
		# regions
		lappend line {*}$extra
		list_foreach {field value regfile} $wtodo {
			if {[lindex $line $field] == "-"} continue
			# if {[lindex $line $field] != ""} continue
			set r [annot_region_get $regfile $chr $begin $end]
			if {$r} {lset line $field $value}
		}
		puts $o [join $line \t]
	}
	if {$addref} {annot_coverage_close $dir}
	close $o
	close $f
	list_foreach {field value regfile} $wtodo {
		annot_region_close $regfile
	}
	file rename -force $outfile.temp $outfile

}


proc chrcomp {chr1 chr2} {
	if {$chr1 eq $chr2} {return 0}
	return [expr {[chr2num $chr1] - [chr2num $chr2]}]
}

proc readgeneset {g} {
	global cache
	set line $cache($g)
	set chr [lindex $line 3]
	set begin [lindex $line 4]
	set list [list $line]
	while {![eof $g]} {
		set line [split [gets $g] \t]
		foreach {nchr nbegin} [list_sub $line {3 4}] break
		if {![llength $line]} continue
		if {$nchr eq $chr && $nbegin eq $begin} {
			lappend list $line
		} else {
			while {![eof $g] && ![llength $line]} {
				set line [split [gets $g] \t]
			}
			set cache($g) $line
			break
		}
	}
	return $list
}

proc var2annotvar {file genefile outfile} {
	global cache

	catch {close $f1}
	catch {close $g}
	catch {close $o}
	unset -nocomplain cache
	set f1 [opencgifile $file header]
	set header1 [split $header \t]
	set newheader {locus chromosome begin end type reference alleleSeq1 alleleSeq2 totalScore1 totalScore2 xRef}
	set extrapos 11
	set extranum 0
	if {$header1 eq "locus haplotype chromosome begin end varType reference alleleSeq totalScore hapLink xRef"} {
		set cgv1 1.4
		set cache(cor1) {}
	} elseif {$header1 eq "locus ploidy haplotype chromosome begin end varType reference alleleSeq totalScore hapLink xRef"} {
		set cgv1 1.7
		set cache(cor1) {0 2 3 4 5 6 7 8 9 10 11}
		set cache($f1) [list_sub $cache($f1) $cache(cor1)]
	} elseif {$header1 eq "locus ploidy allele chromosome begin end varType reference alleleSeq totalScore hapLink xRef"} {
		set cgv1 1.8
		set cache(cor1) {0 2 3 4 5 6 7 8 9 10 11}
		set cache($f1) [list_sub $cache($f1) $cache(cor1)]
	} elseif {$header1 eq "locus ploidy allele chromosome begin end varType reference alleleSeq varScoreVAF varScoreEAF varQuality hapLink xRef"} {
		set cgv1 2.0
		set cache(cor1) {0 2 3 4 5 6 7 8 9 12 13 10 11}
		set cache($f1) [list_sub $cache($f1) $cache(cor1)]
		lappend newheader varScoreEAF1 varScoreEAF2 varQuality1 varQuality2
		set extranum 2
	} else {
		puts stderr "WARNING: unrecognised header (unsupported CG version), trying to interpret anayway"
		set cgv1 ?
		set cor1 [list_cor $header1 {locus allele chromosome begin end varType reference alleleSeq totalScore hapLink xRef}]
		if {[lindex $cor1 8] == -1} {
			set poss [list_find -regexp $header1 {[Ss]core}]
			set pos [lindex [list_lremove $poss $cor1] 0]
			if {[isint $pos]} {
				lset cor1 8 $pos
			}
		}
		if {[lsearch $cor1 -1] != -1} {
			error "Unsupported header error in $file (change in CG format?)"
		}
		set extra [list_remove [list_fill [llength $header1] 0 1] {*}$cor1]
		foreach pos $extra {
			lappend cor1 $pos
			set field [lindex $header1 $pos]
			lappend newheader ${field}1 ${field}2
			incr extranum
		}
		set cache(cor1) $cor1
		set cache($f1) [list_sub $cache($f1) $cache(cor1)]
		
		error "header error in $file (change in CG format?)"
	}
	set cache(extrapos) $extrapos
	set cache(extranum) $extranum
	set g [opencgifile $genefile header]
	set header2 [split $header \t]
	set remheader2 {}
	if {$header2 eq "index locus haplotype chromosome begin end varType reference call xRef geneId mrnaAcc proteinAcc orientation exonCategory exon codingRegionKnown aaCategory nucleotidePos proteinPos aaAnnot aaCall aaRef"} {
		set cgv2 1.4
	} elseif {$header2 eq "index locus haplotype chromosome begin end varType reference call xRef geneId mrnaAcc proteinAcc symbol orientation exonCategory exon codingRegionKnown aaCategory nucleotidePos proteinPos aaAnnot aaCall aaRef"} {
		set cgv2 1.7
	} elseif {$header2 eq "index locus allele chromosome begin end varType reference call xRef geneId mrnaAcc proteinAcc symbol orientation component componentIndex codingRegionKnown impact nucleotidePos proteinPos annotationRefSequence sampleSequence genomeRefSequence"} {
		set cgv2 1.8
		set remheader2 [lrange $header2 10 end]
	} elseif {$header2 eq "index locus allele chromosome begin end varType reference call xRef geneId mrnaAcc proteinAcc symbol orientation component componentIndex hasCodingRegion impact nucleotidePos proteinPos annotationRefSequence sampleSequence genomeRefSequence pfam"} {
		set cgv2 1.9
		set remheader2 [lrange $header2 10 end]
	} else {
		error "header error in $genefile (change in CG format?)"
	}
	lappend newheader {*}$remheader2
	set o [open $outfile w]
	puts $o [join $newheader \t]
	set cur [var2annotvar_readonevar $f1]
	set curgene [readgeneset $g]
	set gchr [lindex $curgene 0 3]
	set gbegin [lindex $curgene 0 4]
	set empty [list_fill [llength $remheader2] {}]
	set next 100000; set num 0
	while {![eof $f1]} {
		set fchr [lindex $cur 1]
		set fbegin [lindex $cur 2]
		incr num
		if {$num >= $next} {putsprogress $num; incr next 100000}
		while {![eof $g]} {
			set chrcomp [chrcomp $gchr $fchr]
			if {$chrcomp > 0} break
			if {($chrcomp == 0) && ($gbegin >= $fbegin)} break
			set curgene [readgeneset $g]
			set gchr [lindex $curgene 0 3]
			set gbegin [lindex $curgene 0 4]
		}
		set annot $empty
		if {($gchr eq $fchr) && ($gbegin == $fbegin)} {
			set fcomp [lrange $cur 3 4]
			foreach gline $curgene {
				set gcomp [lrange $gline 5 6]
				if {$fcomp eq $gcomp} {
					set annot [lrange $gline 10 end]
					break
				}
			}
		}
		puts $o [join $cur \t]\t[join $annot \t]
		set cur [var2annotvar_readonevar $f1]
#if {[lindex $cur 0] == 19411} {error STOP}
	}

	close $o
	close $g
	close $f1

}

if 0 {

	lappend auto_path ~/dev/completegenomics/lib
	package require Extral
package require Tclx
signal -restart error SIGINT
	set base /media/passport/complgen
	set base /complgen
	cd $base


	cd /complgen/projects/ftld1/ftld_e_d5945
	set file svar-ftld_e_d5945.tsv
	set genefile sgene-ftld_e_d5945.tsv
	set outfile annotvar-ftld_e_d5945.tsv
	var2annotvar $file $genefile $outfile

	set file var-test.tsv
	set genefile sgene-test.tsv
	set outfile annotvar-test.tsv
	var2annotvar $file $genefile $outfile

	set file GS102/svar-GS102.tsv
	set genefile GS102/sgene-GS102.tsv
	set outfile GS102/annotvar-GS102.tsv
	var2annotvar $file $genefile $outfile

	set file GS103/svar-GS103.tsv
	set genefile GS103/sgene-GS103.tsv
	set outfile GS103/annotvar-GS103.tsv
	var2annotvar $file $genefile $outfile

	cd /home/complgen/projects/dlb1/dlb_d_d388
	set file svar-dlb_d_d388.tsv
	set genefile sgene-dlb_d_d388.tsv
	set outfile annotvar-dlb_d_d388.tsv

	var2annotvar $file $genefile $outfile
}
