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

proc var2annotvar_readonevar f {
	global cache list
	# join $cache($f,rov) \n
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
			if {$type eq ""} {return {}}
			if {[llength $list] == 1} {
				set line1 [lindex $list 0]
				set line2 {}
			} else {
				set list [lsort -integer -index 3 [lsort -integer -index 4 $list]]
				set keeplist $list
				set poss [list_find -regexp $types {snp|ins|del|delins|sub}]
				set list [list_sub $list $poss]
				set rlist {}
				set line1 [list_shift list]
				set comp1 [lrange $line1 3 5]
				while {[llength $list] || [llength $line1]} {
					set line [list_shift list]
					set comp [lrange $line 3 5]
					if {$comp eq $comp1} {
						lappend rlist $line1 $line
						set line1 [list_shift list]
						set comp1 [lrange $line1 3 5]
					} else {
						set reference [lindex $line1 6]
						set line2 {}
						foreach templine [list_remove $keeplist $line1] {
							foreach {begin end} [lrange $templine 3 4] break
							foreach {b1 e1} $comp1 break
							set overlap [overlap $b1 $e1 $begin $end]
							if {$overlap > 0} {
								set line2 $templine
								break
							}
						}
						if {[llength $line2]} {
							# only change things in line2 that are actually used
							set vartype [lindex $line2 5]
							if {($vartype eq "del") || ($vartype eq "delins")} {
								lset line2 7 {}
							} elseif {($vartype eq "=") || ([inlist {ref-consistent no-call-rc} $vartype])} {
								lset line2 7 $reference
							} elseif {[inlist {ref-inconsistent no-call-ri} $vartype]} {
								lset line2 7 N
							} else {
								set temp [string range [lindex $line2 7] [expr {$begin-$b1}] [expr {$end-$e1}]]
								if {[string length $temp] == [expr {$end-$begin}]} {
									lset line2 7 $temp
								}
							}
							lset line2 5 [lindex $line1 5]
						}
						lappend rlist $line1 $line2
						set line1 $line
						set comp1 $comp
					}
				}
				set line1 [list_shift rlist]
				set line2 [list_shift rlist]
				if {[llength $rlist]} {
					set cache($f,rov) $rlist
				}
			}
			if {![llength $line1]} continue
			break
		}
	}
	if {![llength $line1]} {
		return {}
	}
	foreach {locus haplotype chromosome begin end varType reference alleleSeq totalScore hapLink xRef} $line1 break
	foreach {locus2 haplotype2 chromosome2 begin2 end2 varType2 reference2 alleleSeq2 totalScore2 hapLink2 xRef2} $line2 break
	set type [list [get varType unkown] [get varType2 unknown]]
	set type [list_remdup [list_remove $type = ref-consistent no-call-rc ref-inconsistent no-call-ri ref no-ref no-call no-ref PAR-called-in-X unknown]]
	if {$type ne ""} {
		set alleleSeq2 [get alleleSeq2 $reference]
		set totalScore2 [get totalScore2 ""]
		if {$alleleSeq2 < $alleleSeq} {
			set temp $alleleSeq ; set alleleSeq $alleleSeq2 ; set alleleSeq2 $temp
			set temp $totalScore ; set totalScore $totalScore2 ; set totalScore2 $temp
		}
		if {$alleleSeq eq "?"} {set alleleSeq -}
		if {$alleleSeq2 eq "?"} {set alleleSeq2 -}
		set result [list $locus $chromosome $begin $end [join $type _] $reference $alleleSeq $alleleSeq2 $totalScore $totalScore2 $xRef]
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
		if {![expr {$num%10000}]} {putslog $num}
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
	set header [split [gets $f] \t]
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
		if {![expr $num%10000]} {putslog $num}
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
	} else {
		error "header error in $file"
	}
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
		error "header error in $genefile"
	}
	set o [open $outfile w]
	set newheader {locus chromosome begin end type reference alleleSeq1 alleleSeq2 totalScore1 totalScore2 xRef}
	lappend newheader {*}$remheader2
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
		if {$num >= $next} {putslog $num; incr next 100000}
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
		if {[lindex $cur 4] eq "sub" && [string length [lindex $cur 5]] <= 3} {
			foreach {bin chr start end type ref a1 a2 s1 s2 info} $cur break
			set len [string length $ref]
			if {$a1 eq "?"} {
				set a1 [string_fill ? $len]
			}
			if {$a2 eq "?"} {
				set a2 [string_fill ? $len]
			}
			if {[string length $a1] != $len ||[string length $a2] != $len} {
				puts $o [join $cur \t]\t[join $annot \t]
			} else {
				string_foreach r $ref e1 $a1 e2 $a2 {
					if {$r ne $e1 || $r ne $e2} {
						puts $o [join [list $bin $chr $start [expr {$start+1}] snp $r $e1 $e2 $s1 $s2 $info] \t]\t[join $annot \t]
					}
					incr start
				}
			}
		} else {
			puts $o [join $cur \t]\t[join $annot \t]
		}
		set cur [var2annotvar_readonevar $f1]
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
