proc assert {check message} {
	if {![uplevel expr $check]} {
		error $message
	}
}

proc readlocus {f {pos 0}} {
	global cache
	set line $cache($f)
	set locus [lindex $line $pos]
	set list [list $line]
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		if {[lindex $line $pos] eq $locus} {
			lappend list $line
		} else {
			while {![eof $f] && ![llength $line]} {
				set line [split [gets $f] \t]
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

proc readonevar f {
	global cache list
	if {[info exists cache($f,rov)]} {
		set line1 [list_shift cache($f,rov)]
		set line2 [list_shift cache($f,rov)]
		if {![llength $cache($f,rov)]} {unset cache($f,rov)}
	} else {
		set line1 {}
		set line2 {}
	}
	while {![eof $f]} {
		while {![eof $f]} {
			set list [readlocus $f]
			set type [list_remdup [list_remove [list_subindex $list 5] = ref-consistent]]
			if {$type ne ""} break
		}
		if {$type eq ""} {return {}}
		if {[llength $list] == 1} {
			set line1 [lindex $list 0]
			set line2 {}
		} elseif {[llength $list] == 2} {
			foreach {line1 line2} $list break
		} else {
			set types [list_subindex $list 5]
			set poss [list_find -regexp $types {snp|ins|del}]
			set rlist {}
			foreach lpos $poss {
				set line1 [lindex $list $lpos]
				set pos [lindex $line1 3]
				set type2 unknown
				foreach line [list_sub $list -exclude $lpos] {
					set spos [lindex $line 3]
					set epos [lindex $line 3]
					if {$pos >= $spos && $pos <= $epos} {
						set type2 [lindex $line 5]
						break
					}
				}
				set line2 $line1
				lset line2 5 $type2
				lset line2 7 [lindex $line1 6]
				lset line2 8 {}
				lappend rlist $line1 $line2
			}
			set rlist [lsort -integer -index 3 $rlist]
			set line1 [list_shift rlist]
			set line2 [list_shift rlist]
			if {[llength $rlist]} {
				set cache($f,rov) $rlist
			}
		}
		if {![llength $line1]} continue
		break
	}
	if {![llength $line1]} {
		return {}
	}
	foreach {locus haplotype chromosome begin end varType reference alleleSeq totalScore hapLink xRef} $line1 break
	foreach {locus2 haplotype2 chromosome2 begin2 end2 varType2 reference2 alleleSeq2 totalScore2 hapLink2 xRef2} $line2 break
	set type [list [get varType unkown] [get varType2 unknown]]
	set type [list_remdup [list_remove $type = ref-consistent = ref-inconsistent unknown]]
	if {$type ne ""} {
		set result [list $locus $chromosome $begin $end [join $type _] $alleleSeq [get alleleSeq2 ""] $totalScore [get totalScore2 ""] $xRef]
		return $result
	} else {
		return {}
	}
}

proc readonegene f {
	global cache list
	set list [readlocus $f 1]
	if {![llength $list]} {return {}}
	set result [list [lindex $list 0 1]]
	for {set i 10} {$i < 24} {incr i} {
		set temp [list_remdup [list_remove [list_subindex $list $i] {}]]
		if {$i == 17} {
			set temp [list_remove $temp UNKNOWN NO-CHANGE]
		}
		lappend result $temp
	}
	foreach {index locus haplotype chromosome begin end varType reference call xRef geneId mrnaAcc proteinAcc orientation exonCategory exon codingRegionKnown aaCategory nucleotidePos proteinPos aaAnnot aaCall aaRef} [lindex $list 0] break
	return $result
}

proc var2annotvar {file genefile outfile} {
	global cache
	catch {close $f1}
	catch {close $g}
	catch {close $o}
	unset -nocomplain cache
	set f1 [opencgifile $file header]
	if {[split $header \t] ne "locus haplotype chromosome begin end varType reference alleleSeq totalScore hapLink xRef"} {
		error "header error in $file"
	}
	set g [opencgifile $genefile header]
	if {[split $header \t] ne "index locus haplotype chromosome begin end varType reference call xRef geneId mrnaAcc proteinAcc orientation exonCategory exon codingRegionKnown aaCategory nucleotidePos proteinPos aaAnnot aaCall aaRef"} {
		error "header error in $genefile"
	}
	set o [open $outfile w]
	puts $o [join {locus chromosome begin end type alleleSeq1 alleleSeq2 totalScore1 totalScore2 xRef geneId mrnaAcc proteinAcc orientation exonCategory exon codingRegionKnown aaCategory nucleotidePos proteinPos aaAnnot aaCall aaRef} \t]
	set cur [readonevar $f1]
	set curgene [readonegene $g]
	set num 0
	while {![eof $f1]} {
		incr num
		if {![expr {$num%100000}]} {puts stderr $num}
		set flocus [lindex $cur 0]
		set glocus [lindex $curgene 0]
		while {![eof $g] && [isint $glocus] && ($glocus < $flocus)} {
			set curgene [readonegene $g]
			set glocus [lindex $curgene 0]
		}
		if {$glocus == $flocus} {
			puts $o [join $cur \t]\t[join [lrange $curgene 1 end] \t]
		} else {
			puts $o [join $cur \t]
		}
		set cur [readonevar $f1]
	}
	close $o
	close $g
	close $f1
}

if 0 {

# make sorted first
cd /complgen/GS00102
head -12 gene-GS000000078-ASM.tsv > sgene-GS000000078-ASM.tsv
tail +13 gene-GS000000078-ASM.tsv | sort -nk2 >> sgene-GS000000078-ASM.tsv
cd /complgen/GS00103
head -12 gene-GS000000079-ASM.tsv > sgene-GS000000079-ASM.tsv
tail +13 gene-GS000000079-ASM.tsv | sort -nk2 >> sgene-GS000000079-ASM.tsv

# checks
cd /complgen/GS00102
grep EXON /complgen/GS00102/annotvar-GS000000078-ASM.tsv | cut -f1 > temp1
grep EXON /complgen/GS00102/sgene-GS000000078-ASM.tsv | grep -v ref-consistent | grep -v = | cut -f2 | uniq > temp2
diff temp1 temp2
	
cd /complgen/GS00103
grep EXON /complgen/GS00103/annotvar-GS000000079-ASM.tsv | cut -f1 > temp1
grep EXON /complgen/GS00103/sgene-GS000000079-ASM.tsv | grep -v ref-consistent | grep -v = | cut -f2 | uniq > temp2
diff temp1 temp2
}

if 0 {
	lappend auto_path ~/dev/completegenomics/lib
	package require Extral
package require Tclx
signal -restart error SIGINT
	cd /media/passport/complgen

	#test
	head -9 GS00102/var-GS000000078-ASM.tsv > var-test.tsv
	tail -401 GS00102/var-GS000000078-ASM.tsv >> var-test.tsv

	set file var-test.tsv
	set genefile sgene-test.tsv
	set outfile annotvar-test.tsv
	var2annot $file $genefile $outfile

	set file GS00102/var-GS000000078-ASM.tsv
	set genefile GS00102/sgene-GS000000078-ASM.tsv
	set outfile GS00102/annotvar-GS000000078-ASM.tsv
	var2annot $file $genefile $outfile

	set file GS00103/var-GS000000079-ASM.tsv
	set genefile GS00103/sgene-GS000000079-ASM.tsv
	set outfile GS00103/annotvar-GS000000079-ASM.tsv
	var2annot $file $genefile $outfile
}
