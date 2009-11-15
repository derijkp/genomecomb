proc comparepos {cur1 cur2} {
	array set trans {X 100 Y 101 M 102}
	set chr1 [lindex $cur1 1]
	set chr2 [lindex $cur2 1]
	if {$chr1 eq $chr2} {
		set pos1 [lindex $cur1 2]
		set pos2 [lindex $cur2 2]
		return [expr {$pos1-$pos2}]
	} else {
		set chr1 [get trans($chr1) $chr1]
		set chr2 [get trans($chr2) $chr2]
		return [expr {$chr1-$chr2}]
	}
}

proc sequenced {r1 cur2} {
	if {$r1 eq 1} {return 1}
	array set trans {X 100 Y 101 M 102}
	global cache
	set chr [lindex $cur2 1]
	set chr [get trans($chr) $chr]
	set pos [lindex $cur2 2]
	set line $cache($r1)
	while {[llength $line]} {
		foreach {rchr rstart rend} $line break
		set rchr [get trans($rchr) $rchr]
		if {$rchr > $chr} break
		if {$rchr == $chr} {
			if {$rstart > $pos} break
			if {($rchr == $chr) && ($pos <= $rend) && ($pos >= $rstart)} {
				set cache($r1) $line
				return 1
			}
		}
		set line [gets $r1]
	}
	set cache($r1) $line
	return 0
}

proc readonevar {f} {
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {[llength $line]} break
	}
	return $line
}

proc compare_annot {id1 file1 regfile1 id2 file2 regfile2 output_prefix} {
	global cache
	set outfile ${output_prefix}diff.tsv
	set outmmfile ${output_prefix}mismatch.tsv
	set outsamefile ${output_prefix}same.tsv
	set outuniqfile ${output_prefix}uniqueseq.tsv
	set f1 [open $file1]
	set header [split [gets $f1] \t]
	set annotvarfiels {locus chromosome begin end type alleleSeq1 alleleSeq2 totalScore1 totalScore2 xRef geneId mrnaAcc proteinAcc orientation exonCategory exon codingRegionKnown aaCategory nucleotidePos proteinPos aaAnnot aaCall aaRef}
	if {$header ne $annotvarfiels} {
		puts stderr "header error in annot_varfile1 $file1"
		exit 1
	}
	set f2 [open $file2]
	set header [split [gets $f2] \t]
	if {$header ne $annotvarfiels} {
		puts stderr "header error in annot_varfile2 $file1"
		exit 1
	}
	if {$regfile1 ne ""} {
		set r1 [opencgifile $regfile1 header]
		if {[split $header \t] ne "chromosome begin end ploidy type"} {
			puts stderr "header error in region_file1 $file1"
			exit 1
		}
		set cache($r1) [split [gets $r1] \t]
	} else {
		puts "no region file for file1"
		set r1 1
	}
	if {$regfile2 ne ""} {
		set r2 [opencgifile $regfile2 header]
		if {[split $header \t] ne "chromosome begin end ploidy type"} {
			puts stderr "header error in region_file2 $file2"
			exit 1
		}
		set cache($r2) [split [gets $r2] \t]
	} else {
		puts "no region file for file2"
		set r2 1
	}
	
	set o [open $outfile w]
	set osame [open $outsamefile w]
	set omm [open $outmmfile w]
	set oseq [open $outuniqfile w]
	
	set cur1 [readonevar $f1]
	set cur2 [readonevar $f2]
	set num 1
	while {![eof $f1] || ![eof $f2]} {
		incr num
		if {![expr {$num % 100000}]} {puts $num}
		set d [comparepos $cur1 $cur2]
		if {$d == 0} {
			set type [list_remove [split [lindex $cur1 4] _] ref-consistent ref-inconsistent]
			lset cur1 4 $type
			set type [list_remove [split [lindex $cur2 4] _] ref-consistent ref-inconsistent]
			lset cur2 4 $type
			if {([lrange $cur1 1 4] eq [lrange $cur2 1 4]) && (
				([list_sub $cur1 {5 6}] eq [list_sub $cur2 {5 6}])
				|| ([list_sub $cur1 {6 5}] eq [list_sub $cur2 {5 6}])
				)} {
				puts $osame $id1,$id2\t[join $cur1 \t]
			} else {
				puts $omm $id1,$id2\t[join $cur1 \t]\t[join $cur2 \t]
			}
			set cur1 [readonevar $f1]
			set cur2 [readonevar $f2]
			continue
		} elseif {$d < 0} {
			while {[comparepos $cur1 $cur2] < 0} {
				set s [sequenced $r2 $cur1]
				if {$s} {
					puts $o $id1\t[join $cur1 \t]
				} else {
					puts $oseq $id1\t[join $cur1 \t]
				}
				set cur1 [readonevar $f1]
			}
			continue
		} else {
			while {[comparepos $cur1 $cur2] > 0} {
				set s [sequenced $r1 $cur2]
				if {$s} {
					puts $o $id2\t[join $cur2 \t]
				} else {
					puts $oseq $id2\t[join $cur2 \t]
				}
				set cur2 [readonevar $f2]
			}
		}
	}
	
	close $o
	close $osame
	close $oseq
	close $f1
	close $f2
}

if 0 {

cd /complgen/

set file1 /complgen/GS00102/annotvar-GS000000078-ASM.tsv
set file2 /complgen/GS00103/annotvar-GS000000079-ASM.tsv
set regfile1 /complgen/GS00102/reg-GS000000078-ASM.tsv
set regfile2 /complgen/GS00103/reg-GS000000079-ASM.tsv
set outfile /complgen/diff.tsv
set outsamefile /complgen/same.tsv
set outmmfile /complgen/missmatch.tsv
set outuniqfile /complgen/uniqueseq.tsv
set id1 1
set id2 2

lappend auto_path ~/dev/completegenomics/lib
package require Extral

catch {close $f1}
catch {close $f2}
catch {close $r1}
catch {close $r2}
catch {close $o}
unset -nocomplain cache

wc diff.tsv
906085
grep -v ref-inconsistent diff.tsv | wc
835622
grep -v ref-inconsistent diff.tsv | grep EXON | wc
6084
grep -v ref-inconsistent diff.tsv | grep EXON | grep -P '\tCOMPATIBLE\t' | wc
1912
grep -v ref-inconsistent diff.tsv | grep INTRON | wc
273158

grep -v ref-inconsistent diff.tsv | grep -P '\tsnp\t' | wc
561545 (other = 274077)
grep -v 'ref-inconsistent' diff.tsv | grep -v '\?' | wc
525874
grep -v 'ref-inconsistent' diff.tsv | grep -v '\?' |  grep -P '\tsnp\t' | wc
399994 (other = 125880)
grep -v 'ref-inconsistent' diff.tsv | grep -v '\?' |  grep -P '\tsnp\t' | grep EXON | wc
3480

same
wc same.tsv
2998478
grep -v ref-inconsistent same.tsv | wc
2963002
grep -v ref-inconsistent same.tsv | grep -P '\tsnp\t' | wc
2642820 (other 320182)
grep -v 'ref-inconsistent' same.tsv | grep -v '\?' | wc
2826136
grep -v 'ref-inconsistent' same.tsv | grep -v '\?' |  grep -P '\tsnp\t' | wc
2584269 (other = 241867)	

}
