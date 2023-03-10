#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

if 0 {
	set genome /complgen/refseq/hg19/genome_hg19
	set file /projects/vegdep_data/seq_vegdep.fas
	set resultfile /projects/vegdep/reg_vegdep.map
}

proc make_remap {file genome resultfile} {
	set tempfile [tempfile]
	exec lastal $genome $file > $tempfile
	
	catch {close $o}
	set o [open $resultfile.temp w]
	puts $o [join {name begin end destname destbegin destend} \t]
	set f [open $tempfile]
	while {![eof $f]} {
		set line [gets $f]
		if {[string index $line 0] ne "#"} break
	}
	unset -nocomplain done
	while {![eof $f]} {
		set score [lindex [split $line =] end]
		if {$score < 100} break
		set line1 [gets $f]
		foreach {temp destname deststart destsize deststrand temp destseq} $line1 break
		set destend [expr {$deststart+$destsize}]
		set line2 [gets $f]
		set line1 [lrange $line1 0 end-1]
		set line2 [lrange $line2 0 end-1]
		foreach {temp name start size strand temp seq} $line2 break
		set end [expr {$start+$size}]
		set line [gets $f]
		set line [gets $f]
		if {[regexp -- - $seq] || [regexp -- - $destseq]} continue
		# if {$start == 29190} {error stop}
		set found 0
		list_foreach {b e} [get done($name) ""] {
			if {$b < $end && $e > $start} {
				if {$b > $start} {
					set end $b
					set destend [expr {$deststart+($end-$start)}]
				} elseif {$e < $end} {
					set deststart [expr {$deststart+($e-$start)}]
					set start $e
				} else {
					# puts "Skipping: $name:$b-$e overlaps with $name:$start-$end"
					set found 1; break
				}
			}
		}
		if {$found} continue
		set res [list $name $start $end $destname $deststart $destend]
		puts $o [join $res \t]
		lappend done($name) [list $start $end]
	}
	
	close $o
	catch {exec cg select -s {name begin end}  $resultfile.temp $resultfile.temp2}
	file rename -force -- $resultfile.temp2 $resultfile
	file delete $resultfile.temp
}

proc cg_remap {file remapfile resultfile} {
	catch {close $f}; catch {close $fr}; catch {close $o}
	set f [gzopen $file]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 3]
	foreach {chrpos bpos epos} $poss break
	set fr [gzopen $remapfile]
	set rheader [tsv_open $fr]
	if {
		[lrange $rheader 0 5] ne "name begin end destname destbegin destend"
		&& [lrange $rheader 0 5] ne "chromosome begin end destchromosome destbegin destend"
	} {
		error "remap file $remapfile does not have the required header (name begin end destname destbegin destend) or (chromosome begin end destchromosome destbegin destend)"
	}
	set rline [split [gets $fr] \t]
	foreach {rchr rbegin rend rdestchr rdestbegin rdestend} $rline break
	set o [open $resultfile.temp w]
	set u [open $resultfile.temp.unmapped w]
	puts $u [join $header \t]
	puts $o [join $header \t]\torichr\toribegin\toriend
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set cur [list_sub $line $poss]
		foreach {chr begin end} $cur break
		while 1 {
			set chrcomp [chr_compare $rchr $chr]
			if {$chrcomp > 0} break
			if {($chrcomp == 0) && ($rend > $begin)} break
			if {[eof $fr]} break
			set rline [split [gets $fr] \t]
			foreach {rchr rbegin rend rdestchr rdestbegin rdestend} $rline break
		}
		if {$chrcomp != 0 || $begin < $rbegin || $end > $rend} {
			putslog "cannot remap $line: outside of regions in remap file"
			puts $u [join $line \t]
		} else {
			lset line $chrpos $rdestchr
			set shift [expr {$rdestbegin-$rbegin}]
			lset line $bpos [expr {$begin + $shift}]
			lset line $epos [expr {$end + $shift}]
			lappend line $chr $begin $end
			puts $o [join $line \t]
		}
	}
	close $o
	close $u
	gzclose $f ; gzclose $fr
	cg select -s - $resultfile.temp $resultfile.temp2
	file delete $resultfile.temp
	file rename -force -- $resultfile.temp.unmapped $resultfile.unmapped
	file rename -force -- $resultfile.temp2 $resultfile
}
