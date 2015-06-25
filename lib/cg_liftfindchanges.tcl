package require BioTcl

proc cg_liftfindchanges {args} {
	if {([llength $args] != 3)} {
		exiterror "format is: cg liftfindchanges srcgenome destgenome liftoverfile"
		exit 1
	}
	foreach {srcgenome destgenome liftoverfile} $args break
	set gs [genome_open $srcgenome]
	set gd [genome_open $destgenome]
	catch {close $fl}
	set f [openliftoverfile $liftoverfile header oldref newref]
	puts "\#filetype\ttsv/liftover_refchanges"
	puts "\#fileversion\t[fileversion]"
	puts [join {chromosome begin end ref destchromosome destbegin destend destref destcomplement} \t]
	while 1 {
		if {[gets $f line] == -1} break
		set line [split $line \t]
		foreach {chromosome begin end strand destchromosome destbegin destend deststrand} $line break
		set chromosome [chr_clip $chromosome]
		set destchromosome [chr_clip $destchromosome]
		if {[catch {
			set sseq [string toupper [genome_get $gs $chromosome $begin $end]]
		} msg]} {
			continue
		}
		if {[catch {
			set dseq [string toupper [genome_get $gd $destchromosome $destbegin $destend]]
		} msg]} {
			continue
		}
		if {$strand ne "+"} {exiterror "source strand should allways be +"}
		if {$deststrand eq "-"} {
			set rev 1
			set dseq [seq_complement $dseq]
		} else {
			set rev 0
		}
		if {$sseq ne $dseq} {
 			if {[expr {($end-$begin)-($destend-$destbegin)}]} {
				error "difference in length between matching blocks not supported"
			}
			if {!$rev} {
				set step 1
			} else {
				set destbegin [expr {$destend - 1}]
				set step -1
			}
			foreach s [split $sseq ""] d [split $dseq ""] {
				if {$s ne $d} {
					if {$rev} {set d [seq_complement $d]}
					puts [join [list $chromosome $begin [expr {$begin+1}] $s $destchromosome $destbegin [expr {$destbegin+1}] $d $rev] \t]
				}
				incr begin ; incr destbegin $step
			}
		}
	}	
	genome_close $gs
	genome_close $gd
	closeliftoverfile $f
}
