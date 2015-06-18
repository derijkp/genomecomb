package require BioTcl

proc cg_liftfindchanges {args} {
	if {([llength $args] != 3)} {
		exiterror "format is: cg liftfindchanges srcgenome destgenome liftoverfile"
		exit 1
	}
	foreach {srcgenome destgenome liftoverfile} $args break
	if {[file ext $liftoverfile] eq ".chain"} {
		set useliftoverfile [file root $liftoverfile].tsv
		if {![file exists $useliftoverfile]} {
			cg liftchain2tsv $liftoverfile $useliftoverfile
		}
		set liftoverfile $useliftoverfile
	}
	set gs [genome_open $srcgenome]
	set gd [genome_open $destgenome]
	catch {close $fl}
	set f [gzopen $liftoverfile]
	set header [tsv_open $f]
	if {$header ne {chromosome begin end strand destchromosome destbegin destend deststrand}} {
		exiterror "header of file $liftoverfile should be: chromosome begin end strand destchromosome destbegin destend deststrand"
	}
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
}
