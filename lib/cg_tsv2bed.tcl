proc gatkworkaround_tsv2bed_job {file refseq} {
	upvar job_logdir job_logdir
	job tsv2bed-[file tail $file] -deps {$file $refseq.index} -targets [file root $file].bed -code {
		set f [open $dep2]
		while {![eof $f]} {
			set line [split [gets $f] \t]
			set chr [lindex $line 0]
			set chr [chr_clip $chr]
			set maxa($chr) [lindex $line 1 1]
		}
		close $f
		set f [gzopen $dep]
		set header [tsv_open $f]
		set poss [tsv_basicfields $header 3]
		set temptarget [file_tempwrite $target]
		set o [open $temptarget w]
		while {![eof $f]} {
			set line [split [gets $f] \t]
			set line [list_sub $line $poss]
			foreach {chr begin end} $line break
			set cchr [chr_clip $chr]
			if {[info exists maxa($cchr)]} {
				if {$end > $maxa($cchr)} {set end $maxa($cchr)}
			}
			if {$end == $begin} continue
			puts $o $chr\t$begin\t$end
		}
		close $o
		gzclose $f
		file rename -force $temptarget $target
	}
	return [file root $file].bed
}

proc tsv2bed {file bedfile args} {
	if {[llength $args]} {
		set chromname [list_shift args]
		if {$chromname eq ""} {
			file_write $bedfile.temp \#[join $args \t]\n
			cg select -sh /dev/null -f "$args" $file >> $bedfile.temp
		} else {
			file_write $bedfile.temp \#chrom\t[join $args \t]\n
			cg select -sh /dev/null -f "\{chrom=\"$chromname\"\} $args" $file >> $bedfile.temp
		}
	} else {
		set f [gzopen $file]
		set header [tsv_open $f]
		gzclose $f
		set poss [tsv_basicfields $header 3]
		set fields [list_sub $header $poss]
		file_write $bedfile.temp \#[join $fields \t]\n
		cg select -sh /dev/null -f $fields $file >> $bedfile.temp		
	}
	file rename -force $bedfile.temp $bedfile
}

proc tsv2bed_job {file} {
	upvar job_logdir job_logdir
	job tsv2bed-[file tail $file] -deps $file -targets [file root $file].bed -code {
		tsv2bed $dep $target
	}
	return [file root $file].bed
}

proc cg_tsv2bed {args} {
	set tsvfile {}
	set bedfile {}
	cg_options hsmetrics args {
	} 2 2 {tsvfile bedfile}
	tsv2bed $tsvfile $bedfile
}
