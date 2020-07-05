# job, removes when end=begin, clips chromosome, 
proc gatkworkaround_tsv2bed {reg refseq target} {
	set f [open $refseq.index]
	while {![eof $f]} {
		set line [split [gets $f] \t]
		set chr [lindex $line 0]
		set chr [chr_clip $chr]
		set maxa($chr) [lindex $line 1 1]
	}
	close $f
	set f [gzopen $reg]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 3]
	set temptarget [filetemp $target]
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
	file rename -force -- $temptarget $target
}

# job, removes when end=begin, clips chromosome
proc gatkworkaround_tsv2bed_job {file refseq} {
	upvar job_logdir job_logdir
	set target [file root $file].bed
	job [job_relfile2name tsv2bed- $file] -optional 1 -deps {$file $refseq} -targets {$target} -code {
		gatkworkaround_tsv2bed $dep $dep2 $target
	}
	return $target
}

proc tsv2bed {tsvfile bedfile {fields {}}} {
	if {$tsvfile eq ""} {set tsvfile -}
	# puts "[list ../bin/tsv2bed $tsvfile {*}$fields]"
	if {$bedfile eq ""} {
		exec tsv2bed $tsvfile {*}$fields >@ stdout
	} else {
		exec tsv2bed $tsvfile {*}$fields > $bedfile
	}
}

proc tsv2bed_job {tsvfile {bedfile {}} {fields {}}} {
	upvar job_logdir job_logdir
	if {$bedfile eq ""} {set bedfile [file root $file].bed}
	job tsv2bed-[job_relfile2name $file] -deps {$file} -targets {$bedfile} -code {
		tsv2bed $dep $target
	}
	return $bedfile
}

proc cg_tsv2bed {args} {
	set tsvfile {}
	set bedfile {}
	set fields {}
	cg_options tsv2bed args {
		-f - -fields {
			set fields $value
		}
	} {tsvfile bedfile} 0 2
	tsv2bed $tsvfile $bedfile $fields
}
