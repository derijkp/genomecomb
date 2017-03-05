proc bam_chrs {bamfile} {
	set bamheader [exec samtools view -H $bamfile]
	list_unmerge [regexp -all -inline {SN:([^\t]+)} $bamheader] 1 result
	return $result
}

proc bam2reg_job {bamfile {mincoverage 5} {compress 0}} {
	upvar job_logdir job_logdir
	set bamfile [file_absolute $bamfile]
	set pre [lindex [split $bamfile -] 0]
	set dir [file dir $bamfile]
	set file [file tail $bamfile]
	set root [join [lrange [split [file root $file] -] 1 end] -]
#	job bam2coverage-$root -deps $bamfile -targets {$dir/coverage-$root $dir/coverage-$root/coverage-$root.FINISHED} -vars {root} -code {
#		cg bam2coverage $dep $target/coverage-$root
#	}
	job cov$mincoverage-$root -optional 1 -deps $bamfile -targets $dir/sreg-cov$mincoverage-$root.tsv -vars {mincoverage compress} -code {
		set temptarget [filetemp $target]
		cg regextract -min $mincoverage $dep > $temptarget
		file rename -force $temptarget $target
	}
	return $dir/sreg-cov$mincoverage-$root.tsv
}

