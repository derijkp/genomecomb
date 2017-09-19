proc bam_chrs {bamfile} {
	set bamheader [exec samtools view -H $bamfile]
	list_unmerge [regexp -all -inline {SN:([^\t]+)} $bamheader] 1 result
	return $result
}

proc bam2reg_job {args} {
	upvar job_logdir job_logdir
	set mincoverage 5
	set compress 1
	set skip {}
	cg_options bam2reg args {
		-mincoverage {
			set mincoverage $value
		}
		-compress {
			set compress $value
		}
		-skip {
			set skip $value
		}
	} {bamfile mincoverage} 1 2
	set bamfile [file_absolute $bamfile]
	set pre [lindex [split $bamfile -] 0]
	set dir [file dir $bamfile]
	set file [file tail $bamfile]
	set root [join [lrange [split [file root $file] -] 1 end] -]
#	job bam2coverage-$file -deps {$bamfile} -targets {$dir/coverage-$root $dir/coverage-$root/coverage-$root.FINISHED} -vars {root} -code {
#		cg bam2coverage $dep $target/coverage-$root
#	}
	set target $dir/sreg-cov$mincoverage-$root.tsv
	if {$compress} {append target .lz4}
	job cov$mincoverage-$root -optional 1 -deps {$bamfile} -targets {$target} -vars {mincoverage compress} \
	-skip $skip -code {
		set compress [compresspipe $target]
		set temptarget [filetemp $target]
		exec cg regextract -min $mincoverage $dep {*}$compress > $temptarget
		file rename -force $temptarget $target
		cg lz4index $target
	}
	return $target
}

