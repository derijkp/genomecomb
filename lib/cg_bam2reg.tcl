proc bam2reg_job {args} {
	upvar job_logdir job_logdir
	set compress 1
	set skips {}
	set distrreg 0
	cg_options bam2reg args {
		-mincoverage {
			set mincov $value
		}
		-compress {
			set compress $value
		}
		-distrreg {
			set distrreg [distrreg_checkvalue $value]
		}
		-refseq {
			set refseq $value
		}
		-skip {
			lappend skips -skip $value
		}
	} {bamfile mincoverage target} 1 3
	set bamfile [file_absolute $bamfile]
	set pre [lindex [split $bamfile -] 0]
	set dir [file dir $bamfile]
	set file [file tail $bamfile]
	set root [file_rootname $file]
	if {![info exists target] && [info exists mincoverage] && [info exists mincov]} {
		set target $mincoverage
		set mincoverage $mincov
	}
	if {![info exists mincoverage]} {
		if {[info exists mincov]} {set mincoverage $mincov} else {set mincoverage 5}
	}
	if {![info exists target]} {
		set target $dir/sreg-cov$mincoverage-$root.tsv
		if {$compress} {append target .zst}
	}
	if {![info exists job_logdir]} {
		set_job_logdir [file dir $target]/log_jobs
	}
	bam_index_job {*}$skips $bamfile
	set bamindexfile [index_file $bamfile]
	if {$distrreg in {0 {}}} {
		job cov$mincoverage-$root -optional 1 {*}$skips -deps {
			$bamfile ($bamindexfile)
		} -targets {
			$target
		} -vars {
			bamfile mincoverage refseq
		} -code {
			analysisinfo_write $bamfile $target regextract genomecomb regextract_version [version genomecomb] regextrac_samtools [version samtools]
			set compress [compresspipe $target]
			set temptarget [filetemp $target]
			cg regextract -refseq $refseq -min $mincoverage $dep {*}$compress > $temptarget
			file rename -force -- $temptarget $target
			if {[file extension $target] eq ".zst"} {zstindex $target}
		}
	} else {
		if {![info exists refseq]} {
			error "-distrreg cannot be used without giving a refseq using the -refseq option"
		}
		set regions [distrreg_regs $distrreg $refseq]
		set workdir [workdir $target]
		set todo {}
		foreach region $regions {
			set subtarget $workdir/sreg-cov$mincoverage-$root-$region.tsv.zst
			lappend todo $subtarget
			job cov$mincoverage-$root-$region -optional 1 {*}$skips -skip $target -deps {
				$bamfile ($bamindexfile)
			} -targets {
				$subtarget
			} -vars {
				bamfile mincoverage subtarget region refseq
			} -code {
				set bamheader [catch_exec samtools view --no-PG -H $bamfile]
				set chr [lindex [split $region :-] 0]
				if {![regexp \tSN:$chr\t $bamheader]} {
					file_write $target ""
				} else {
					set compress [compresspipe $subtarget 1]
					set temptarget [filetemp $subtarget]
					cg regextract -stack 1 -region $region -refseq $refseq -min $mincoverage $bamfile {*}$compress > $temptarget
					file rename -force -- $temptarget $subtarget
				}
			}
		}
		job cov$mincoverage-$root-merge -optional 1 {*}$skips -deps $todo -targets {
			$target
		} -vars {
			bamfile
		} -code {
			analysisinfo_write $bamfile $target regextract genomecomb regextract_version [version genomecomb] regextrac_samtools [version samtools]
			set compress [compresspipe $target]
			cg cat -c 0 {*}$deps {*}$compress > $target.temp
			file rename -- $target.temp $target
			if {[file extension $target] eq ".zst"} {zstindex $target}
			foreach file $deps {
				file delete $file [analysisinfo_file $file] [index_file $file]
			}
			catch {file delete [file dir $file]}
		}
	}
	return $target
}

proc cg_bam2reg {args} {
	set args [job_init {*}$args]
	unset job_logdir
	set result [bam2reg_job {*}$args]
	job_wait
	return $result
}
