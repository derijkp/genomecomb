proc process_multicompar_job {experiment dbdir todo args} {
	upvar job_logdir job_logdir
	set skipincomplete 1
	set split 0
	set dbfiles {}
	set addtargets 0
	foreach {key value} $args {
		if {$key eq "-skipincomplete"} {
			set skipincomplete $value
		} elseif {$key eq "-split"} {
			set split $value
		} elseif {$key eq "-dbfiles"} {
			set dbfiles $value
		} elseif {$key eq "-targetsfile"} {
			set addtargets 1
			set targetsfile $value
		} else {
			lappend opts $key $value
		}
	}
	file mkdir compar
	#
	# multicompar
	if {[file exists samples]} {set sampledir samples/} else {set sampledir {}}
	if {[catch {cg select -n compar/compar-$experiment.tsv} done]} {set done {}}
	set done [split $done \n]
	set stilltodo {}
	set deps {}
	foreach sample $todo {
		set name [lindex [split $sample -] end]
		if {![inlist $done $sample]} {
			lappend stilltodo $sampledir$name/var-$sample.tsv
			lappend deps \($sampledir$name/sreg-$sample.tsv\) \($sampledir$name/varall-$sample.tsv\)
			lappend deps \($sampledir$name/coverage/coverage-*.bcol\) \($sampledir$name/coverage/refScore-*.bcol\)
			lappend deps \($sampledir$name/coverage/coverage-*.tsv\)
			lappend deps \($sampledir$name/reg_refcons-$sample.tsv\) \($sampledir$name/reg_nocall-$sample.tsv\) \($sampledir$name/reg_cluster-$sample.tsv\)
		}
	}
	if {$addtargets} {
		if {[catch {cg select -n compar/compar-$experiment.tsv} header]} {set header {}}
		if {![llength $stilltodo] && [inlist $header [lindex [split $targetsfile -] end]]} {
			set addtargets 0
		} else {
			lappend deps $targetsfile
		}
	}
	if {[llength $stilltodo] || $addtargets} {
		file delete compar/compar-$experiment.tsv.temp
		if {[file exists compar/compar-$experiment.tsv]} {
			file rename -force compar/compar-$experiment.tsv compar/compar-$experiment.tsv.temp
		}
		job multicompar-$experiment -deps [list_concat $stilltodo $deps] -targets compar/compar-$experiment.tsv \
		-vars {stilltodo skipincomplete split addtargets targetsfile} -code {
			# should maybe better recheck todo here
			if {$addtargets} {
				cg multicompar -split $split -targetsfile $targetsfile $target.temp {*}$stilltodo
			} else {
				cg multicompar -split $split $target.temp {*}$stilltodo
			}
			if {$skipincomplete} {
				cg multicompar_reannot -paged 100 $target.temp skipincomplete
			} else {
				cg multicompar_reannot -paged 100 $target.temp
			}
			file rename -force $target.temp $target
		}
	}
	job annotcompar-$experiment -deps [list compar/compar-$experiment.tsv {*}$dbfiles] \
	-targets compar/annot_compar-$experiment.tsv -vars {dbdir dbfiles} -code {
		cg annotate $dep $target.temp $dbdir {*}$dbfiles
		file rename -force $target.temp $target
	}
	job indexannotcompar-$experiment \
	-deps compar/annot_compar-$experiment.tsv \
	-targets compar/annot_compar-$experiment.tsv.index/info.tsv -vars dbdir -code {
		cg index -colinfo $dep
	}
	#
	# multi sreg
	if {[catch {cg select -n compar/sreg-$experiment.tsv} done]} {set done {}}
	set stilltodo {}
	foreach sample $todo {
		set name [lindex [split $sample -] end]
		if {![inlist $done $sample]} {
			lappend stilltodo \($sampledir$name/sreg-$sample.tsv\)
		}
	}
	if {[llength $stilltodo]} {
		file delete compar/sreg-$experiment.tsv.temp
		if {[file exists compar/sreg-$experiment.tsv]} {
			file rename -force compar/sreg-$experiment.tsv compar/sreg-$experiment.tsv.temp
		}
		job sreg-$experiment -deps $stilltodo -targets compar/sreg-$experiment.tsv -vars stilltodo -code {
			cg multireg $target.temp {*}[list_remove $deps {}]
			file rename -force $target.temp $target
		}
		job sreg-index-$experiment -deps compar/sreg-$experiment.tsv -targets compar/sreg-$experiment.tsv.index -code {
			cg index $dep
		}
	}
	#
	# cgsv
	# ----
	set files [jobglob $sampledir*/cgsv-*.tsv]
	if {[llength $files]} {
		set target compar/cgsv-${experiment}.tsv
		set names [list_regsub {.*/cgsv-(.*)\.tsv.*} $files {\1}]
		testmultitarget $target $names "$sampledir\$name/cgsv-\$name.tsv"
		job cgsv_multicompar -deps $files -targets {$target} -code {
			puts "Checking $target"
			if {[file exists $target.temp]} {
				set done [cg select -n $target.temp]
			} else {
				set done {}
			}
			set todo {}
			foreach file $deps {
				regexp {cgsv-(.*)\.tsv*} [file tail $file] temp name
				if {![inlist $done $name]} {
					lappend todo $file
				}
			}
			if {[llength $done]} {
				puts "Multicgsv already present: $done"
			}
			if {[llength $todo]} {
				cg svmulticompar $target.temp {*}$todo
			}
			file rename -force $target.temp $target
		}
		job cgsv_annotate \
		-deps {compar/cgsv-$experiment.tsv} \
		-targets {compar/annot_cgsv-$experiment.tsv} -vars {refseqdir build} -code {
			cg annotate $dep $target.temp $refseqdir/$build
			file rename -force $target.temp $target
		}
		job cgsv_annotate_index \
		-deps {compar/annot_cgsv-$experiment.tsv} \
		-targets {compar/annot_cgsv-$experiment.tsv.index/info.tsv} -code {
			cg index $dep
		}
	}
	# cgcnv
	# ----
	set files [jobglob $sampledir*/cgsv-*.tsv]
	if {[llength $files]} {
		set target compar/cgcnv-${experiment}.tsv
		set names [list_regsub {.*/cgsv-(.*)\.tsv.*} $files {\1}]
		testmultitarget $target $names "$sampledir\$name/cgcnv-\$name.tsv"
		job cgcnv_multicompar -deps $files -targets {compar/cgcnv-${experiment}.tsv} -code {
			puts "Checking $target"
			if {[file exists $target.temp]} {
				set done [cg select -n $target.temp]
			} else {
				set done {}
			}
			set todo {}
			foreach file $deps {
				regexp {cgsv-(.*)\.tsv*} [file tail $file] temp name
				if {![inlist $done $name]} {
					lappend todo $file
				}
			}
			if {[llength $done]} {
				puts "Multicgcnv already present: $done"
			}
			if {[llength $todo]} {
				cg svmulticompar $target.temp {*}$todo
			}
			file rename -force $target.temp $target
		}
		job cgcnv_annotate \
		-deps {compar/cgcnv-$experiment.tsv} \
		-targets {compar/annot_cgcnv-$experiment.tsv} -vars {refseqdir build} -code {
			cg annotate $dep $target.temp $refseqdir/$build
			file rename -force $target.temp $target
		}
		job cgcnv_annotate_index \
		-deps {compar/annot_cgcnv-$experiment.tsv} \
		-targets {compar/annot_cgcnv-$experiment.tsv.index/info.tsv} -code {
			cg index -colinfo $dep
		}
	}
}

proc process_multicompar {args} {
	set dbdir {}
	set dbfiles {}
	set realign 1
	set paired 1
	set adapterfile {}
	set conv_nextseq 0
	cg_options process_multicompar args {
		-dbdir {
			set dbdir $value
		}
		-split {
			set split $value
		}
		-dbfile {
			lappend dbfiles $value
		}
	} 1 2
	set len [llength $args]
	if {$len == 1} {
		set destdir [lindex $args 0]
	} elseif {$len == 2} {
		foreach {destdir dbdir} $args break
	}
	set destdir [file_absolute $destdir]
	# check projectinfo
	projectinfo $destdir dbdir split
	set refseq [glob $dbdir/genome_*.ifas]
	set samples {}
	set experiment [file tail $destdir]
	if {[file exists $destdir/samples]} {
		set sampledir $destdir/samples
	} else {
		set sampledir $destdir
	}
	set todo {}
	foreach file [jobglob ${sampledir}/*/var-*.tsv] {
		lappend todo [string range [file root [file tail [gzroot $file]]] 4 end]
	}
	job_logdir $destdir/log_jobs
	set keeppwd [pwd]
	cd $destdir
	set todo [list_remdup $todo]
	process_multicompar_job $experiment $dbdir $todo -skipincomplete 1 -split $split -dbfiles $dbfiles
	cd $keeppwd

}

proc cg_process_multicompar {args} {
	set args [job_init {*}$args]
	if {[llength $args] < 1} {
		errorformat process_multicompar
		exit 1
	}
	process_multicompar {*}$args
	job_wait
}
