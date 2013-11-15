#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require Extral

proc mcompar_samples {file} {
	set header [cg select -h $file]
	set poss [list_find -glob $header sequenced-*]
	set done {}
	foreach el [list_sub $header $poss] {
		lappend done [string range $el 10 end]
	}
	return $done
}

proc testmultitarget {target names args} {
	file delete $target.temp
	if {[file exists $target]} {
		# test if existing target is already ok
		set done [mcompar_samples $target]
		foreach name $done {
			foreach pattern $args {
				set testfile [subst $pattern]
				if {![file exists $testfile] || ([file mtime $target] < [file mtime $testfile])} {
					file rename -force $target $target.old
					break
				}
			}
		}
		if {[file exists $target]} {
			set names [list_lremove $names $done]
			if {[llength $names]} {
				file rename -force $target $target.temp
			}
		}
	}
}

proc project {args} {
	set chrs {0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 M X Y}
	set knownactions {samples compar compar_regonly bam sv clean users}
	if {([llength $args] < 3)} {
		puts stderr "format is: $::base project_file refseqdir action/option ..."
		puts stderr " - processes project according to project file"
		puts stderr " - actions: [join $knownactions {, }]"
		exit 1
	}
	foreach {projectfile refseqdir} $args break
	if {![file isdir $refseqdir]} {
		error "$refseqdir directory does not exist"
	}
	set actions [lrange $args 2 end]
	if {[llength [list_lremove $actions $knownactions]]} {
		error "unkown action(s) [join [list_lremove $actions $knownactions] ,], must be one or more of: samples compar sv clean users"
	}
	set projectfile [file normalize $projectfile]
	if {[file isdir $projectfile]} {set projectfile [lindex [glob $projectfile/*.cgprj] 0]}
	puts "Projectfile $projectfile"
	set projectdir [file dir $projectfile]
	job_logdir [file dir $projectfile]/log_jobs
	set project [file root [file tail $projectfile]]
	set resultdir [file dir [file normalize $projectfile]]
	set c [split [file_read $projectfile] \n]
	set pos [lsearch $c {}]
	array set a [list_concat [lrange $c 0 [expr {$pos-1}]]]
	if {![info exists a(build)]} {
		error "no build given in project file"
	}
	set build $a(build)
	set cdata [list_remove [lrange $c [expr {$pos+1}] end] {}]
	set poss [list_find -glob $cdata {#*}]
	set data {}
	set names {}
	list_foreach {cgdir name} [list_sub $cdata -exclude $poss] {
		lappend data $resultdir/oricg/$cgdir $name
		lappend names $name
	}
	set names [lsort $names]
	puts "data:\n$data"

	puts "Resultdir:\n $resultdir"
	file mkdir $resultdir
	cd $resultdir
	job_logdir $resultdir/log_jobs
	# samples
	# =======
	if {[inlist $actions samples] || [inlist $actions sample]} {
		foreach {cgdir sample} $data {
			set dir [file dir $cgdir]
			file mkdir $sample
			# process_sample already uses jobs
			process_sample $cgdir $sample
		}
	}
	
	# multicompar
	# ===========
	if {[inlist $actions compar] || [inlist $actions compar_regonly]} {
		file mkdir compar
		if {[inlist $actions compar_regonly]} {
			set reannot -reannotregonly
		} else {
			set reannot -reannot
		}
		set deps {}
		foreach name $names {
			lappend deps $name/fannotvar-$name.tsv $name/sreg-$name.tsv
		}
		set target compar/compar-${project}.tsv
		testmultitarget $target $names {$name/fannotvar-$name.tsv} {$name/sreg-$name.tsv}
		job multicompar -deps {$_deps} \
		-targets {$target} \
		-vars {names reannot} -code {
			if {[file exists $target.temp]} {
				set done [mcompar_samples $target.temp]
				set names [list_lremove $names $done]
			} else {
				set done {}
			}
			puts "Already done: [join $done {, }]"
			if {[llength $names]} {
				cg multicompar $reannot $target.temp {*}$names
			}
			file rename -force $target.temp $target
		}
		job compar_annotate \
		-deps {compar/compar-$project.tsv} \
		-targets {compar/annot_compar-$project.tsv} -vars {refseqdir build} -code {
			cg annotate $dep $target.temp $refseqdir/$build
			file rename -force $target.temp $target
		}
		job compar_annotate_index \
		-deps {compar/annot_compar-$project.tsv} \
		-targets {compar/annot_compar-$project.tsv.index/info.tsv} -code {
			cg index $dep
		}

		# multireg
		# --------
		set target compar/reg-${project}.tsv
		testmultitarget $target $names {$name/sreg-$name.tsv}
		job multireg -deps [lforeach name $names $name/sreg-$name.tsv] \
		-targets {$target} \
		-vars {names} -code {
			if {[file exists $target.temp]} {
				set done [mcompar_samples $target.temp]
				set names [list_lremove $names $done]
			} else {
				set done {}
			}
			puts "Already done: [join $done {, }]"
			if {[llength $names]} {
				foreach name $names {
					lappend files $name/sreg-$name.tsv
				}
				cg multireg $target.temp {*}$files
			}
			file rename -force $target.temp $target
		}
		job multireg_index \
		-deps {compar/reg-$project.tsv} \
		-targets {compar/reg-$project.tsv.index/info.tsv} -code {
			cg index $dep
		}

		# cgsv
		# ----
		set target compar/cgsv-${project}.tsv
		testmultitarget $target $names {$name/cgsv-$name.tsv}
		job cgsv_multicompar -deps [lforeach name $names $name/cgsv-$name.tsv] \
		-targets {$target} \
		-vars {names data} -code {
			puts "Checking $target"
			if {[file exists $target.temp]} {
				set done [cg select -n $target.temp]
			} else {
				set done {}
			}
			set files {}
			foreach {cgdir name} $data {
				if {![inlist $done $name] && [file exists $name/cgsv-$name.tsv]} {
					lappend files $name/cgsv-$name.tsv
				}
			}
			if {[llength $done]} {
				puts "Multicgsv already present: $done"
			}
			if {[llength $files]} {
				cg svmulticompar $target.temp {*}$files
			}
			file rename -force $target.temp $target
		}
		job cgsv_annotate \
		-deps {compar/cgsv-$project.tsv} \
		-targets {compar/annot_cgsv-$project.tsv} -vars {refseqdir build} -code {
			cg annotate $dep $target.temp $refseqdir/$build
			file rename -force $target.temp $target
		}
		job cgsv_annotate_index \
		-deps {compar/annot_cgsv-$project.tsv} \
		-targets {compar/annot_cgsv-$project.tsv.index/info.tsv} -code {
			cg index $dep
		}

		# cgcnv
		# ----
		set target compar/cgcnv-${project}.tsv
		testmultitarget $target $names {$name/cgcnv-$name.tsv}
		job cgcnv_multicompar -deps [lforeach name $names $name/cgcnv-$name.tsv] \
		-targets {compar/cgcnv-${project}.tsv} \
		-vars {names data} -code {
			puts "Checking $target"
			if {[file exists $target.temp]} {
				set done [cg select -n $target.temp]
			} else {
				set done {}
			}
			set files {}
			foreach {cgdir name} $data {
				if {![inlist $done $name] && [file exists $name/cgcnv-$name.tsv]} {
					lappend files $name/cgcnv-$name.tsv
				}
			}
			if {[llength $done]} {
				puts "Multicgcnv already present: $done"
			}
			if {[llength $files]} {
				cg svmulticompar $target.temp {*}$files
			}
			file rename -force $target.temp $target
		}
		job cgcnv_annotate \
		-deps {compar/cgcnv-$project.tsv} \
		-targets {compar/annot_cgcnv-$project.tsv} -vars {refseqdir build} -code {
			cg annotate $dep $target.temp $refseqdir/$build
			file rename -force $target.temp $target
		}
		job cgcnv_annotate_index \
		-deps {compar/annot_cgcnv-$project.tsv} \
		-targets {compar/annot_cgcnv-$project.tsv.index/info.tsv} -code {
			cg index $dep
		}
	}
	# make bam files
	# ==================
	if {[inlist $actions bam]} {
		# untested ! 
		cd $resultdir
		foreach {cgdir name} $data {
			file mkdir $resultdir/$name/bam
			set destprefix $resultdir/$name/bam/bam_$name
			job cg2bam-$name -deps {$cgdir} -targets {$destprefix-chr$_chrs.bam} -vars {cgdir destprefix refseqdir build} {
				cg cg2bam $cgdir $destprefix $refseqdir/$build
			}
		}
	}

	# sv (our algorithm)
	# ==================
	if {[inlist $actions sv]} {
		# untested
		cd $resultdir
		set svjobs {}
		foreach {cgdir name} $data {
			set dir [file dir $cgdir]
			job process_sv-$name -deps {$cgdir $cgdir/MAP} -targets {$name/sv-$name.tsv $name/sv} -vars {name refseqdir build} -code {
				cg process_sv $dep $name $refseqdir/$build
			}
		}
		#
		# sv compar
		# ---------
		set target compar/sv-${project}.tsv
		testmultitarget $target $names {$name/sv-$name.tsv}
		job sv_multicompar -deps [lforeach name $names $name/sv-$name.tsv] \
		-targets {compar/sv-${project}.tsv} \
		-vars {names data} -code {
			puts "Checking $target"
			if {[file exists $target.temp]} {
				set done [cg select -n $target.temp]
			} else {
				set done {}
			}
			set files {}
			foreach {cgdir name} $data {
				if {![inlist $done $name] && [file exists $name/sv-$name.tsv]} {
					lappend files $name/sv-$name.tsv
				}
			}
			if {[llength $done]} {
				puts "Multisv already present: $done"
			}
			if {[llength $files]} {
				cg svmulticompar $target.temp {*}$files
			}
			file rename -force $target.temp $target
		}
		job sv_annotate \
		-deps {compar/sv-$project.tsv} \
		-targets {compar/annot_sv-$project.tsv} -vars {refseqdir build} -code {
			cg annotate $dep $target.temp $refseqdir/$build
			file rename -force $target.temp $target
		}
		job sv_annotate_index \
		-deps {compar/annot_sv-$project.tsv} \
		-targets {compar/annot_sv-$project.tsv.index/info.tsv} -code {
			cg index $dep
		}
	}

	# cleanup/compress
	# ================
	if {[inlist $actions clean]} {
		cd $resultdir
		foreach {cgdir name} $data {
			puts "Cleanup $name"
			# set host [lindex [file split $cgdir] 2]
			# remove temp files
			set files [glob -nocomplain $name/*.temp $name/*/*.temp $name/*/*/*.temp $name/*.old $name/*/*.old]
			if {[llength $files]} {
				puts "Deleting $files"
				file delete {*}$files
			}
			# compress sample results
			set files [glob -nocomplain $name/*.tsv $name/*/*.tsv $name/*/*/*.tsv]
			foreach file $files {
				puts "Compressing $file"
				set job [submit -deps $alljobs bgzip $file]
			}
			# compress sv data
			set dir [file dir $cgdir]
			set files [glob -nocomplain $name/sv/*-paired.tsv]
			foreach file $files {
				puts "Compressing $file"
				set job [submit -deps $alljobs bgzip $file]
			}
			set files [glob -nocomplain compar/*.temp compar/*.old]
			if {[llength $files]} {
				puts "Deleting $files"
				file delete {*}$files
			}
		}
	}

	# users
	# =====
	if {[inlist $actions users]} {
		cd $resultdir
		exec chmod g-w {*}[glob /complgen/projects/*]
		foreach user [get a(users) ""] {
			puts "$project -> $user"
			if {![file exists /home/MOLGEN/$user/complgen]} {
				catch {file mkdir /home/MOLGEN/$user/complgen}
				exec chown "$user.domain users" /home/MOLGEN/$user
				exec chown "$user.domain users" /home/MOLGEN/$user/complgen
			}
			cd /home/MOLGEN/$user/complgen
			if {![file exists docs]} {	
				exec ln -s /complgen/projects/docs .
			}
			if {![catch {file link $destdir}]} {
				file delete $destdir
			}
			cg_cplinked $projectdir .
			exec chown -R "$user.domain users" .
		}
	}
}

proc cg_project {args} {
	set args [job_init {*}$args]
	project {*}$args
	job_wait
}

if 0 {
	foreach {cgdir name} $data {
		set dir [file dir $cgdir]
		cg select -q {$size > 150} $name/sv/$name-1-paired-sv.tsv $name/sv-$name.tsv
		foreach chr {2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 M X Y} {
			catch {file delete temp}
			catch {
				cg select -q {$size > 150} $name/sv/$name-$chr-paired-sv.tsv temp
				exec tail -n +2 temp >> $name/sv-$name.tsv
			} e
			puts $e
		}
	}
	set names {}
	foreach {cgdir name} $data {
		lappend names $name
	}
	file mkdir compar
	if {[info exists svjobs] && [llength $svjobs]} {
		set temp [eval {exec submit.tcl -deps [join $svjobs ,] cg svmulticompar compar/${project}_svcompar.tsv} $names]
	} else {
		set temp [eval {exec submit.tcl cg svmulticompar compar/compar-${project}.tsv} $names]
	}

}
