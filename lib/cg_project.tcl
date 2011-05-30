#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require Extral

proc submit {args} {
	set host {}
	set deps {}
	set submitargs {}
	set direct 0
	while 1 {
		set key [lindex $args 0]
		if {$key eq "-host"} {
			set host [lindex $args 1]
			if {[llength $host]} {
				lappend submitargs -host $host
			}
			set args [lrange  $args 2 end]
		} elseif {$key eq "-deps"} {
			set deps [lindex $args 1]
			if {[llength $deps]} {
				lappend submitargs -deps [join $deps ,]
			}
			set args [lrange  $args 2 end]
		} elseif {$key eq "-direct"} {
			set direct [lindex $args 1]
			set args [lrange  $args 2 end]
		} else {
			break
		}
	}
	if {$direct} {
		puts "Processing: $args"
		exec {*}$args >@ stdout 2>@stderr
		return ""
	}
	puts "Submitting: $submitargs $args"
	set error [catch {
		set temp [eval exec submit.tcl $submitargs $args]
	} result]
	if {$error} {
		if {[regexp {is running, skipping} $result]} {
			puts $result
			return ""
		} else {
			error $result $::errorInfo
		}
	}
	puts $result
	if {[regexp {[0-9]+} $temp job]} {
		return $job
	} else {
		return ""
	}
}

proc mcompar_samples {file} {
	set header [cg select -h $file]
	set poss [list_find -glob $header sequenced-*]
	set done {}
	foreach el [list_sub $header $poss] {
		lappend done [string range $el 10 end]
	}
	return $done
}

proc cg_project {args} {
	set knownactions {samples compar compar_regonly sv clean users}
	if {([llength $args] < 3)} {
		puts stderr "format is: $::base project_file refseqdir action/option ..."
		puts stderr " - processes project according to project file"
		puts stderr " - actions: [join $knownactions {, }]"
		puts stderr " - options: direct"
		exit 1
	}
	foreach {projectfile refseqdir} $args break
	if {![file isdir $refseqdir]} {
		error "$refseqdir directory does not exist"
	}
	set actions [lrange $args 2 end]
	set pos [lsearch $actions direct]
	set direct 0
	if {$pos != -1} {
		set actions [list_remove $actions direct]
		set direct 1
	}
	if {[llength [list_lremove $actions $knownactions]]} {
		error "unkown action(s) [join [list_lremove $actions $knownactions] ,], must be one or more of: samples compar sv clean users"
	}
	set projectfile [file normalize $projectfile]
	if {[file isdir $projectfile]} {set projectfile [lindex [glob $projectfile/*.cgprj] 0]}
	puts "Projectfile $projectfile"
	set projectdir [file dir $projectfile]
	set project [file root [file tail $projectfile]]
	set resultdir [file dir [file normalize $projectfile]]
	set c [split [file_read $projectfile] \n]
	set pos [lsearch $c {}]
	array set a [list_concat [lrange $c 0 [expr {$pos-1}]]]
	set build [get a(build) hg18]
	set data [list_remove [lrange $c [expr {$pos+1}] end] {}]
	set poss [list_find -glob $data {#*}]
	set data [join [list_sub $data -exclude $poss] \n]
	puts "data:\n$data"

	# samples
	# =======
	puts "Resultdir: $resultdir"
	file mkdir $resultdir
	cd $resultdir
	set alljobs {}
	set jobs {}
	if {[inlist $actions samples] || [inlist $actions sample]} {
		foreach {cgdir name} $data {
			set dir [file dir $cgdir]
			file mkdir $name
			set host [lindex [file split $cgdir] 2]
			set job [submit -direct $direct -host $host cg process_sample $cgdir $name $refseqdir/$build]
			if {[isint $job]} {lappend jobs $job}
		}
		lappend alljobs {*}$jobs
	}
#	cd $resultdir
#	foreach {cgdir name} $data {
#		set dir [file dir $cgdir]
#		catch {exec ln -s $dir/$name .}
#		file mkdir ori
#		catch {exec ln -s $cgdir ori}
#		catch {exec ln -s $cgdir ori/$name.ori}
#	}
	
	# multicompar
	# ===========
	if {[inlist $actions compar] || [inlist $actions compar_regonly]} {
		if {[inlist $actions compar_regonly]} {
			set reannot -reannotregonly
		} else {
			set reannot -reannot
		}
		set names {}
		foreach {cgdir name} $data {
			lappend names $name
		}
		if {[file exists compar/${project}_compar.tsv]} {
			set done [mcompar_samples compar/${project}_compar.tsv]
			set names [list_lremove $names $done]
		} else {
			set done {}
		}
		puts "Already done: [join $done {, }]"
		if {[llength $names]} {
			file mkdir compar
			catch {file delete [gzfile compar/annot${project}_compar.tsv]}
			set cjob [submit -direct $direct -host lungo -deps $jobs cg multicompar $reannot compar/${project}_compar.tsv {*}$names]
			if {[isint $cjob]} {lappend alljobs $cjob}
		}
		if {![file exists [gzfile compar/annot${project}_compar.tsv]]} {
			set ajob [submit -direct $direct -host lungo -deps [get cjob {}] cg annotate compar/${project}_compar.tsv compar/annot${project}_compar.tsv $refseqdir/$build $refseqdir/$build/annovar]
	 		if {[isint $ajob]} {lappend alljobs $ajob}
		}
		# multireg
		# --------
		set resultfile compar/reg-${project}.tsv
		puts "Checking [file normalize $resultfile]"
		set done {}
		if {[file exists $resultfile]} {
			set list [cg select -h $resultfile]
			set poss [list_find -glob $list sreg-*]
			set done [list_sub $list $poss]
		}
		set files {}
		set names {}
		foreach {cgdir name} $data {
			if {![inlist $done sreg-$name]} {
				lappend files $name/sreg-$name.tsv
				lappend names $name
			}
		}
		if {[llength $done]} {
			puts "Multireg already done: $done"
		}
		if {[llength $files]} {
			set crjob [submit -direct $direct -deps $jobs cg multireg $resultfile {*}$files]
	 		if {[isint $crjob]} {lappend alljobs $crjob}
		}
		# cgsv
		# ----
		set done {}
		set resultfile cgsv-${project}.tsv
		puts "Checking [file normalize compar/$resultfile]"
		set names {}
		if {[file exists compar/$resultfile]} {
			set list [cg select -h compar/$resultfile]
			set poss [list_find -glob $list start1-*]
			set done [list_sub $list $poss]
			set done [list_regsub -all {^start1-} $done {}]
		}
		set files {}
		set names {}
		foreach {cgdir name} $data {
			if {![inlist $done $name]} {
				lappend files $name/cgsv-$name.tsv
				lappend names $name
			}
		}
		if {[llength $done]} {
			puts "Multicgsv already present: $done"
		}
		if {[llength $files]} {
			set cgsvjob [submit -direct $direct -deps $jobs cg svmulticompar compar/$resultfile {*}$files]
	 		if {[isint $cgsvjob]} {lappend alljobs $cgsvjob}
		}
		if {![file exists [gzfile compar/annot$resultfile]] || [llength $files]} {
			set cmd [list cg annotate compar/$resultfile compar/annot$resultfile {*}[glob -nocomplain $refseqdir/$build/reg_*.tsv $refseqdir/$build/gene_*.tsv]]
			set ajob [submit -direct $direct -host lungo -deps [get cgsvjob {}] {*}$cmd]
	 		if {[isint $ajob]} {lappend alljobs $ajob}
		}
	}

	# sv (our algorithm)
	# ==================
	if {[inlist $actions sv]} {
		cd $resultdir
		set svjobs {}
		foreach {cgdir name} $data {
			set dir [file dir $cgdir]
			set host [lindex [file split $cgdir] 2]
			set job [submit -direct $direct -host $host cg process_sv $cgdir $name $refseqdir]
			if {[isint $job]} {lappend svjobs $job}
		}
		lappend alljobs {*}$svjobs
		# qalter -u peter -p 1
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
				set job [submit -direct $direct -deps $alljobs bgzip $file]
			}
			# compress sv data
			set dir [file dir $cgdir]
			set files [glob -nocomplain $name/sv/*-paired.tsv]
			foreach file $files {
				puts "Compressing $file"
				set job [submit -direct $direct -deps $alljobs bgzip $file]
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
			if {![file exists $project]} {
				exec ln -s $projectdir .
			}
		}
	}
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
		set temp [eval {exec submit.tcl -host lungo -deps [join $svjobs ,] cg svmulticompar compar/${project}_svcompar.tsv} $names]
	} else {
		set temp [eval {exec submit.tcl -host lungo cg svmulticompar compar/${project}_compar.tsv} $names]
	}

}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base $scriptname
	cg_project {*}$argv
}

