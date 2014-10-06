#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc multidb_merge_job {varsfile files {split 1}} {
	global multi_merge_num
	upvar job_logdir job_logdir
	set newfiles {}
	if {[llength $files] < 2} {
		mklink [lindex $files 0] $varsfile
		return
	} 
	foreach {file1 file2} $files {
		if {$file2 eq ""} continue
		incr multi_merge_num
		job multi_merge-$multi_merge_num -deps [list $file1 $file2] -vars split -targets $varsfile.$multi_merge_num -code {
			set f [gzopen $dep1]
			set header [tsv_open $f]
			set poss1 [tsv_basicfields $header 6 1]
			lappend poss1 [lsearch $header id]
			close $f
			set f [gzopen $dep2]
			set header [tsv_open $f]
			set poss2 [tsv_basicfields $header 6 1]
			lappend poss2 [lsearch $header id]
			close $f
			exec multi_merge $dep1 {*}$poss1 $dep2 {*}$poss2 $split > $target.temp
			file rename -force $target.temp $target
		}
		lappend newfiles $varsfile.$multi_merge_num
	}
	if {$file2 eq ""} {
		lappend newfiles $file1
	}
	if {[llength $newfiles] > 1} {
		multidb_merge_job $varsfile $newfiles $split
	} else {
		incr multi_merge_num
		job multi_merge-$multi_merge_num -deps $newfiles -vars split -targets $varsfile -code {
			file rename -force $dep1 $target
		}
	}
	
}

proc multidb_getfileinfo {dirs targetsfield aVar datafilesVar genofieldsVar ngssequpdate_file} {
	upvar $aVar a
	upvar $datafilesVar datafiles
	upvar $genofieldsVar genofields
	set datafiles {}
	foreach dir $dirs {
		set dir [file_absolute $dir]
		if {[file isdir $dir]} {
			# a directory is given, look for the variant file(s)
			set name [file tail $dir]
			if {[file exists $dir/fannotvar-$name.tsv]} {
				set file [gzfile $dir/fannotvar-$name.tsv]
				lappend datafiles $file
				if {[info exists a(file,$name)]} {
					error "sample \"$name\" in \"$file\" was already in \"$a(file,$name)\""
				}
				set a(file,$name) $file
			} else {
				foreach file [gzfiles $dir/var-*-$name.tsv] {
					set base [file root [file tail [gzroot $file]]]
					set file [gzfile $file]
					set name [sourcename $base]
					lappend datafiles $file
					if {[info exists a(file,$name)]} {
						error "sample \"$name\" in \"$file\" was already in \"$a(file,$name)\""
					}
					set a(file,$name) $file
				}
			}
		} elseif {[llength [cg select -n $dir]]} {
			set samples [cg select -n $dir]
			set file $dir
			lappend datafiles $file
			foreach name $samples {
				if {[info exists a(file,$name)]} {
					error "sample \"$name\" in \"$file\" was already in \"$a(file,$name)\""
				}
				set a(file,$name) $file
			}
		} else {
			set base [file root [file tail [gzroot $dir]]]
			set name [sourcename $base]
			set file [gzfile $dir]
			lappend datafiles $file
			if {[info exists a(file,$name)]} {
				error "sample \"$name\" in \"$file\" was already in \"$a(file,$name)\""
			}
			set a(file,$name) $file
		}
	}
	# check samples, find all fields
	set fs [open $ngssequpdate_file w]
	set varcols {id chomosome begin end type ref alt id}
	set genofields {var ngsseq sequenced zyg alleleSeq1 alleleSeq2 genotypes quality coverage}
	set sampleid 1
	foreach file $datafiles {
		set f [gzopen $file]
		set header [tsv_open $f]
		close $f
		set file [file_absolute $file]
		set a(header,$file) $header
		set dir [file dir $file]
		set filebase [file root [file tail [gzroot $file]]]
		set basicposs [tsv_basicfields $header]
		set samples [samples $header]
		if {[llength $samples] > 0} {
			set a(samples,$file) $samples
			foreach sample $samples {
				set poss [list_find -glob $header *-$sample]
				set len [string length $sample]
				incr len
				set fields {}
				foreach field [list_sub $header $poss] {
					lappend fields [string range $field 0 end-$len]
				}
				lappend genofields {*}$fields
				set a(id,$sample) $sampleid
				set a(file,$sample) $file
				set a(fields,$sample) $fields
				set a(poss,$sample) $poss
				# in a multicompar file, sequenced columns should already be present
				# put a value not -1 or -2, so no new one will be added
				set a(seqpos,$sample) 0
				puts $fs [join [list $sampleid $sample] \t]
				incr sampleid
			}
		} else {
			set sample [sourcename $filebase]
			set a(samples,$file) [list $sample]
			set seqpos [lsearch $header sequenced]
			set fields [list_sub $header -exclude $basicposs]
			if {[string match fannotvar-* [file tail $file]]} {
				set fields [list_lremove $fields $mergefields]
			}
			set fields [list_remove $fields $targetsfield]
			set poss [list_cor $header $fields]
			lappend genofields {*}$fields
			set a(id,$sample) $sampleid
			set a(file,$sample) $file
			set a(fields,$sample) $fields
			set a(poss,$sample) $poss
			set a(seqpos,$sample) $seqpos
			puts $fs [join [list $sampleid $sample] \t]
			incr sampleid
		}
	}
	close $fs
	set genofields [list_remdup $genofields]
	return $datafiles
}

proc cg_multidb {args} {

	set args [job_init -silent {*}$args]
	set reannot 0
	set regonly 0
	set split 0
	set listfields {}
	set targetsfile {}
	set targetsfield {}
	set pos 0
	set varidstart 1
	while 1 {
		set key [lindex $args $pos]
		switch -- $key {
			-varidstart {
				incr pos
				set varidstart [lindex $args $pos]
			}
			-reannot {
				putslog "Also reannot"
				set reannot 1
			}
			-reannotregonly {
				putslog "Also reannot"
				set reannot 1
				set regonly 1
			}
			-split {
				incr pos
				set split [true [lindex $args $pos]]
			}
			-listfields {
				incr pos
				set listfields [lindex $args $pos]
			}
			-targetsfile {
				incr pos
				set targetsfile [lindex $args $pos]
				set targetsfield [lindex [split [file root [file tail $targetsfile]] -] end]
			}
			-- break
			default {
				if {[string index $key 0] eq "-"} {error "unknown option \"$key\""}
				break
			}
		}
		incr pos
	}
	set args [lrange $args $pos end]
	if {([llength $args] < 1)} {
		puts "Wrong number of arguments"
		errorformat multidb
		exit 1
	}
	foreach {compar_dir} $args break
	if {[file exists $compar_dir] && ![file isdir $compar_dir]} {file delete $compar_dir}
	file mkdir $compar_dir
	set compar_file $compar_dir/vars.tsv
	set workdir $compar_dir/work
	# should take into account existing instead of deleting and starting all over -> not now
	if {[file exists $workdir]} {file delete -force $workdir}
	file mkdir $workdir
	set dirs [lrange $args 1 end]
	#
	# get all datafiles to be processed from $dirs
	# store file and sample info in array a
	# genofields contains all fields used in any of the data files
	# The file $compar_dir/ngsseq.tsv.update contain the ids of all ngsseq and additional data
	unset -nocomplain a
	multidb_getfileinfo $dirs $targetsfield a datafiles genofields $compar_dir/ngsseq.tsv.update
	#
	# merge variants
	# todo: check for concurrency
	job_logdir $workdir/log_jobs
	set multi_merge_num 0
	set files $datafiles
	# compar_file and targetsfile will only add variants (not data), 
	# so they are used as a variant file in merge, but not later
	lappend files $compar_file
	if {$targetsfile ne ""} {lappend files $targetsfile}
	multidb_merge_job $workdir/vars.tsv $files $split
	#
	# make todo file
	set todofile $workdir/multitodo.txt
	set ftodo [open $todofile w]
	# line with:
	#     filename chrpos beginpos endpos typepos refpos altpos idpos
	#     maxcoltoread todoseqpos sizeofnextline
	puts $ftodo [join [list $workdir/vars.tsv 0 1 2 3 4 5 6   6 0 0] \t]
	# line with: data positions that must go in the multicompar
	puts $ftodo ""
	set reannotheader {chromosome begin end type ref alt}
	set addfields [list_remove $genofields var ngsseq]
	foreach file $datafiles {
		set keepposs {}
		foreach sample $a(samples,$file) {
			if {$a(seqpos,$sample) != -1} {
				set cor [list_cor $a(fields,$sample) $addfields]
			} else {
				set cor [list_cor $a(fields,$sample) [lrange $addfields 1 end]]
			}
			set poss [list_sub $a(poss,$sample) $cor]
			set poss [list_change $poss {{} -1}]
			lappend keepposs -2 $a(id,$sample) {*}$poss
		}
		set max [lmath_max [list_concat $keepposs $basicposs]]
		# line with: filename chrpos beginpos endpos typepos refpos altpos maxcoltoread todoseqpos sizeofnextline
		puts $ftodo [join [list $file {*}$basicposs -1 $max $a(seqpos,$sample) [llength $keepposs]] \t]
		# line with: data positions that must go in the multicompar
		puts $ftodo [join $keepposs \t]
	}
	set len [expr {[llength $datafiles]+1}]
	close $ftodo
	#
	# write headers of targets
	file_write $compar_dir/vars.tsv.new [join {chromosome begin end type ref alt id} \t]\n
	file_write $compar_dir/vars.tsv.update [join {chromosome begin end type ref alt id} \t]\n
	file_write $compar_dir/geno.tsv.update [join $genofields \t]\n
	# run
	set deps [list $todofile $workdir/vars.tsv {*}$datafiles]
	set targets [list $compar_dir/vars.tsv.new $compar_dir/vars.tsv.update $compar_dir/geno.tsv.update]
	job multi_join -deps $deps -vars {split len varidstart} -targets $targets -code {
		exec multidb_join [lindex $deps 0] $len $split {*}$targets $varidstart
		file rename -force $target.temp $target
	}
	if {$reannot} {
		job multicompar_reannot -deps {$compar_file} -vars {regonly} -targets {$compar_file $compar_file.reannot} -code {
			putslog "Reannotating $dep"
			multicompar_reannot $dep 0 $regonly
			file_write $dep.reannot [timestamp]
		}
	}
	job_wait
}
