#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc multidb_analysisinfo {sampleid file sample iscompar} {
	foreach {experiment ngsseq mapper varcall gentli_ngsseq reference individual samplenum} {? ? ? ? ? ? ? ?} break
	set dirlist [file split [file dir $file]]
	set filebase [file root [file tail $file]]
	if {$iscompar} {
		set experiment [lindex [split $filebase -] end]
	} else {
		set experiment [lindex $dirlist end-1]
		if {$experiment eq "samples"} {
			set experiment [lindex $dirlist end-2]
		}
	}
	set split [split $sample -]
	if {[llength $split] == 3} {
		foreach {varcall mapper sample} $split break
	}
	regsub {_[0-9]+$} $sample {} individual
	list $sampleid $experiment $sample $mapper $varcall $gentli_ngsseq $reference $individual $samplenum
}

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
			set temptarget [filetemp $target]
			exec multi_merge $split $dep1 $dep2 > $temptarget
			file rename -force $temptarget $target
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

proc multidb_getfileinfo {dirs aVar datafilesVar genofieldsVar compar_dir} {
	upvar $aVar a
	upvar $datafilesVar datafiles
	upvar $genofieldsVar genofields
	set mergefields {xRef geneId mrnaAcc proteinAcc symbol orientation component componentIndex hasCodingRegion impact nucleotidePos proteinPos annotationRefSequence sampleSequence genomeRefSequence pfam}
	set datafiles {}
	if {![info exists genofields]} {
		set genofields {var analysis sequenced zyg alleleSeq1 alleleSeq2 quality coverage}
	}
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
		} elseif {[file exists $dir] && [llength [cg select -n $dir]]} {
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
	set fs [open $compar_dir/analysis.tsv.insert w]
	puts $fs [join {id experiment ngsseq mapper varcall gentli_ngsseq reference individual sample} \t]
	set varcols {id chomosome begin end type ref alt id}
	# get next sampleid
	if {[file exists $compar_dir/analysis.tsv.maxid]} {
		set sampleid [file_read $compar_dir/analysis.tsv.maxid]
	} elseif {[file exists $compar_dir/analysis.tsv]} {
		set sampleid [lindex [cg select -g all -gc {max(id)} $compar_dir/analysis.tsv] end]
		if {$sampleid eq "all"} {set sampleid 0}
	} else {
		set sampleid 0
	}
	incr sampleid
	set insertcount 0
	foreach file $datafiles {
		set f [gzopen $file]
		set header [tsv_open $f]
		gzclose $f
		set file [file_absolute $file]
		set a(header,$file) $header
		set dir [file dir $file]
		set filebase [file root [file tail [gzroot $file]]]
		set basicposs [tsv_basicfields $header]
		set a(basicposs,$file) $basicposs
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
				puts $fs [join [multidb_analysisinfo $sampleid $file $sample 1] \t]
				incr insertcount
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
			set poss [list_cor $header $fields]
			lappend genofields {*}$fields
			set a(id,$sample) $sampleid
			set a(file,$sample) $file
			set a(fields,$sample) $fields
			set a(poss,$sample) $poss
			set a(seqpos,$sample) $seqpos
			puts $fs [join [multidb_analysisinfo $sampleid $file $sample 0] \t]
			incr insertcount
			incr sampleid
		}
	}
	gzclose $fs
	file_write $compar_dir/analysis.tsv.insert.maxid $sampleid
	file_write $compar_dir/analysis.tsv.insert.count $insertcount
	set genofields [list_remdup $genofields]
	return $datafiles
}

proc multidb_job {args} {
	set reannot 0
	set targetsfile {}
	set targetsfield {}
	set dbdir /complgen/refseq/hg19
	cg_options multidb args {
		-monetdb {
			set monetdb $value
		}
		-dbdir {
			set dbdir $value
		}
		-split {
			set split [true $value]
		}
		-targetsfile {
			set targetsfile $value
			set targetsfield [lindex [split [file root [file tail $targetsfile]] -] end]
		}
	} compar_dir 1
	set dirs [lrange $args 1 end]

	if {[file exists $compar_dir/vars.tsv.insert]
		|| [file exists $compar_dir/analysis.tsv.insert]
		|| [file exists $compar_dir/geno.tsv.insert]
		|| [file exists $compar_dir/analysis.tsv.update]
		|| [file exists $compar_dir/geno.tsv.update]
	} {
		error "unimported updates/inserts present in $compar_dir"
	}
	#
	set compar_dir [file_absolute $compar_dir]
	if {[file exists $compar_dir] && ![file isdir $compar_dir]} {file delete $compar_dir}
	if {![file exists $compar_dir]} {
		file mkdir $compar_dir
	}
	projectinfo $compar_dir {monetdb {}} {split 1}
	if {$monetdb ne ""} {
		set genofields [multidb_monet_open $compar_dir $monetdb]
	} else {
		set genofields [multidb_dir_open $compar_dir $monetdb]
	}
	set compar_file $compar_dir/vars.tsv
	set workdir $compar_dir/work
	# should take into account existing instead of deleting and starting all over -> not now
	if {[file exists $workdir]} {file delete -force $workdir}
	file mkdir $workdir
	#
	# get all datafiles to be processed from $dirs
	# store file and sample info in array a
	# genofields contains all fields used in the database (from previous call) or in any of the data files
	# The file $compar_dir/analysis.tsv.insert contain the ids of all analysis and additional data
	unset -nocomplain a
	multidb_getfileinfo $dirs a datafiles genofields $compar_dir
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
	set addfields [list_remove $genofields var analysis]
	foreach file $datafiles {
		set basicposs $a(basicposs,$file)
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
	# run
	set deps [list $todofile $workdir/vars.tsv {*}$datafiles]
	set targets [list $compar_dir/vars.tsv.new $compar_dir/vars.tsv.insert $compar_dir/geno.tsv.insert]
	job multidb_join -deps $deps -vars {split len genofields compar_dir} -targets $targets -code {
		set temptargets {}
		if {[file exists $dep.maxid]} {
			set varidstart [file_read $dep.maxid]
		} else {
			set varidstart [lindex [cg select -g all -gc {max(id)} $dep] end]
			if {$varidstart eq "all"} {set varidstart 1}
		}
		incr varidstart
		# write headers of targets
		file_write $compar_dir/vars.tsv.new.temp [join {chromosome begin end type ref alt id} \t]\n
		file_write $compar_dir/vars.tsv.insert.temp [join {chromosome begin end type ref alt id} \t]\n
		file_write $compar_dir/geno.tsv.insert.temp [join $genofields \t]\n
		foreach t $targets {lappend temptargets $t.temp}
		exec multidb_join [lindex $deps 0] $len $split {*}$temptargets $varidstart
		foreach t $targets {
			file rename -force $t.temp $t
		}
		file rename $compar_dir/vars.tsv.new.temp.maxid $compar_dir/vars.tsv.new.maxid
		file rename $compar_dir/vars.tsv.new.temp.count $compar_dir/vars.tsv.new.count
		file rename $compar_dir/vars.tsv.insert.temp.count $compar_dir/vars.tsv.insert.count
		file rename $compar_dir/geno.tsv.insert.temp.count $compar_dir/geno.tsv.insert.count
	}
	job multidb_annot -deps $compar_dir/vars.tsv.insert -vars {compar_dir dbdir} -targets $compar_dir/annot.tsv.insert -code {
		cg annotate -multidb 1 $compar_dir/vars.tsv.insert $compar_dir/annot.tsv.insert.temp \
			$dbdir {*}[glob -nocomplain $dbdir/extra/var_*_evs.tsv $dbdir/extra/var_*_dbnsfp.tsv]
		file rename $compar_dir/annot.tsv.insert.temp $compar_dir/annot.tsv.insert
		file copy -force $compar_dir/vars.tsv.insert.count $compar_dir/annot.tsv.insert.count
	}
}

proc cg_multidb {args} {
	set args [job_init {*}$args]
	multidb_job {*}$args
	job_wait
}

proc cg_multidb_import {args} {
	cg_options multidb_import args {
		-monetdb {
			set monetdb $value
		}
	} compar_dir
	projectinfo $compar_dir {monetdb {}} {split 1}
	if {$monetdb ne ""} {
		multidb_monet_import $compar_dir $monetdb
	} else {
		multidb_dir_import $compar_dir
	}
}
