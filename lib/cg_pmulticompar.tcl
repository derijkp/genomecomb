#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc multi_merge_job {varsfile files args} {
	set workdir [workdir $varsfile]/paste
	file mkdir $workdir
	upvar job_logdir job_logdir
	set optional 0
	set split 1
	set force 1
	set limitreg {}
	foreach {k v} $args {
		switch $k {
			-split {set split $v}
			-force {set force $v}
			-limitreg {set limitreg $v}
			default {error "Unkown option $k"}
		}
	}
	set varsfile [file_absolute $varsfile]
	set maxfiles [maxopenfiles]
	set len [llength $files]
	if {$len <= $maxfiles} {
		# merge in one go
		set target $varsfile
		job [job_relfile2name multi_merge- $varsfile] -force $force -deps $files -targets {
			$target
		} -vars {
			split limitreg workdir
		} -code {
			# puts [list ../bin/multi_merge $split {*}$deps]
			if {$limitreg ne ""} {
				set templist {}
				foreach file $deps {
					set tempfile [tempfile]
					lappend templist $tempfile
					if {[catch {
						exec cg regselect $file $limitreg > $tempfile
					} msg]} {
						error "error selection region $limitreg from file $file: $msg"
					}
				}
				set deps $templist
			}
			exec multi_merge $split {*}$deps > $workdir/target.merge.temp 2>@ stderr
			file rename -force -- $workdir/target.merge.temp $target
		}
		return
	}
	catch {file delete {*}[glob -nocomplain $workdir/paste.temp*]}
	set todo $files
	set delete 0
	set num 1
	while 1 {
		if {$len <= $maxfiles} {
			# final merge (we have less than maxfiles files to merge left)
			set target $varsfile
			job [job_relfile2name multi_merge- $varsfile] -optional $optional -force $force -deps $todo -targets {
				$target
			} -vars {
				split delete workdir limitreg
			} -code {
				# puts [list ../bin/multi_merge $split {*}$deps]
				exec multi_merge $split {*}$deps > $target.merge.temp 2>@ stderr
				file rename -force -- $target.merge.temp $target
				if {$delete} {file delete {*}$deps}
			}
			break
		}
		set pos 0
		set newtodo {}
		while {$pos < $len} {
			# go over sets (size maxfiles) and merge to temporary files that will be merged in the final merge
			set target $workdir/paste.temp$num
			incr num
			lappend newtodo $target
			set deps [lrange $todo $pos [expr {$pos+$maxfiles-1}]]
			incr pos $maxfiles
			job [job_relfile2name multi_merge- $target] -optional $optional -force $force -deps $deps -targets {
				$target
			} -vars {
				split delete limitreg
			} -code {
				set deps [gzfile_multi $deps]
				if {[llength $deps] > 1} {
					# puts [list ../bin/multi_merge $split {*}$deps]
					if {!$delete && $limitreg ne ""} {
						set templist {}
						foreach file $deps {
							set tempfile [tempfile]
							lappend templist $tempfile
							if {[catch {
								exec cg regselect $file $limitreg > $tempfile
							} msg]} {
								error "error selection region $limitreg from file $file: $msg"
							}
						}
						set deps $templist
					}
					exec multi_merge $split {*}$deps > $target.merge.temp 2>@ stderr
					if {$delete} {file delete {*}$deps}
				} elseif {!$delete} {
					# one file left, does not have to be merged
					if {$limitreg ne ""} {
						exec cg regselect $dep $limitreg > $target.merge.temp
					} else {
						mklink $dep $target.merge.temp
					}
				} else {
					file rename -force -- $dep $target.merge.temp
				}
				file rename -force -- $target.merge.temp $target
			}
		}
		set delete 1
		set todo $newtodo
		set len [llength $todo]
	}
}

proc pmulticompar_findsample {basedir sample args} {
	if {![llength $args]} {set args [list {}]}
	set sampledir [lindex [split $sample -] end]
	foreach startdir [list $basedir/samples $basedir [file dir $basedir]/samples [file dir $basedir]] {
		foreach pattern $args {
			foreach usesample [list $sample $sampledir {}] {
				set test [file join $startdir $usesample $pattern]
				if {[jobfileexists $test]} {
					return $test
				}
			}
		}
	}
	return {}
}

proc multicompar_tcl_addvars_coverageRefScore {data restVar pos chromosome begin} {
	upvar $restVar rest
	foreach {covpos refscorepos cr_file f poss rchr rpos rcoverage rrefscore} $data break
	if {$rchr ne $chromosome} {
		if {$f ne ""} {gzclose $f}
		set tail [split [file tail $cr_file] -]
		set newfile [file dir $cr_file]/[join [lrange $tail 0 end-2] -]-$chromosome-[lindex $tail end]
		if {[file exists $newfile]} {
			set f [gzopen $newfile]
			set rchr $chromosome
			set h [tsv_open $f]
			set poss [list_cor $h {offset uniqueSequenceCoverage refScore}]
			set l [split [gets $f] \t]
			foreach {cr_prevpos cr_prevcov cr_prevrefscore} $l break
		} else {
			if {![inlist {Y chrY} $chromosome]} {
				putslog "coverage(RefScore) file not found ($newfile))"
			}
			set f ""
			set rchr ""
		}
	}
	if {$f ne ""} {
		set match 0
		while 1 {
			if {$rpos > $begin} break
			if {$rpos == $begin} {
				set match 1 ; break
			}
			if {[gets $f line] == -1} break
			set line [split $line \t]
			foreach {rpos rcoverage rrefscore} [list_sub $line $poss] break
		}
		if {$match} {
			if {$covpos != -1} {lset rest $covpos $rcoverage}
			if {$refscorepos != -1} {lset rest $refscorepos $rrefscore}
		} else {
			lset rest $covpos 0
		}
	}
	list $covpos $refscorepos $cr_file $f $poss $rchr $rpos $rcoverage $rrefscore
}

proc multicompar_tcl_addvars {sample target split allvarsfile samplevarsfile sregfile varallfile bcolannot oldbcolannot regfiles coverageRefScorefiles keepfields} {
#puts "Using Tcl version"
# putsvars sample target split allvarsfile samplevarsfile sregfile varallfile bcolannot oldbcolannot regfiles coverageRefScorefiles keepfields
	catch {close $allvars} ; set allvars [gzopen $allvarsfile]
	set allheader [tsv_open $allvars]
	if {$allheader ne "chromosome begin end type ref alt"} {error "internal error: index vars.tsv in wrong format"}
	set allposs [tsv_basicfields $allheader 6 0]
	#
	catch {close $orivars} ; set orivars [gzopen $samplevarsfile]
	set oriheader [tsv_open $orivars]
	set keepposs [list_cor $oriheader $keepfields]
	set cposs [list_cor $oriheader {sequenced zyg alleleSeq1 alleleSeq2}]
	foreach {oriseqpos orizygpos oria1pos oria2pos} $cposs break
	set oriposs [tsv_basicfields $oriheader 6 0]
	lappend oriposs {*}$cposs
	# blanks will be used as a start for rest, see which sources we have to fil in the blanks
	set blanks [list_fill [llength $keepposs] ?]
	set todo {}
	unset -nocomplain todoa
	foreach {pos bcolfile} $bcolannot {
		if {[info exists todoa($pos)]} continue
		lappend todo $pos bcol
		set todoa($pos) [bcol_open $bcolfile]
	}
	foreach {pos bcolfile} $oldbcolannot {
		if {[info exists todoa($pos)]} continue
		lappend todo $pos oldbcol
		set todoa($pos) [bcol_open $bcolfile]
	}
	foreach {pos regfile scorepos} $regfiles {
		if {[info exists todoa($pos)]} continue
		set f [gzopen $regfile]
		set header [tsv_open $f]
		gzclose $f
		set poss [tsv_basicfields $header 3]
		lappend poss $scorepos
		lappend todo $pos regfile
		set line [gets $f]
		foreach {rchr rbegin rend rscore} [list_sub $line $poss] break
		set todoa($pos) [list $f $poss $scorepos $rchr $rbegin $rend $rscore]
	}
	if {$coverageRefScorefiles ne ""} {
		set covpos [lsearch $keepfields coverage]
		if {[info exists todoa($covpos)]} {set covpos -1}
		set refscorepos [lsearch $keepfields refscore]
		if {[info exists todoa($refscorepos)]} {set refscorepos -1}
		if {$covpos != -1} {
			lappend todo $covpos coverageRefScore
			set todoa($covpos) [list $covpos $refscorepos $coverageRefScorefiles]
		} elseif {$refscorepos != -1} {
			lappend todo $refscorepos coverageRefScore
			set todoa($refscorepos) [list $covpos $refscorepos $coverageRefScorefiles]
		}
	}
	# first oriline
	set oriline [split [gets $orivars] \t]
	set oriloc [list_sub $oriline $oriposs]
	foreach {orichromosome oribegin oriend oritype oriref orialt oriseq orizyg oria1 oria2} [list_sub $oriline $oriposs] break
	set orichromosome [chr_clip $orichromosome]
	if {$oritype eq "del"} {
		set oridelchromosome $orichromosome; set oridelbegin $oribegin ; set oridelend $oriend ; set oridelzyg $orizyg
	} else {
		set oridelchromosome "" ; set oridelbegin -1 ; set oridelend -1 ; set oridelzyg m;
	}
	#
	if {[file exists $sregfile]} {
		catch {close $sreg} ; set sreg [gzopen $sregfile]
		set sregheader [tsv_open $sreg]
		set sregposs [tsv_basicfields $sregheader 3 0]
		set sregline [split [gets $sreg] \t]
		foreach {sregchromosome sregbegin sregend} [list_sub $sregline $sregposs] break
		set sregchromosome [chr_clip $sregchromosome]
	} else {
		set sreg ""
	}
	#
	catch {close $o} ; set o [open "[compresspipe $target 1 1] > $target" w]
	set newheader [list sequenced-$sample]
	if {$orizygpos != -1} {lappend newheader zyg-$sample}
	if {$oria1pos != -1} {lappend newheader alleleSeq1-$sample}
	if {$oria2pos != -1} {lappend newheader alleleSeq2-$sample}
	foreach field $keepfields {
		lappend newheader $field-$sample
	}
	puts $o [join $newheader \t]
	#
	set prevcomp -2
	while 1 {
		if {[gets $allvars line] == -1} break
		set line [split $line \t]
		foreach {chromosome begin end type ref alt} $line break
		set chromosome [chr_clip $chromosome]
# putsvars chromosome begin end type orichromosome oribegin oriend oritype
		if {$orichromosome eq $chromosome && $oribegin == $begin && $oriend == $end && $oritype == $type} {
			if {$orialt eq $alt} {
				set comp 0
			} else {
				set comp -1
			}
		} else {
			set comp -2
		}
		if {$comp == 0 || (!$split && $comp == -1)} {
			# match, output
			if {$oriseq eq ""} {set oriseq v}
			set result $oriseq
			if {$orizygpos != -1} {lappend result $orizyg}
			if {$oria1pos != -1} {lappend result $oria1}
			if {$oria2pos != -1} {lappend result $oria2}
			lappend result {*}[list_sub $oriline $keepposs]
			puts $o [join $result \t]
			# next oriline
			set prevline $oriline ; set prevcomp $comp ; set prevchromosome $chromosome ; set prevbegin $oribegin ; set prevend $oriend ; set prevtype $oritype ; set preva1 $oria1 ; set preva2 $oria2
			if {[gets $orivars oriline] != -1} {
				set oriline [split $oriline \t]
				foreach {orichromosome oribegin oriend oritype oriref orialt oriseq orizyg oria1 oria2} [list_sub $oriline $oriposs] break
				set orichromosome [chr_clip $orichromosome]
				if {$oritype eq "del"} {
					set oridelchromosome $orichromosome; set oridelbegin $oribegin ; set oridelend $oriend ; set oridelzyg $orizyg
				} elseif {$orichromosome ne $oridelchromosome} {
					set oridelchromosome "" ; set oridelbegin -1 ; set oridelend -1 ; set oridelzyg m;
				}
			} else {
				set orichromosome {}
				set oridelchromosome "" ; set oridelbegin -1 ; set oridelend -1 ; set oridelzyg m;
			}
			continue
		}
		set rest $blanks
		if {$split && $comp == -1} {
			set out_seq r
			set out_a1 $oria1
			set out_a2 $oria2
			if {$out_a1 ne $ref || $out_a2 ne $ref} {
				set out_zyg o
			} else {
				set out_zyg r
			}
			# do not take values from ori (some values are likely allele specific)
			# set rest [list_sub $oriline $keepposs]
		} elseif {$split && $prevcomp == 0 && $prevchromosome eq $chromosome && $prevbegin == $begin && $prevend == $end && $prevtype == $type} {
			# if split and match to previous (except allele), take data of previous
			# putsvars split prevcomp chromosome begin end type prevchromosome prevbegin prevend prevtype
			set out_seq r
			set out_a1 $preva1
			set out_a2 $preva2
			set out_zyg ?
			if {$out_a1 ne $ref || $out_a2 ne $ref} {
				set out_zyg o
			} else {
				set out_zyg r
			}
			# do not take values from ori (some values are likely allele specific)
			# set rest [list_sub $prevline $keepposs]
		} elseif {$sreg eq ""} {
			set out_seq ?
			set out_zyg ?
			set out_a1 ?
			set out_a2 ?
		} else {
			while 1 {
				set chrcomp [loc_compare $sregchromosome $chromosome]
				if {$chrcomp > 0} break
				if {$chrcomp == 0 && $sregend > $begin} break
				if {[gets $sreg sregline] == -1} {set chrcomp 2; break}
				set sregline [split $sregline \t]
				foreach {sregchromosome sregbegin sregend} [list_sub $sregline $sregposs] break
				set sregchromosome [chr_clip $sregchromosome]
			}
			if {$chrcomp == 0 && $begin >= $sregbegin && $end <= $sregend} {
				# putsvars chromosome begin end oridelchromosome oridelend oridelbegin
				if {$oridelchromosome eq $chromosome && $begin < $oridelend && $end > $oridelbegin} {
					set out_seq r
					set out_zyg o
					if {$oridelzyg == "m"} {
						set out_a1 @ ; set out_a2 @
					} elseif {$oridelzyg in "t c"} {
						set out_a1 $ref ; set out_a2 @
					} else {
						set out_a1 $ref ; set out_a2 $ref
					}
				} else {
					set out_seq r
					set out_zyg r
					set out_a1 $ref
					set out_a2 $ref
				}
			} else {
				set out_seq u
				set out_zyg u
				set out_a1 ?
				set out_a2 ?
			}
		}
		foreach {pos type} $todo {
			switch $type {
				bcol {
					set bcol $todoa($pos)
					set value [bcol_get $bcol $begin $begin $chromosome]
					lset rest $pos $value
				}
				oldbcol {
					set bcol $todoa($pos)
					if {[lindex [dict get $bcol tablechr] 0] ne $chromosome} {
						set oldfile [dict get $bcol file]
						set tail [split [file tail $oldfile] -]
						set newfile [file dir $oldfile]/[join [lrange $tail 0 end-2] -]-$chromosome-[lindex $tail end]
						bcol_close $bcol
						if {[file exists $newfile]} {
							set bcol [bcol_open $newfile]
						} elseif {[inlist {Y chrY} $chr]} {
							set bcol [dict create tablechr $chromosome notpresent 1]
						} else {
							putslog "bcol file not found ($newfile)"
							set bcol [dict create tablechr $chromosome notpresent 1]
						}
						set todoa($pos) $bcol
					}
					if {![dict exists $bcol notpresent]} {
						set value [bcol_get $bcol $begin $begin $chromosome]
						lset rest $pos $value
					}
				}
				regfile {
					foreach {f poss scorepos rchr rbegin rend rscore} $todoa($pos) break
					# putsvars chromosome begin end rchr rbegin rend rscore
					set match 0
					while 1 {
						set chrcomp [loc_compare $rchr $chromosome]
						if {$chrcomp > 0} break
						if {$chrcomp == 0} {
							if {$rbegin >= $end} break
							if {$rend > $begin} {
								set match 1 ; break
							}
						}
						if {[gets $f line] == -1} break
						set line [split $line \t]
						foreach {rchr rbegin rend rscore} [list_sub $line $poss] break
					}
					if {$scorepos == -1} {set rscore 1}
					if {$match} {
						lset rest $pos $rscore
					} else {
						lset rest $pos {}
					}
					set todoa($pos) [list $f $poss $scorepos $rchr $rbegin $rend $rscore]
				}
				coverageRefScore {
					set todoa($pos) [multicompar_tcl_addvars_coverageRefScore $todoa($pos) rest $pos $chromosome $begin]
				}
			}
		}
		set result $out_seq
		if {$orizygpos != -1} {lappend result $out_zyg}
		if {$oria1pos != -1} {lappend result $out_a1}
		if {$oria2pos != -1} {lappend result $out_a2}
		lappend result {*}$rest
		puts $o [join $result \t]
		set prevcomp $comp
	}
	foreach {pos type} $todo {
		switch $type {
			bcol {
				bcol_close $todoa($pos)
			}
			oldbcol {
				bcol_close $todoa($pos)
			}
			regfile {
				foreach {f poss scorepos rchr rbegin rend rscore} $todoa($pos) break
				gzclose $f
			}
		}
	}
	close $o; catch {close $allvars} ; catch {close $orivars}; catch {close $sreg}
}

proc pmulticompar_job {args} {
	upvar job_logdir job_logdir
	set cmdline [list cg pmulticompar {*}$args]
	set regonly 0
	set split 1
	set erroronduplicates 0
	set targetvarsfile {}
	set skipincomplete 1
	set optkeepfields *
	set force [job_force_get 1]
	set limitreg {}
	set type var
	cg_options pmulticompar args {
		-r - -reannotregonly {
			putslog "Reannot reg only"
			set regonly $value
		}
		-s - -split {
			set split $value
		}
		-e - -erroronduplicates {
			set erroronduplicates $value
		}
		-t - -targetvarsfile {
			set targetvarsfile $value
		}
		-i - -skipincomplete {
			set skipincomplete $value
		}
		-limitreg {
			set limitreg $value
		}
		-keepfields {
			set optkeepfields $value
		}
		-type {
			if {$value ni "var meth"} {
				error "pmulticompar error: unsupported value \"$value\" for option -type"
			}
			set type $value
		}
		-m - -maxopenfiles {
			maxopenfiles [expr {$value - 4}]
		}
	} compar_file 1
	set dirs $args
# putsvars compar_file dirs regonly split targetvarsfile erroronduplicates
	job_logfile [file dir $compar_file]/[job_relfile2name pmulticompar- $compar_file] [file dir $compar_file] $cmdline \
		{*}[versions dbdir zst os]
	if {[jobfileexists $compar_file]} {
		set dirs [list $compar_file {*}$dirs]
	}
	set compar_file [file_absolute $compar_file]
	set basedir [file dir $compar_file]
	set compar_file_root [gzroot $compar_file]
	unset -nocomplain samplesa
	set files {}
	set samples {}
	foreach dir $dirs {
		set dir [file_absolute $dir]
		if {![jobfileexists $dir]} {error "$dir does not exist (and will not be made by previous jobs))"}
		if {[file isdir $dir]} {
			# a directory is given, look for the variant file(s)
			set sample [file tail $dir]
			set file [lindex [jobglob $dir/fannotvar-$sample.tsv] 0]
			if {$file ne ""} {
				lappend files $file
				if {[info exists samplesa($sample)]} {
					puts stderr "sample \"$sample\" in \"$file\" was already in \"$samplesa($sample)\""
					if {$erroronduplicates} {exit 1} else continue
				}
				set samplesa($sample) $file
				lappend samples $sample
			} else {
				set tempfiles [jobglob $dir/var-*-$sample.tsv]
				if {![llength $tempfiles]} {
					set tempfiles [jobglob $dir/var-$sample.tsv]
				}
				if {![llength $tempfiles]} {error "no variant files found in dir $dir"}
				foreach file $tempfiles {
					set base [file root [file tail [gzroot $file]]]
					set file [gzfile $file]
					set sample [sourcename $base]
					if {[info exists samplesa($sample)]} {
						puts stderr "sample \"$sample\" in \"$file\" was already in \"$samplesa($sample)\""
						if {$erroronduplicates} {exit 1} else continue
					}
					lappend files $file
					set samplesa($sample) $file
					lappend samples $sample
				}
			}
		} else {
			if {[file exists $dir]} {
				set filesamples [cg select -a $dir]
			} else {
				set filesamples {}
			}
			if {[llength $filesamples]} {
				# is a multicompar file
				set file $dir
				set newsample 0
				foreach sample $filesamples {
					if {[info exists samplesa($sample)]} {
						puts stderr "sample \"$sample\" in \"$file\" was already in \"$samplesa($sample)\""
						if {$erroronduplicates} {exit 1} else continue
					}
					# find ori variant file (and thus sampledir)
					set samplefile [pmulticompar_findsample $basedir $sample fannotvar-$sample.tsv var-$sample.tsv var-*-$sample.tsv]
					if {$samplefile eq ""} {
						error "Variant file for sample $sample not found"
					}
					set samplesa($sample) $samplefile
					lappend samples $sample
					set newsample 1
				}
				if {$newsample} {lappend files $file}
			} else {
				set base [file root [file tail [gzroot $dir]]]
				set sample [sourcename $base]
				set file $dir
				if {[info exists samplesa($sample)]} {
					puts stderr "sample \"$sample\" in \"$file\" was already in \"$samplesa($sample)\""
					if {$erroronduplicates} {exit 1} else continue
				}
				lappend files $file
				set samplesa($sample) $file
				lappend samples $sample
			}
		}
	}
	#
	# merge variants
	# todo: check for concurrency
	set workdir [workdir $compar_file]/multicompar
	file mkdir $workdir
	set real_compar_file [jobglob $compar_file]
	if {$real_compar_file ne ""} {
		set allfiles [list $real_compar_file {*}$files]
	} else {
		set allfiles $files
	}
	# analysisinfo
	set deps {}
	foreach file [bsort $files] {
		lappend deps "([analysisinfo_file $file])"
	}
	set target [analysisinfo_file $compar_file]
	job [job_relfile2name multicompar_analysisinfo- $compar_file] -deps $deps -optional 1 \
	-targets {$target} -code {
		set deps [list_remove $deps {}]
		if {[llength $deps]} {
			cg cat -c 0 -m 1 {*}$deps > $target.pmulticompar.temp
			file rename -- $target.pmulticompar.temp $target
		} else {
			file_write $target ""
		}
	}
	# for calculating the varlines needed, we can treat targetvarsfile as just another variant file
	set files $allfiles
	if {$targetvarsfile ne ""} {lappend files $targetvarsfile}
	set allvarsfile $workdir/vars.tsv
	multi_merge_job $allvarsfile $files -limitreg $limitreg -split $split -force $force
	# 
	# add extra var lines to each sample file to get all vars in vars.tsv
	set pastefiles [list $workdir/vars.tsv]
	foreach sample $samples {
		set samplevarsfile $samplesa($sample)
		set sampledir [file dir $samplesa($sample)]
		set sregfile $sampledir/sreg-$sample.tsv
		set varallfile $sampledir/varall-$sample.tsv
		if {![jobfileexists $varallfile]} {
			set varallfile $sampledir/varall-$sample.gvcf
			if {![jobfileexists $varallfile]} {
				# back to original for error
				set varallfile $sampledir/varall-$sample.tsv
			}
		}
		set target $workdir/avars-$sample.tsv.zst
		lappend pastefiles $target
		if {![jobfileexists $sregfile] && ![jobfileexists $varallfile]} {
			set msg "no sorted region file ($sregfile) or varallfile ($varallfile) found: not properly processed sample"
			if {!$skipincomplete} {
				error $msg
			} else {
				putslog "warning: $msg"
			}
		}
		set deps [list $allvarsfile $samplevarsfile]
		if {[jobfileexists $sregfile]} {lappend deps $sregfile}
		if {[jobfileexists $varallfile]} {lappend deps $varallfile}
		# add all possibly usefull files to deps, which ones are actually used is decided in the job
		lappend deps {*}[jobglob $sampledir/*-$sample.bcol]
		lappend deps {*}[jobglob $sampledir/coverage*/*-$sample.bcol]
		lappend deps {*}[jobglob $sampledir/coverage*/*-$sample.tsv]
		lappend deps {*}[jobglob $sampledir/reg_*-$sample.tsv]
		job multicompar_addvars-$sample -force $force -deps $deps -targets {
			$target
		} -vars {
			allvarsfile samplevarsfile sregfile varallfile sample split sampledir optkeepfields limitreg type
		} -code {
			set allvarsfile [gzfile $allvarsfile]
			set samplevarsfile [gzfile $samplevarsfile]
			if {$limitreg ne ""} {
				set tempfile [tempdir]/[gzroot [file tail $samplevarsfile]]
				cg regselect $samplevarsfile $limitreg > $tempfile
				set samplevarsfile $tempfile
			}
			# if the file does not exist, the following will make it empty
			set sregfile [lindex [gzfiles $sregfile] 0]
			set varallfile [lindex [gzfiles $varallfile] 0]
			set vf [gzopen $allvarsfile]
			set vheader [tsv_open $vf]
			gzclose $vf
			if {$vheader ne "chromosome begin end type ref alt"} {error "internal error: index vars.tsv in wrong format"}
			set f [gzopen $samplevarsfile]
			set header [tsv_open $f]
			gzclose $f
			set basicposs [tsv_basicfields $header 6 0]
			set seqpos [lsearch $header sequenced]
			set zygpos [lsearch $header zyg]
			set a1pos [lsearch $header alleleSeq1]
			set a2pos [lsearch $header alleleSeq2]
			set keepfields [list_sub $header -exclude $basicposs]
			if {$type eq "meth"} {
				set mergefields {num_motifs_in_group group_sequence}
				set keepfields [list_lremove $keepfields $mergefields]
			} elseif {[string match fannotvar-* [file tail $samplevarsfile]]} {
				set mergefields {xRef geneId mrnaAcc proteinAcc symbol orientation component componentIndex hasCodingRegion impact nucleotidePos proteinPos annotationRefSequence sampleSequence genomeRefSequence pfam}
				set keepfields [list_lremove $keepfields $mergefields]
			}
			set keepfields [list_remove $keepfields sequenced zyg alleleSeq1 alleleSeq2]
			if {$optkeepfields ne "*"} {
				set keepfields [list_common $keepfields $optkeepfields]
			}
			set keepposs [list_cor $header $keepfields]
			if {$seqpos == -1} {
				lappend reannotheader sequenced-$sample
			}
			set bcolannot {}
			set oldbcolannot {}
			set regfiles {}
			set coverageRefScorefiles {}
			set allfound 1
			set pos 0
			foreach field $keepfields {
				if {$field eq "refscore"} {set field refScore}
				set file $sampledir/${field}-$sample.bcol
				if {[file exists $file]} {
					lappend bcolannot $pos $file
					incr pos ; continue
				}
				set file $sampledir/bcolall/${field}-$sample.bcol
				if {[file exists $file]} {
					lappend bcolannot $pos $file
					incr pos ; continue
				}
				set allfound 0
				set file [lindex [glob -nocomplain $sampledir/coverage*/${field}-*-$sample.bcol] 0]
				if {$file ne ""} {
					lappend oldbcolannot $pos $file
					incr pos ; continue
				}
				set file [lindex [gzfiles $sampledir/reg_${field}-$sample.tsv] 0]
				if {$file ne ""} {
					set f [gzopen $file]
					set header [tsv_open $f]
					gzclose $f
					set regfield [lsearch $header $field]
					if {$regfield == -1} {
						set regfield [lsearch $header score]
						if {$regfield == -1} {
							set regfield [lsearch $header score]
							if {$regfield == -1} {
								set regfield [lsearch $header name]
							}
						}
					}
					lappend regfiles $pos $file $regfield
					incr pos ; continue
				}
				if {$field in {coverage refscore}} {
					set coverageRefScorefiles [lindex [gzfiles $sampledir/coverage*/coverageRefScore-*-$sample.tsv] 0]
				}
				incr pos
			}
			set numbcolannot [expr {[llength $bcolannot]/2}]
			set numregfiles [expr {[llength $regfiles]/3}]
			set newheader [list sequenced-$sample]
			if {$zygpos != -1} {lappend newheader zyg-$sample}
			if {$a1pos != -1} {lappend newheader alleleSeq1-$sample}
			if {$a2pos != -1} {lappend newheader alleleSeq2-$sample}
			foreach field $keepfields {
				lappend newheader ${field}-$sample
			}
			if {[file extension [gzroot $varallfile]] eq ".gvcf"} {
				# puts [list cg vcf2tsv -refout 1 -sort 0 -keepfields $keepfields $varallfile | ../bin/multicompar_addvars $split [join $newheader \t] $allvarsfile $samplevarsfile $sregfile - $numbcolannot $numregfiles {*}$bcolannot {*}$regfiles {*}$keepposs]
				set temp [filetemp $target 1 1]
				exec cg vcf2tsv -refout 1 -sort 0 -keepfields $keepfields $varallfile | multicompar_addvars $split [join $newheader \t] $allvarsfile $samplevarsfile $sregfile - $numbcolannot $numregfiles {*}$bcolannot {*}$regfiles {*}$keepposs {*}[compresspipe $target 1 1] > $temp
				file rename -force -- $temp $target
			} elseif {$varallfile ne "" || $allfound || (![llength $oldbcolannot] && ![llength $coverageRefScorefiles])} {
				# puts stderr [list ../bin/multicompar_addvars $split [join $newheader \t] $allvarsfile $samplevarsfile $sregfile $varallfile $numbcolannot $numregfiles {*}$bcolannot {*}$regfiles {*}$keepposs]
				set temp [filetemp $target 1 1]
				exec multicompar_addvars $split [join $newheader \t] $allvarsfile $samplevarsfile $sregfile $varallfile $numbcolannot $numregfiles {*}$bcolannot {*}$regfiles {*}$keepposs {*}[compresspipe $target 1 1] > $temp
				file rename -force -- $temp $target
			} else {
				set temp [filetemp $target 1 1]
				multicompar_tcl_addvars $sample $temp $split $allvarsfile $samplevarsfile $sregfile $varallfile $bcolannot $oldbcolannot $regfiles $coverageRefScorefiles $keepfields
				file rename -force -- $temp $target
			}
		}
	}
	if {$targetvarsfile ne ""} {
		# add targetvarsfile annotation
		set target $workdir/targets_annot.tsv
		lappend pastefiles $target
		job multicompar_targets -force $force -deps {$allvarsfile $targetvarsfile} -targets {$target} \
		  -vars {allvarsfile targetvarsfile split} -code {
			set targetsfield [lindex [split [file root [file tail $targetvarsfile]] -] end]
			set f [gzopen $targetvarsfile]
			set header [tsv_open $f]
			gzclose $f
			set tempdbposs [tsv_basicfields $header 6 0]
			set dbposs [lrange $tempdbposs 0 2]
			set type2pos [lindex $tempdbposs 3]
			if {$type2pos == -1} {
				error "$targetvarsfile has no type field"
			}
			set alt2pos [lindex $tempdbposs 5]
			if {$alt2pos == -1} {
				error "$targetvarsfile has no alt field"
			}
			set keeppos [lsearch $header name]
			if {$keeppos == -1} {set keeppos {}}
			file_write $target.pmulticompar.temp $targetsfield\n
			# puts stderr [list ../bin/var_annot $allvarsfile 0 1 2 3 5 $targetvarsfile {*}$dbposs $type2pos $alt2pos "" {*}$keeppos]
			exec var_annot $allvarsfile 0 1 2 3 5 $targetvarsfile {*}$dbposs $type2pos $alt2pos "" {*}$keeppos >> $target.pmulticompar.temp 2>@ stderr
			file rename -force -- $target.pmulticompar.temp $target
		}
	}
	#
	# paste all together for final multicompar
	# putsvars pastefiles
	set endcommand [list file delete {*}$pastefiles]
	set endcommand {}
	if {$type eq "meth"} {
		append endcommand "set temp \[filetemp_ext [list $compar_file]\]\n"
		append endcommand "cg select -overwrite 1 -rf {type ref alt sequenced-*} [list $compar_file] \$temp\n"
		append endcommand "file rename -force \$temp [list $compar_file]\n"
	}
	tsv_paste_job -o $compar_file -forcepaste 1 -endcommand $endcommand {*}$pastefiles
}

proc cg_pmulticompar {args} {
	set args [job_init {*}$args]
	pmulticompar_job {*}$args
	job_wait
}
