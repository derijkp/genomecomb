#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc process_compare_checkfield {file field} {
	set f [open $file]
	set header [tsv_open $f]
	close $f
	inlist $header $field
}

proc process_indexcompress {file} {
	set ext [file extension $file]
	if {$ext eq ".gz"} {
		gunzip $file
		set file [file root $file]
	}
	set f [gzopen $file]
	set header [tsv_open $f]
	foreach field {offset end1} {
		set fpos [lsearch $header $field]
		if {$fpos != -1} break
	}
	if {$fpos == -1} {error "no column offset, end1 in file $file"}
	if {([gets $f] eq "") && [eof $f]} return
	gzclose $f
	if {![file exists $file.${field}_index]} {
		tsv_index $field $file
	}
	putslog "Compressing $file"
	if {![inlist {.rz .bgz} $ext]} {
		exec bgzip -c $file > $file.gz.temp
		file rename -force $file.gz.temp $file.gz
		file delete $file
	}
}

proc cg_process_indexcompress {args} {
	global scriptname action
	if {[llength $args] != 1} {
		error "format is: $scriptname $action file\n - makes index, and compresses to bgzip"
	}
	foreach {file} $args break
	process_indexcompress $file
}

proc process_sample {args} {
	set pos 0
	set split 1
	foreach {key value} $args {
		switch -- $key {
			-split - -s {
				set split $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	set srcdir ""
	if {[llength $args] == 1} {
		foreach {destdir} $args break
		set workdir [file_absolute $destdir]
		putslog "Making $workdir"
		file mkdir $workdir
	} elseif {[llength $args] == 2} {
		foreach {srcdir destdir} $args break
		set workdir [file_absolute $destdir]
		putslog "Making $workdir"
		file mkdir $workdir
		set srcdir [file_absolute $srcdir]
		if {[file exists $srcdir]} {
			set srcdir [file_absolute $srcdir]
			file delete $workdir/ori
			cplinked_file $srcdir $workdir/ori
		}
	} else {
		errorformat process_sample
	}
	set sample [file tail $workdir]
	# process_bam2cg $srcdir
	set keepdir [pwd]
	cd $workdir
	job_logdir $workdir/log_jobs
	set chromosomes {}
	set files [gzfiles $workdir/ori/ASM/REF/coverage*.tsv]
	if {[llength $files]} {
		foreach file $files {
			lappend chromosomes [chr_clip [lindex [split [file tail $file] -] 1]]
		}
		set chromosomes [ssort -natural [list_remdup $chromosomes]]
	} else {
		set files [gzfiles $workdir/coverage/coverage*.tsv]
		foreach file $files {
			lappend chromosomes [chr_clip [lindex [split [file tail $file] -] end-1]]
		}
		set chromosomes [ssort -natural [list_remdup $chromosomes]]
	}

	# start from CGI data
	job cg_svar-$sample -optional 1 -deps {ori/ASM/var-*-ASM*.tsv} -targets {svar-$sample.tsv} \
	-skip {fannotvar-$sample.tsv sreg-$sample.tsv reg_refcons-$sample.tsv reg_nocall-$sample.tsv reg_cluster-$sample.tsv reg_ns-$sample.tsv reg_lowscore-$sample.tsv} \
	-code {
		set varfile $dep
		set f [gzopen $varfile]
		set info {}
		while {![eof $f]} {
			set line [gets $f]
			if {[string index $line 0] ne "#"} break
			lappend info $line
		}
		catch {gzclose $f}
		if {[file exists info.txt]} {
			set test [split [file_read info.txt] \n]
			if {$info ne $test} {
				error "$target1 already has info.txt that differs from data in the source $dir"
			}
		}
		file_write info.txt [join $info \n]
		putslog "Sort var file ($varfile)"
		cg select -s "chromosome begin end varType" $varfile $target.temp
		file rename -force $target.temp $target
	}

	job cg_sgene-$sample -optional 1 {ori/ASM/gene-*-ASM*.tsv} {sgene-$sample.tsv} {
		set genefile $dep
		if {[llength $genefile] != 1} {error "could not identify genefile"}
		putslog "Sort gene file ($genefile)"
		cg select -s "chromosome begin end" $genefile $target.temp
		file rename -force $target.temp $target
	}

	# annotated vars file
	job cg_annotvar-$sample -optional 1 -vars {split} \
	-deps {svar-$sample.tsv (sgene-$sample.tsv)} -targets {annotvar-$sample.tsv} \
	-skip {fannotvar-$sample.tsv reg_cluster-$sample.tsv reg_ns-$sample.tsv reg_lowscore-$sample.tsv} -code {
		putslog "Create annotated varfile $target"
		if {[file exists $dep2]} {
			cg cg2tsv -split $split -sorted 1 $dep1 $dep2 $target.temp
		} else {
			cg cg2tsv -split $split -sorted 1 $dep1 $target.temp
		}
		file rename -force $target.temp $target
	}

	# if not from CGI, we do not have svar and sgene, take from first var_* file found
	job cg_annotvar_var-$sample -optional 1 {vars_*.tsv} {annotvar-$sample.tsv} {
		gzmklink $dep $target
	}

	# if we also do not find a var_* file, take from first variant* file found
	job cg_annotvar_variant-$sample -optional 1 {variant*.tsv} {annotvar-$sample.tsv} {
		set file [gzfile $dep]
		gzmklink $file $target
	}

	# only if reg file exists, otherwise extract from svar (next)
	job cg_sreg-$sample -optional 1 {ori/ASM/reg-*-ASM*.tsv} {sreg-$sample.tsv} {
		set regfile [gzfile $dep]
		putslog "Sort region file ($regfile)"
		cg select -s "chromosome begin end" $regfile $target.temp
		file rename -force $target.temp $target
	}

	job cg_regfromsvar-$sample -optional 1 {svar-$sample.tsv} {sreg-$sample.tsv} {
		set svarfile $dep
		putslog "Extract $target from $svarfile"
		cg select -q {$varType != "no-call" && $varType != "no-ref"} -f "chromosome begin end" $svarfile $target.temp
		cg regjoin $target.temp > $target.temp2
		file rename -force $target.temp2 $target
		file delete $target.temp
	}

	# multiarch
	# if we do not have an svar (CGI), try getting a region file from the coverage
	# currently hardcoded at coverage > 7
	job bam_regsfromscoverage-$sample -optional 1 {coverage/coverage-*-coverage-*.bcol} {sreg-$sample.tsv} {
		set files [ssort -natural $deps]
		cg regextract -above 1 7 {*}$files > $target.temp
		cg select -s {chromosome begin end} $target.temp $target.temp2
		file rename -force $target.temp2 $target
		file delete $target.temp
	}

	job cg_reg_refcons-$sample -optional 1 {svar-$sample.tsv} {reg_refcons-$sample.tsv} {
		putslog "Find refcons regions for $dep"
		cg refconsregions $dep > $target.temp
		file rename -force $target.temp $target
	}

	job cg_reg_nocall-$sample -optional 1 {svar-$sample.tsv} {reg_nocall-$sample.tsv} {
		putslog "Find partial no-call regions for dep"
		if {[catch {
			nocallregions $dep $target.temp
		}]} {
			puts stderr "Could not make $job (old version files ?)"
		} else {
			file rename -force $target.temp $target
		}
	}

	set files [jobglob ori/ASM/REF/coverage*-chr*]
	if {[llength $files]} {
		# this will only work if ori/ASM/REF/coverage*-chr* already exist from the start
		# maybe later make more flexible
		set fields [list_remove [cg select -h [lindex $files 0]] offset uniqueSequenceCoverage coverage]
		lappend fields coverage
		set files [glob -nocomplain ori/ASM/REF/coverage*-chr*-*]
		foreach file $files {
			regexp {^ori/ASM/REF/coverage.*-chr([^-]*)-.*$} $file temp chr
			job cg_coverage-$sample-$chr -deps $file -vars {sample chr} \
			   -targets {coverage/$_fields-$chr-$sample.bcol} -code {
				# make coverage files
				set file $dep
				file mkdir coverage
				set chr [chr_clip $chr]
				set header [cg select -h $file]
				foreach posfield {offset pos} {
					if {[lsearch $header $posfield] != -1}  break
				}
				if {$posfield == -1} {
					exiterror "No position/offset field found in $file"
				}
				foreach covfield {uniqueSequenceCoverage coverage} {
					if {[lsearch $header $covfield] != -1}  break
				}
				if {$covfield == -1} {
					exiterror "No coverage/uniqueSequenceCoverage field found in $file"
				}
				set other [list_remove $header $posfield $covfield]
				foreach field $other {
					set base coverage/$field-$chr-$sample
					if {![file exists $base.bcol]} {
						putslog "Making $base.bcol"
						if {[catch {
							exec [catprog $file] $file | cg bcol make -p $posfield -t s -n -1 $base $field
						} e]} {
							exec [catprog $file] $file | cg bcol make -p $posfield -t i -n -1 $base $field >@ stdout 2>@ stderr
						}
					}
				}
				set base coverage/coverage-$chr-$sample
				if {![file exists $base.bcol]} {
					putslog "Making $base.bcol"
						if {[catch {
							exec [catprog $file] $file | cg bcol make -p $posfield -t su $base $covfield
						} e]} {
							exec [catprog $file] $file | cg bcol make -p $posfield -t iu $base $covfield >@ stdout 2>@ stderr
						}
				}
			}
		}
		foreach field $fields {
			set target coverage/bcol_[string tolower $field]-$sample.tsv
			job cg_bcol_coverage-$field-$sample -deps {coverage/$field-*-$sample.bcol} -vars {sample field} \
			-targets $target -code {
				set deps [ssort -natural $deps]
				set f [open $target.temp w]
				puts $f "chromosome\tfile"
				foreach dep $deps {
					set tail [file tail $dep]
					regexp "^$field-(\[^-\]+)" $tail temp chr
					puts $f $chr\t[file tail $dep]
				}
				close $f
				file rename -force $target.temp $target
			}
		}
	}

	# multiarch
	# if we are coming from bams, coverage file name looks different, use these by making link
	job cg_bamcoverage-$sample -optional 1 {^coverage/coverage-(.*)-coverage-(.*)\.bcol ^coverage/coverage-(.*)-coverage-(.*)\.bcol\.bin$} {coverage/coverage-\1-\2.bcol coverage/coverage-\1-\2.bcol.bin} {
		gzmklink [lindex $dep 0] [lindex $target 0]
		gzmklink [lindex $dep 1] [lindex $target 1]
	}

	job cg_cpSV-$sample -optional 1 {ori/ASM/SV ^ori/ASM/SV/(.*)$} {SV/\_} {
		putslog "Copying SV"
		set targetdir [file dir $target]
		file delete -force $targetdir
		file delete -force $targetdir.temp
		file copy $dep $targetdir.temp
		file rename -force $targetdir.temp $targetdir
	}

	job cg_cgsv-$sample -optional 1 {SV/allJunctionsBeta-*.tsv} {cgsv-$sample.tsv} {
		cg convcgsv $dep $target
	}

	job cg_cgsv_alpha-$sample -optional 1 {SV/annotatedJunctionsAlpha-*.tsv} {cgsv-$sample.tsv} {
		cg convcgsv $dep $target
	}

	job cg_cpCNV-$sample -optional 1 {ori/ASM/CNV ^ori/ASM/CNV/(.*)$} {CNV/\_} {
		putslog "Copying CNV"
		set targetdir [file dir $target]
		file delete -force $targetdir
		file delete -force $targetdir.temp
		file copy $dep $targetdir.temp
		file rename -force $targetdir.temp $targetdir
	}

	job cg_cgcnv -optional 1 {CNV/cnvSegmentsBeta-*.tsv} {cgcnv-$sample.tsv} {
		cg convcgcnv $dep $target
	}

	job cg_cgcnv_diploid-$sample -optional 1 {CNV/cnvSegmentsDiploidBeta-*.tsv} {cgcnv-$sample.tsv} {
		cg convcgcnv $dep $target
	}

	job cg_cgcnv_alpha-$sample -optional 1 {CNV/cnvSegmentsAlpha-*.tsv} {cgcnv-$sample.tsv} {
		cg convcgcnv [gzfile $dep] $target
	}

	# multiarch
	job reg_cluster-$sample -optional 1 {annotvar-$sample.tsv} {reg_cluster-$sample.tsv} {
		cg clusterregions < $dep > $target.temp
		file rename -force $target.temp $target
	}

	job reg_ns-$sample -optional 1 {annotvar-$sample.tsv} {reg_ns-$sample.tsv} {
		putslog "Find regions with N's for $dep"
		cg select -f {chromosome begin end} -q {$alleleSeq1 ~ /[N?]/ || $alleleSeq2 ~ /[N?]/} < $dep > $target.temp
		file rename -force $target.temp $target
	}

	job reg_lowscore-$sample -optional 1 {annotvar-$sample.tsv} {reg_lowscore-$sample.tsv} {
		set header [cg select -h $dep]
		if {[llength [list_common $header {totalScore1 totalScore2}]] == 2} {
			putslog "Find regions with lowscores for $dep"
			cg select -f {chromosome begin end} -q {$totalScore1 < 60 || $totalScore2 < 60} < $dep > $target.temp
			file rename -force $target.temp $target
		}
	}

	job cg_fannotvar-$sample -optional 1 {annotvar-$sample.tsv (reg_refcons-$sample.tsv) (reg_cluster-$sample.tsv) (coverage/bcol_coverage-$sample.tsv) (coverage/bcol_refscore-$sample.tsv)} {fannotvar-$sample.tsv} {
		cg annotate $dep $target {*}[list_remove [lrange $deps 1 end] {}]
	}

	job cg_multitechlink_var-$sample -optional 1 {fannotvar-$sample.tsv} {var-cg-cg-$sample.tsv} {
		gzmklink $dep $target
	}

	job cg_multitechlink_sreg-$sample -optional 1 {sreg-$sample.tsv} {sreg-cg-cg-$sample.tsv} {
		gzmklink $dep $target
	}

	job cg_multitechlink_coverage-$sample -optional 1 {coverage} {coverage-cg-$sample} {
		gzmklink $dep $target
	}

	job reg_covered-$sample {sreg-$sample.tsv} -optional 1 {reg-$sample.covered} {
		putslog "Genomic coverage of sequenced regions"
		cg covered $dep > $target.temp
		file rename -force $target.temp $target
	}

	job cg_filteredrefcons-$sample -optional 1 -vars sample {sreg-$sample.tsv reg_refcons-$sample.tsv} {filteredrefcons-$sample.tsv filteredrefcons-$sample.covered} {
		putslog "Coverage of refcons region"
		cg regsubtract $dep1 $dep2 > $target1.temp
		file rename -force $target1.temp $target1
		cg covered $target1 > $target2.temp
		file rename -force $target2.temp $target2
	}

	job cg_filteredns-$sample -optional 1 {sreg-$sample.tsv reg_ns-$sample.tsv} {filteredns-$sample.tsv} {
		putslog "Coverage of ns region"
		cg regsubtract $dep1 $dep2 > $target.temp
		file rename -force $target.temp $target
	}

	job cg_filteredns_covered-$sample -optional 1 {filteredns-$sample.tsv} {filteredns-$sample.covered} {
		putslog "Making $target"
		cg covered $dep > $target.temp
		file rename -force $target.temp $target
	}

	job cg_filteredlowscore-$sample -optional 1 {sreg-$sample.tsv reg_lowscore-$sample.tsv} {filteredlowscore-$sample.tsv} {
		cg regsubtract $dep1 $dep2 > $target.temp
		file rename -force $target.temp $target
	}

	job cg_filteredlowscore_covered-$sample -optional 1 {filteredlowscore-$sample.tsv} {filteredlowscore-$sample.covered} {
		cg covered $dep > $target.temp
		file rename -force $target.temp $target
	}

	job cg_refcons_histo-$sample -optional 1 {reg_refcons-$sample.tsv} {histo-refcons-$sample.tsv} {
		putslog "Making $target"
		cg reghisto $dep > $target.temp
		file rename -force $target.temp $target
	}

	job cg_filteredcluster-$sample -optional 1 {sreg-$sample.tsv reg_cluster-$sample.tsv} {filteredcluster-$sample.tsv} {
		putslog "Coverage of clusters region"
		cg regsubtract $dep1 $dep2 > $target.temp
		file rename -force $target.temp $target
	}

	job cg_filteredcluster_covered-$sample -optional 1 {filteredcluster-$sample.tsv} {filteredcluster-$sample.covered} {
		putslog "Making $target"
		cg covered $dep > $target.temp
		file rename -force $target.temp $target
	}
	job cg_process_cleanup-$sample -optional 1 -deps {(svar-$sample.tsv) (annotvar-$sample.tsv) fannotvar-$sample.tsv sreg-$sample.tsv reg_refcons-$sample.tsv reg_nocall-$sample.tsv reg_cluster-$sample.tsv reg_ns-$sample.tsv reg_lowscore-$sample.tsv} \
		-vars {sample} -rmtargets {svar-$sample.tsv annotvar-$sample.tsv} -code {
			catch {file delete svar-$sample.tsv}
			catch {file delete annotvar-$sample.tsv}
	}
	job cg_process_summary-$sample -deps {
		annotvar-$sample.tsv
		sreg-$sample.tsv
		reg-$sample.covered
		coverage/coverage-$_chromosomes-$sample.bcol
		fannotvar-$sample.tsv
		(reg_refcons-$sample.tsv)
		(reg_nocall-$sample.tsv)
		(SV) (CNV)
		(cgsv-$sample.tsv)
		(cgcnv-$sample.tsv)
		(reg_cluster-$sample.tsv)
		(reg_ns-$sample.tsv)
		(reg_lowscore-$sample.tsv)
		(reg-$sample.covered)
		(filteredrefcons-$sample.covered)
		(filteredns-$sample.covered)
		(filteredlowscore-$sample.covered)
		(filteredcluster-$sample.covered)
		(histo-refcons-$sample.tsv)
	} -targets {summary-$sample.txt} -vars {sample} -code {
		set f [open $target.temp w]
		puts $f "finished\t\t[job_timestamp]"
		set c [split [exec cg select -f chromosome $dep | uniq -c] \n]
		puts $f ""
		puts $f "Variants:\nchromosome\tvalue"
		set total 0
		foreach line [lrange $c 1 end] {
			foreach {num chr} $line break
			puts $f "$chr\t$num"
			incr total $num
		}
		puts $f "total\t$total"
		puts $f ""
		puts $f "covered:\nchromosome\tvalue"
		set c [split [file_read $dep3] \n]
		set total 0
		foreach line [lrange $c 1 end] {
			if {![llength $line]} break
			foreach {chr num} $line break
			puts $f "$chr\t$num"
			incr total $num
		}
		puts $f "total\t$total"
		set header [cg select -h $dep]
		if {[inlist $header impact]} {
			puts $f ""
			puts $f "impact:\nimpact\tnumber_of_variants"
			set c [split [exec cg select -s impact -f impact $dep | uniq -c] \n]
			foreach line [lrange $c 1 end] {
				foreach {num impact} $line break
				puts $f "$impact\t$num"
			}
		}
		close $f
		file rename -force $target.temp $target
	}
	cd $keepdir
}

proc cg_process_sample {args} {
	set args [job_init {*}$args]
	process_sample {*}$args
	job_wait
}
