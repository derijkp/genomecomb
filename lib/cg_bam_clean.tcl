proc bam_clean_job {args} {
	set sort 1
	set removeduplicates 1
	set realign 1
	set realignopts {}
	set realigndeps {}
	set clipamplicons {}
	set cleanup 1
	set regionfile {}
	set threads 2
	cg_options process_sample args {
		-sort {set sort $value}
		-removeduplicates {set removeduplicates $value}
		-realign {set realign $value}
		-realignopts {set realignopts $value}
		-realigndeps {set realigndeps $value}
		-clipamplicons {set clipamplicons $value}
		-threads {
			set threads $value
			# not used yet
		}
		-cleanup {set cleanup $value}
		-regionfile {set regionfile $value}
	} {bamfile refseq sample}
	if {$regionfile ne ""} {
		lappend realigndeps $regionfile
	}
	if {[llength $args]} {
		array set opt $args
	}
	upvar job_logdir job_logdir
	set dir [file dir $bamfile]
	set file [file tail $bamfile]
	set pre [lindex [split $file -] 0]
	set root [file_rootname $file]
	# make gatk refseq
	set gatkrefseq [gatk_refseq_job $refseq]
	set dict [file root $gatkrefseq].dict
	# precalc skips and cleanuplist
	set skips {}
	set cleanuplist {}
	set temproot $root
	if {$sort} {
		lappend cleanuplist $dir/$pre-$temproot.bam $dir/$pre-$temproot.bam.bai
		set temproot s$temproot
		lappend skips -skip [list $dir/$pre-$temproot.bam $dir/$pre-$temproot.bam.analysisinfo]
	}
	if {$removeduplicates ne "0"} {
		lappend cleanuplist $dir/$pre-$temproot.bam $dir/$pre-$temproot.bam.bai
		set temproot d$temproot
		lappend skips -skip [list $dir/$pre-$temproot.bam $dir/$pre-$temproot.bam.analysisinfo]
	}
	if {$realign ne "0"} {
		lappend cleanuplist $dir/$pre-$temproot.bam $dir/$pre-$temproot.bam.bai
		set temproot r$temproot
		lappend skips -skip [list $dir/$pre-$temproot.bam $dir/$pre-$temproot.bam.analysisinfo]
	}
	if {$clipamplicons ne ""} {
		lappend cleanuplist $dir/$pre-$temproot.bam $dir/$pre-$temproot.bam.bai
		set temproot c$temproot
		lappend skips -skip [list $dir/$pre-$temproot.bam $dir/$pre-$temproot.bam.analysisinfo]
	}
	# start jobs
	# sort using default
	if {$sort} {
		list_pop skips 0; list_pop skips 0
		if {$sort == 2} {
			bam_sort_job -method alreadysorted {*}$skips $bamfile $dir/$pre-s$root.bam
		} else {
			bam_sort_job {*}$skips $bamfile $dir/$pre-s$root.bam
		}
		set root s$root
	}
	if {$removeduplicates ne "0"} {
		list_pop skips 0; list_pop skips 0;
		bam_markduplicates_job {*}$skips $dir/$pre-$root.bam $dir/$pre-d$root.bam
		set root d$root
	}
	# index intermediate result
	job bamindex-$pre-$root -optional 1 -deps {$dir/$pre-$root.bam} -targets {$dir/$pre-$root.bam.bai} {*}$skips -code {
		exec samtools index $dep >@ stdout 2>@ stderr
		puts "making $target"
	}
	if {$realign ne "0"} {
		if {[isint $regionfile]} {
			# extract regions with coverage >= $regionfile (for cleaning)
			set regionfile [bam2reg_job -mincoverage $regionfile \
				-distrreg $distrreg -refseq $refseq \
				-skip [list $resultbamfile $resultbamfile.analysisinfo] \
				$sampledir/map-${aligner}-$sample.bam]
		}

		list_pop skips 0; list_pop skips 0;
		# realign around indels
		set deps [list $dir/$pre-$root.bam $dir/$pre-$root.bam.bai $dict $gatkrefseq $refseq {*}$realigndeps]
		if {$realign eq "srma"} {
			set srma [srma]
			job bamrealign-$root -deps $deps \
			-targets {$dir/$pre-r$root.bam $dir/$pre-r$root.bam.analysisinfo} \
			-vars {gatkrefseq refseq srma pre realignopts regionfile} {*}$skips -code {
				analysisinfo_write $dep $target realign srma realign_version [version srma]
				if {$regionfile ne ""} {
					set bedfile [tempbed $regionfile $refseq]
					lappend realignopts RANGES=$bedfile
				}
				exec java -XX:ParallelGCThreads=1 -jar $srma I=$dep O=$target.temp R=$gatkrefseq {*}$realignopts 2>@ stderr >@ stdout
				catch {file rename -force $target.temp.bai $target.bai}
				catch {file delete $target.intervals}
				file rename -force $target.temp $target
			}
		} else {
			job bamrealign-$root -mem [job_mempercore 10G 2] -cores 2 -deps $deps \
			-targets {$dir/$pre-r$root.bam $dir/$pre-r$root.bam.analysisinfo} \
			-vars {gatkrefseq refseq gatk realignopts regionfile} {*}$skips -code {
				analysisinfo_write $dep $target realign gatk realign_version [version gatk3]
				if {$regionfile ne ""} {
					set bedfile [tempbed $regionfile $refseq]
					lappend realignopts -L $bedfile
				}
				gatk3exec {-XX:ParallelGCThreads=1 -Xms512m -Xmx8g} RealignerTargetCreator -R $gatkrefseq -I $dep -o $target.intervals {*}$realignopts 2>@ stderr >@ stdout
				if {[loc_compare [version gatk3] 2.7] >= 0} {
					set extra {--filter_bases_not_stored}
				} else {
					set extra {}
				}
				lappend extra --filter_mismatching_base_and_quals
				gatk3exec {-XX:ParallelGCThreads=1 -Xms512m -Xmx8g} IndelRealigner -R $gatkrefseq \
					-targetIntervals $target.intervals -I $dep \
					-o $target.temp {*}$extra 2>@ stderr >@ stdout
				catch {file rename -force $target.temp.bai $target.bai}
				catch {file delete $target.intervals}
				file rename -force $target.temp $target
			}
		}
		set root r$root
		job bamrealign_index-$root -optional 1 -deps {$dir/$pre-$root.bam} {*}$skips -targets {$dir/$pre-$root.bam.bai} -code {
			exec samtools index $dep >@ stdout 2>@ stderr
			puts "making $target"
		}
	}
	if {$clipamplicons ne ""} {
		list_pop skips 0; list_pop skips 0;
		job bamclean_clipamplicons-$root -deps {
			$dir/$pre-$root.bam $clipamplicons
		} -targets {
			$dir/$pre-c$root.bam $dir/$pre-c$root.bam.analysisinfo
		} {*}$skips -vars {clipamplicons} -code {
			analysisinfo_write $dep $target clipamplicons genomecomb clipamplicons_version [version genomecomb] clipampliconsfile [filename $clipamplicons]
			cg sam_clipamplicons $dep2 $dep $target.temp
			file rename $target.temp $target
		}
		set root c$root
		job bamclean_clipamplicons_index-$root -optional 1 -deps {$dir/$pre-$root.bam} -targets {$dir/$pre-$root.bam.bai} -code {
			exec samtools index $dep >@ stdout 2>@ stderr
			puts "making $target"
		}
	}
	if {$cleanup} {
		cleanup_job bamclean_remtemp-$root $cleanuplist [list $dir/$pre-$root.bam]
	}
	return $dir/$pre-$root.bam
}

