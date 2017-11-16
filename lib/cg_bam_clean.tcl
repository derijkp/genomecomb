proc bam_clean_job {args} {
	set sort 1
	set removeduplicates 1
	set realign 1
	set realignopts {}
	set realigndeps {}
	set clipamplicons {}
	set cleanup 1
	set regionfile {}
	cg_options process_sample args {
		-sort {set sort $value}
		-removeduplicates {set removeduplicates $value}
		-realign {set realign $value}
		-realignopts {set realignopts $value}
		-realigndeps {set realigndeps $value}
		-clipamplicons {set clipamplicons $value}
		-cleanup {set cleanup $value}
		-regionfile {set regionfile $value}
	} {bamfile refseq sample}
	if {$regionfile ne ""} {
		lappend realigndeps $regionfile
	}
	if {[llength $args]} {
		array set opt $args
	}
	set gatk [gatk]
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
		lappend skips -skip [list $dir/$pre-$temproot.bam]
	}
	if {$removeduplicates} {
		lappend cleanuplist $dir/$pre-$temproot.bam $dir/$pre-$temproot.bam.bai
		set temproot d$temproot
		lappend skips -skip [list $dir/$pre-$temproot.bam]
	}
	if {$realign ne "0"} {
		lappend cleanuplist $dir/$pre-$temproot.bam $dir/$pre-$temproot.bam.bai
		set temproot r$temproot
		lappend skips -skip [list $dir/$pre-$temproot.bam]
	}
	if {$clipamplicons ne ""} {
		lappend cleanuplist $dir/$pre-$temproot.bam $dir/$pre-$temproot.bam.bai
		set temproot c$temproot
		lappend skips -skip [list $dir/$pre-$temproot.bam]
	}
	# start jobs
	# sort using default
	bam_sort_job {*}$skips $bamfile $dir/$pre-s$root.bam
	set root s$root
	list_pop skips 0; list_pop skips 0;
	if {$removeduplicates eq "picard"} {
		list_pop skips 0; list_pop skips 0;
		job bamremdup-$root -mem 7G -cores 2 -deps {$dir/$pre-$root.bam} \
		-targets {$dir/$pre-d$root.bam $dir/$pre-d$root.bam.analysisinfo} \
		-vars {sample} {*}$skips -code {
			analysisinfo_write $dep $target removeduplicates picard removeduplicates_version [version picard]
			puts "removing duplicates"
			file mkdir [scratchdir]/picard
			picard MarkDuplicates	I=$dep	O=$target.temp METRICS_FILE=$target.dupmetrics TMP_DIR=[scratchdir]/picard 2>@ stderr >@ stdout
			file rename -force $target.temp $target
		}
		set root d$root
	} elseif {$removeduplicates} {
		list_pop skips 0; list_pop skips 0;
		job bamremdup-$root -deps {$dir/$pre-$root.bam} \
		-targets {$dir/$pre-d$root.bam $dir/$pre-d$root.bam.analysisinfo} \
		-vars {sample} {*}$skips -code {
			analysisinfo_write $dep $target removeduplicates biobambam removeduplicates_version [version biobambam]
			biobambam bammarkduplicates2 I=$dep	O=$target.temp M=$target.dupmetrics rmdup=0 markthreads=1 tmpfile=[scratchfile] 2>@ stderr >@ stdout
			file rename -force $target.temp $target
		}
		set root d$root
	}
	# index intermediate result
	job bamindex-$pre-$root -optional 1 -deps {$dir/$pre-$root.bam} -targets {$dir/$pre-$root.bam.bai} {*}$skips -code {
		exec samtools index $dep >@ stdout 2>@ stderr
		puts "making $target"
	}
	if {$realign ne "0"} {
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
			job bamrealign-$root -mem 10G -cores 2 -deps $deps \
			-targets {$dir/$pre-r$root.bam $dir/$pre-r$root.bam.analysisinfo} \
			-vars {gatkrefseq refseq gatk realignopts regionfile} {*}$skips -code {
				analysisinfo_write $dep $target realign gatk realign_version [version gatk]
				if {$regionfile ne ""} {
					set bedfile [tempbed $regionfile $refseq]
					lappend realignopts -L $bedfile
				}
				exec [gatkjava] -XX:ParallelGCThreads=1 -Xms512m -Xmx8g -jar $gatk -T RealignerTargetCreator -R $gatkrefseq -I $dep -o $target.intervals {*}$realignopts 2>@ stderr >@ stdout
				if {[loc_compare [version gatk] 2.7] >= 0} {
					set extra {--filter_bases_not_stored}
				} else {
					set extra {}
				}
				lappend extra --filter_mismatching_base_and_quals
				exec [gatkjava] -XX:ParallelGCThreads=1 -Xms512m -Xmx8g -jar $gatk -T IndelRealigner -R $gatkrefseq \
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

