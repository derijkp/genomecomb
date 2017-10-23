proc bam_clean_job {args} {
	oargs bam_clean_job {bamfile refseq sample
		{sort 1}
		{removeduplicates 1}
		{realign 1}
		{realignopts {}}
		{realigndeps {}}
		{clipamplicons {}}
		{cleanup 1}
		{regionfile {}}
		args
	} $args
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
	set root [join [lrange [split [file root $file] -] 1 end] -]
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
	bam_sort_job -infostring infostring {*}$skips $bamfile $dir/$pre-s$root.bam
	set root s$root
	list_pop skips 0; list_pop skips 0;
	if {$removeduplicates eq "picard"} {
		list_pop skips 0; list_pop skips 0;
		job bamremdup-$root -mem 7G -cores 2 -deps {$dir/$pre-$root.bam} -targets {$dir/$pre-d$root.bam} \
		-vars {sample} {*}$skips -code {
			puts "removing duplicates"
			file mkdir [scratchdir]/picard
			picard MarkDuplicates	I=$dep	O=$target.temp METRICS_FILE=$target.dupmetrics TMP_DIR=[scratchdir]/picard 2>@ stderr >@ stdout
			file rename -force $target.temp $target
		}
		set root d$root
	} elseif {$removeduplicates} {
		list_pop skips 0; list_pop skips 0;
		job bamremdup-$root -deps {$dir/$pre-$root.bam} -targets {$dir/$pre-d$root.bam} \
		-vars {sample} {*}$skips -code {
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
			job bamrealign-$root -deps $deps -targets {$dir/$pre-r$root.bam} \
			-vars {gatkrefseq refseq srma pre realignopts regionfile} {*}$skips -code {
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
			job bamrealign-$root -mem 10G -cores 2 -deps $deps -targets {$dir/$pre-r$root.bam} \
			-vars {gatkrefseq refseq gatk pre realignopts regionfile} {*}$skips -code {
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
		job bamclean_clipamplicons-$root -deps {$dir/$pre-$root.bam $clipamplicons} -targets {$dir/$pre-c$root.bam} {*}$skips -code {
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

