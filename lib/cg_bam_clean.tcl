proc bam_clean_job {args} {
	oargs bam_clean_job {bamfile refseq sample
		{removeduplicates 1}
		{realign 1}
		{realignopts {}}
		{realigndeps {}}
		{clipamplicons {}}
		{cleanup 1}
		{bed {}}
		args
	} $args
	set regionfile {}
	if {$bed ne ""} {
		set regionfile $bed
		lappend realigndeps $bed
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
	lappend cleanuplist $dir/$pre-$root.bam $dir/$pre-$root.bam.bai
	set temproot s$root
	lappend skips -skip [list $dir/$pre-$temproot.bam]
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
	# sort using picard
	job bamsort-$root -deps {$bamfile} -targets {$dir/$pre-s$root.bam} \
	-vars {removeduplicates sample} {*}$skips -code {
		file delete $target.temp
		picard SortSam	I=$dep	O=$target.temp	SO=coordinate 2>@ stderr >@ stdout
		file rename -force $target.temp $target
	#	# picard AddOrReplaceReadGroups	I=$src	O=$target.temp3	RGID=$sample	RGLB=solexa-123	RGPL=illumina	RGPU=$sample RGSM=$sample 2>@ stderr >@ stdout
	#	# file delete $target.temp $target.temp2
	#	# file rename -force $target.temp3 $target
	}
	set root s$root
	list_pop skips 0; list_pop skips 0;
	if {$removeduplicates} {
		list_pop skips 0; list_pop skips 0;
		job bamremdup-$root -deps {$dir/$pre-$root.bam} -targets {$dir/$pre-d$root.bam} \
		-vars {sample} {*}$skips -code {
			puts "removing duplicates"
			picard MarkDuplicates	I=$dep	O=$target.temp METRICS_FILE=$target.dupmetrics 2>@ stderr >@ stdout
			file rename -force $target.temp $target
		}
		set root d$root
	}
	# index intermediate result
	job bamindex-$pre-$root -optional 1 -deps $dir/$pre-$root.bam -targets $dir/$pre-$root.bam.bai {*}$skips -code {
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
				exec java -jar $srma I=$dep O=$target.temp R=$gatkrefseq {*}$realignopts 2>@ stderr >@ stdout
				catch {file rename -force $target.temp.bai $target.bai}
				catch {file delete $target.intervals}
				file rename -force $target.temp $target
			}
		} else {
			job bamrealign-$root -deps $deps -targets {$dir/$pre-r$root.bam} \
			-vars {gatkrefseq refseq gatk pre realignopts regionfile} {*}$skips -code {
				if {$regionfile ne ""} {
					set bedfile [tempbed $regionfile $refseq]
					lappend realignopts -L $bedfile
				}
				exec [gatkjava] -jar $gatk -T RealignerTargetCreator -R $gatkrefseq -I $dep -o $target.intervals {*}$realignopts 2>@ stderr >@ stdout
				if {[loc_compare [version gatk] 2.7] >= 0} {
					set extra {--filter_bases_not_stored}
				} else {
					set extra {}
				}
				lappend extra --filter_mismatching_base_and_quals
				exec [gatkjava] -jar $gatk -T IndelRealigner -R $gatkrefseq \
					-targetIntervals $target.intervals -I $dep \
					-o $target.temp {*}$extra 2>@ stderr >@ stdout
				catch {file rename -force $target.temp.bai $target.bai}
				catch {file delete $target.intervals}
				file rename -force $target.temp $target
			}
		}
		set root r$root
		job bamrealign_index-$root -optional 1 -deps $dir/$pre-$root.bam {*}$skips -targets $dir/$pre-$root.bam.bai -code {
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
		job bamclean_clipamplicons_index-$root -optional 1 -deps $dir/$pre-$root.bam -targets $dir/$pre-$root.bam.bai -code {
			exec samtools index $dep >@ stdout 2>@ stderr
			puts "making $target"
		}
	}
	if {$cleanup} {
		cleanup_job bamclean_remtemp-$root $cleanuplist [list $dir/$pre-$root.bam]
	}
	return $dir/$pre-$root.bam
}

