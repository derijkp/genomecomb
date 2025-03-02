proc validate_meth_remora {refseq preset distrreg} {
	if {[catch {exec which modkit}]} {
		error "command \"modkit\" not available (needed for remora analysis), make sure it is installed, e.g. using \"cg install\""
	}
	# could check whether remora/methylation data is in ubams, but ubams are not passed to validation (yet)
	return
}

proc modkit_job {args} {
	set refseq {}
	set methodused remora
	cg_options modkit args {
		-methodused {set methodused $value}
		-refseq {set refseq $value}
	} bamfile 1 1
	set bamfile [file_absolute $bamfile]
	set destdir [file dir $bamfile]
	set refseq [refseq $refseq]
	set rootname [file_rootname $bamfile]
	set target $destdir/meth-$methodused-$rootname.bed.gz
	if {[file extension $bamfile] eq ".cram"} {
		# modkit does not work well with cram (# "detected non-BAM input format, please consider using BAM, CRAM may be unstable")
		set tempfile [tempfile].bam
		catch_exec samtools view -b -1 -T $refseq $bamfile > $tempfile
		catch_exec samtools index $tempfile
		set bamfile $tempfile
	}
	exec modkit pileup $bamfile $target.temp.gz \
		--ref $refseq \
		--preset traditional \
		>@ stdout 2>@ stderr
	file rename $target.temp.gz $target
	if {[file extension $bamfile] eq ".cram"} {
		file delete $tempfile
	}
}

proc meth_remora_job {args} {
	# putslog [list meth_nanopolish_job {*}$args]
	global appdir
	upvar job_logdir job_logdir
	set cmdline [clean_cmdline cg meth_nanopolish {*}$args]
	set refseq {}
	set skips {}
	set resultfile {}
	set basecaller {}
	set threads 1
	set distrreg 5000000
	set opts {}
	set callerroot remora
	cg_options meth_nanopolish args {
		-preset {
			# not supported
		}
		-refseq {
			set refseq $value
		}
		-skip {
			lappend skips -skip $value
		}
		-threads {
			set threads $value
		}
		-distrreg {
			# not supported
			set distrreg $value
		}
		-opts {
			set opts $value
		}
	} {fast5dir fastqdir bamfile resultfile} 3 4
	set bamfile [file_absolute $bamfile]
	set bamtail [file tail $bamfile]
	# fast5dir is not used for this
	# set fast5dir [file_absolute $fast5dir]
	# fastqdir is not used for this
	# set fastqdir [file_absolute $fastqdir]
	set refseq [refseq $refseq]
	if {$resultfile eq ""} {
		set root $callerroot-[file_rootname $bamtail]
		set resultfile [file dir $bamfile]/meth-$root.bed.gz
	} else {
		set resultfile [file_absolute $resultfile]
	}
	set destdir [file dir $resultfile]
	job_logfile $destdir/meth_nanopolish_[file root $bamtail] $destdir $cmdline \
		{*}[versions samtools gnusort8 modkit os]
	# start
	set keeppwd [pwd]
	# set regions [list_remove [distrreg_regs $distrreg $refseq] unaligned]
	set bamfileindex [index_file $bamfile]
	bam_index_job $bamfile

	# distrmethod regions
	# -------------------
	set target $resultfile
	set target2 [file root [gzroot $resultfile]].tsv.zst
	set bamfileindex $bamfile.[indexext $bamfile]
	job [job_relfile2name meth_remora_modkit- $bamfile] {*}$skips -cores $threads -deps {
		$bamfile $bamfileindex
	} -targets {
		$target $target2
	} -vars {
		bamfile refseq
	} -code {
		file delete -force $target.temp
		if {[file extension $bamfile] eq ".cram"} {
			# modkit does not work well with cram (# "detected non-BAM input format, please consider using BAM, CRAM may be unstable")
			set tempfile [tempfile].bam
			catch_exec samtools view -b -1 -T $refseq $bamfile > $tempfile
			catch_exec samtools index $tempfile
			set bamfile $tempfile
		}
		exec modkit pileup $bamfile $target.temp \
			--ref $refseq \
			--preset traditional \
			--only-tabs \
			>@ stdout 2>@ stderr
		cg gzip $target.temp
		file rename $target.temp.gz $target
		set o [wgzopen $target2]
		puts $o [join {chromosome begin end modification_code score strand thickstart thickend color coverage modified_frequency n_mod n_canonical n_othermod n_delete n_fail n_diff n_nocall} \t]
		set f [gzopen $target]
		while {[gets $f line] != -1} {
			regsub -all " " $line \t line
			puts $o $line
		}
		gzclose $f
		gzclose $o
	}
	return [list $target2 $resultfile]
}
