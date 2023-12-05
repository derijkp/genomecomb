# install modekit
proc install_modkit {bindir} {
	mkdir $bindir/modkit-0.2.1-linux-x86_64
	cd $bindir/modkit-0.2.1-linux-x86_64
	wgetfile https://github.com/nanoporetech/modkit/releases/download/v0.2.1/modkit_v0.2.1_centos7_x86_64.tar.gz
	exec tar xvzf modkit_v0.2.1_centos7_x86_64.tar.gz
	file delete modkit_v0.2.1_centos7_x86_64.tar.gz
	foreach file [glob dist/*] {
		file rename $file .
	}
	file delete dist
	cd $bindir
	exec ln -s modkit-0.2.1-linux-x86_64/modkit modkit
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
	exec modkit pileup $bamfile $target.temp.gz \
		--ref $refseq \
		--preset traditional \
		>@ stdout 2>@ stderr
	file rename $target.temp.gz $target
	
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
	file mkdir $resultfile.temp
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
	job [job_relfile2name seqmeth_nanopolish_index- $bamfile] {*}$skips -cores $threads -deps $deps -targets {
		$target $target2
	} -vars {
		bamfile refseq
	} -code {
		exec modkit pileup $bamfile $target.temp \
			--ref $refseq \
			--preset traditional \
			>@ stdout 2>@ stderr
		cg gzip $target.temp
		file rename $target.temp.gz $target
		set o [wgzopen $target2]
		puts $o [join {chromsome begin end modification_code score strand thickstart thickend color coverage modified_frequency n_mod n_canonical n_othermod n_delete n_fail n_diff n_nocall} \t]
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
