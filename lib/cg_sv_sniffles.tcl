proc sv_sniffles_job {args} {
	upvar job_logdir job_logdir
	set pre ""
	set opts {}
	set split 1
	set deps {}
	set regionfile {}
	set threads 2
	set cleanup 1
	set regmincoverage 3
	set resultfiles 0
	set rootname {}
	set skips {}
	set min_support 2
	set min_seq_size 300
	cg_options sv_sniffles args {
		-L - -deps {
			lappend deps $value
		}
		-pre {
			set pre $value
		}
		-split {
			set split $value
		}
		-threads - -t {
			set threads $value
		}
		-cleanup {
			set cleanup $value
		}
		-resultfiles {
			set resultfiles $value
		}
		-rootname {
			set rootname $value
		}
		-skip {
			lappend skips -skip $value
		}
		-maxdist {
			lappend opts -d $value
		}
		-min_support {
			set min_support $value
		}
		-min_seq_size {
			set min_seq_size $value
		}
		default {
			if {[regexp {^-..} $key]} {set key -$key}
			lappend opts $key $value
		}
	} {bamfile refseq}
	set bamfile [file_absolute $bamfile]
	set refseq [file_absolute $refseq]
	set destdir [file dir $bamfile]
	set file [file tail $bamfile]
	if {$rootname eq ""} {
		set root sniffles-[file_rootname $file]
	} else {
		set root $rootname
	}
	# resultfiles
	set varfile ${pre}var-$root.tsv
	set sregfile ${pre}sreg-$root.tsv
	set varallfile ${pre}varall-$root.tsv
	set resultlist [list $destdir/$varfile.lz4]
	if {$resultfiles} {
		return $resultlist
	}
	# logfile
	set cmdline [list cg sv_sniffles]
	foreach option {
		split deps pre maxdist
	} {
		if {[info exists $option]} {
			lappend cmdline -$option [get $option]
		}
	}
	lappend cmdline {*}$opts $bamfile $refseq
	job_logfile $destdir/sv_sniffles_[file tail $bamfile] $destdir $cmdline \
		{*}[versions sniffles gnusort8 lz4 os]
	# start
	## Produce sniffles SNP calls
	set keeppwd [pwd]
	cd $destdir
	set deps [list $file $refseq $file.bai {*}$deps]
	job ${pre}var-$root {*}$skips -mem 1G -cores $threads \
	-deps $deps \
	-targets {${pre}varall-$root.vcf ${pre}varall-$root.vcf.analysisinfo} \
	-skip [list $varallfile $varallfile.analysisinfo] \
	-vars {sniffles opts regionfile refseq threads root} \
	-code {
		analysisinfo_write $dep $target sample $root varcaller sniffles varcaller_version [version sniffles] varcaller_cg_version [version genomecomb]
		exec sniffles {*}$opts -threads $threads --genotype \
			--min_support $min_support --min_seq_size $min_seq_size
			-m $dep -v $target.temp 2>@ stderr >@ stdout
		file rename -force $target.temp $target
	}
	job ${pre}var-sniffles2tsv-$root {*}$skips -deps [list ${pre}var-$root.vcf] \
	-targets {$varfile.lz4 $varfile.analysisinfo} \
	-vars {sample split} \
	-code {
		analysisinfo_write $dep $target
		cg vcf2tsv -split $split -removefields {name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR} $dep $target.temp.lz4
		file rename -force $target.temp.lz4 $target
	}
	# lz4_job $varfile -i 1
	lz4index_job {*}$skips $varfile.lz4
#	# make sreg
#	sreg_sniffles_job ${pre}sreg-$root $varallfile $sregfile.lz4 $skips
	# cleanup
	cd $keeppwd
	return $resultlist
	# return [file join $destdir $varfile]
}

proc cg_sv_sniffles {args} {
	set args [job_init {*}$args]
	sv_sniffles_job {*}$args
	job_wait
}
