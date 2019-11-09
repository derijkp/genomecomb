proc cg_bcl2fastq {rundir outdir {rtr 6} {dtr 6} {ptr 6} {wtr 6} } {
	#-r, --loading-threads Number of threads used for loading BCL data.
	#-d, --demultiplexing-threads Number of threads used for demultiplexing.
	#-p, --processing-threads  Number of threads used for processing demultiplexed data.
	#-w, --writing-threads Number of threads used for writing FASTQ data. This must not be higher than number of samples.
	exec [bcl2fastq] --create-fastq-for-index-reads -r $rtr -d $dtr -p $ptr -w $wtr --runfolder-dir $rundir --output-dir $outdir 2>@ stderr >@ stdout
}

proc cg_process_conv_illnextseq {illsrc destdir} {
	set illsrc [file_absolute $illsrc]
	set destdir [file_absolute $destdir]
	file mkdir $destdir
	set keeppwd [pwd]
	cd $destdir
	# copy files from illumina dir
	set files [glob $illsrc/*_R*.fastq*]
	foreach file $files {
		set sample [file tail $file]
		regsub {_[^_]+_[^_]+_[^_]+_[^_]+\.fastq.*} $sample {} sample
		regsub -all -- - $sample _ sample 
		if {$sample != "Undetermined"} {
			file mkdir $destdir/$sample/ori
			file mkdir $destdir/$sample/ori/fastq
			file mkdir $destdir/$sample
			file mkdir $destdir/$sample/fastq
			hardlink $file $destdir/$sample/ori/fastq
			cplinked $destdir/$sample/ori/fastq $destdir/$sample/fastq
		}
	}
	cd $keeppwd
}

proc cg_process_illumina {args} {
	error "This command has been depricated, use cg process_project instead"
}
