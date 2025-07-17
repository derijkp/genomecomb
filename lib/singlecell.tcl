proc sc_pre_job {args} {
	upvar job_logdir job_logdir
	cg_options sc_pre args {
	} {fastqdir sampledir} 1 2 {
		creates summary files needed for further sc processing (if not there already)
	}
	set fastqfiles [gzfiles $fastqdir/*.fq $fastqdir/*.fastq $fastqdir/*.bam $fastqdir/*.cram $fastqdir/*.sam]
	set sample [file tail $sampledir]
	job singlecell_pre-$sample \
	-deps $fastqfiles \
	-targets {
		$sampledir/reads_per_cell_raw.tsv
		$sampledir/umis_per_cell_raw.tsv
	} -vars {
		sampledir fastqfiles sample
	} -code {
		set o [wgzopen $sampledir/reads_per_cell_raw.tsv.temp.zst]
		puts $o cellbarcode\tumi
		foreach fastq $fastqfiles {
			if {[file ext $fastq] in ".bam .cram .sam"} {
				set usefastq [tempfile].fastq.gz
				catch_exec samtools fastq -T "RG,CB,QT,MI,MM,ML,Mm,Ml" $fastq | gzip > $usefastq
				set ubams 1
			} else {
				set usefastq $fastq
				set ubams 0
			}
			set f [gzopen $usefastq]
			while 1 {
				if {[gets $f read] == -1} break
				if {[regexp {([A-Z]+)_([A-Z]+)#} $read temp cellbarcode umi]} {
					puts $o $cellbarcode\t$umi
				}
				if {[gets $f line] == -1} break
				if {[gets $f line] == -1} break
				if {[gets $f line] == -1} break
			}
			gzclose $f
		}
		gzclose $o
		# read_counts
		cg select -overwrite 1 -g cellbarcode -s -count \
			$sampledir/reads_per_cell_raw.tsv.temp.zst $sampledir/reads_per_cell_raw.tsv.temp2
		file rename -force $sampledir/reads_per_cell_raw.tsv.temp2 $sampledir/reads_per_cell_raw.tsv
		# umi_counts
		cg select -overwrite 1 -optim memory -g {cellbarcode * umi *} $sampledir/reads_per_cell_raw.tsv.temp.zst \
			| cg select -g cellbarcode -s -count \
			> $sampledir/umi_per_cell_raw.tsv.temp2
		file delete $sampledir/reads_per_cell_raw.tsv.temp.zst
		file rename -force $sampledir/umi_per_cell_raw.tsv.temp2 $sampledir/umis_per_cell_raw.tsv
	}
}