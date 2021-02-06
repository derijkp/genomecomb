proc mergesamples_job {result args} {
	set cmdline "[list cd [pwd]] \; [list cg mergesamples {*}$args]"
	cg_options multidb args {
	} result 1 ... {
		merge different samples into one sample. some files, like bams, are actually merged, some, like variant calls, have to be recalled on the merged samples
	}
	set result [file_absolute $result]
	putsvars result
	file mkdir $result
	set destsample [file tail $result]
	job_logfile $result/mergesamples_$destsample $result \
		$cmdline \
		{*}[versions genomecomb samtools]
	file mkdir $result/orisamples
	file mkdir $result/fastq
	foreach src $args {
		putsvars src
		mklink $src $result/orisamples/[file tail $src]
		puts "Linking fastqs from $src"
		foreach file [glob $src/fastq/*] {
			set dest $result/fastq/[file tail $file]
			if {![file exists $dest]} {
				mklink $file $dest
			}
		}
		set fast5s [glob -nocomplain $src/fast5/*]
		if {[llength $fast5s]} {
			puts "Linking fast5s from $src"
			file mkdir $result/fast5
			foreach file $fast5s {
				set dest $result/fast5/[file tail $file]
				if {![file exists $dest]} {
					mklink $file $dest
				}
			}
		}
	}
	if {[llength $args] == 1} {
		set src [lindex $args 0]
		puts "Linking results from $src"
		foreach file [list_remove [glob $src/*] $src/fastq $src/fast5] {
			set srcsample [file tail [file dir $file]]
			set tail [file tail $file]
			regsub -- -$srcsample\\. $tail -$destsample. tail
			mklink $file $result/$tail
		}
		return
	}
	set bams [glob -nocomplain $result/orisamples/*/*.bam]
	if {[llength $bams]} {
		puts "merge bams"
		unset -nocomplain a
		foreach bam $bams {
			set type [join [lrange [split [file tail $bam] -] 0 1] -]
			lappend a($type) $bam
		}
		foreach type [array names a] {
			set bams $a($type)
			set target $result/$type-$destsample.bam
			job mergebams-$type-$destsample -deps $bams -targets {
				$target
			} -vars {
				type destsample bams result
			} -code {
				file copy [analysisinfo_file [lindex $bams 0]] [analysisinfo_file $result/$type-$destsample.bam]
				exec samtools merge -c -p --threads 8 $result/temp$type-$destsample.bam {*}$bams
				file rename -force $result/temp$type-$destsample.bam $target
			}
		}
	}
	set smethfiles [gzfiles $result/orisamples/*/smeth-*.tsv]
	if {[llength $smethfiles]} {
		puts "merge meth"
		set extension [file extension [lindex $smethfiles 0]]
		set target $result/smeth-nanopolish-$destsample.tsv$extension
		job mergesmeth-$destsample -deps $smethfiles -targets {
			$target
		} -vars {
			smethfiles destsample
		} -code {
			exec devcg mergesorted {*}$smethfiles | cg [string range $extension 1 end] > $target.temp$extension
			file rename $target.temp$extension $target
			file copy -force [analysisinfo_file [lindex $smethfiles 0]]] [analysisinfo_file $target]
			exec cg meth_nanopolish_freqs $target meth-nanopolish-$destsample$extension
		}
	}
}

proc cg_mergesamples {args} {
	set args [job_init {*}$args]
	mergesamples_job {*}$args
	job_wait
}
