proc cg_realign_srma {args} {
	set regionfile {}
	set threads 1
	set refseq {}
	set bamfile -
	set sourcefile -
	set resultfile -
	set inputformat -
	set outputformat -
	cg_options realign_srma args {
		-regionfile {
			set regionfile $value
		}
		-refseq {
			set refseq $value
		}
		-inputformat - -if {
			set inputformat $value
		}
		-outputformat - -of {
			set outputformat $value
		}
		-threads - -t {
			set threads $value
		}
	} {sourcefile resultfile refseq} 0 3 {
		realign around indels using srma
	}
	if {$inputformat eq "-"} {set inputformat [ext2format $sourcefile bam {bam sam}]}
	if {$outputformat eq "-"} {set outputformat [ext2format $resultfile bam {bam sam}]}
	set refseq [refseq $refseq]
	set srma [findjar srma]
	set opts {}
	set optsio {}
	if {$regionfile eq "" || ![file exists $regionfile]} {
		if {$sourcefile eq "-"} {
			set tempfile [tempfile].bam
			exec samtools view --threads $threads --no-PG -b -u -o $tempfile <@ stdin
			exec samtools index $tempfile
			set sourcefile $tempfile
		}
		putslog "making regionfile"
		set regionfile [cg_bam2reg -mincoverage 3 $sourcefile]
		lappend opts I=$sourcefile
	} elseif {$sourcefile eq "-"} {
		lappend optsio <@ stdin
	} else {
		lappend opts I=$sourcefile
	}
	if {$resultfile eq "-"} {
		lappend optsio >@ stdout
		lappend opts QUIET_STDERR=true QUIET=true
	} else {
		set tempresult [filetemp $resultfile 0 1]
		lappend opts O=$tempresult
	}
	lappend opts COMPRESSION_LEVEL=[defcompressionlevel 5]
	set gatkrefseq [gatk_refseq $refseq]
	putslog "making $resultfile"
	analysisinfo_write $sourcefile $resultfile realign srma realign_version [version srma]
	if {![file exists $sourcefile.[indexext $sourcefile]]} {exec samtools index $sourcefile}
	if {$regionfile ne ""} {
		set regionfile [tempbed $regionfile $refseq]
	}
	set workingdir [scratchdir]/work
	file mkdir $workingdir
	# puts stderr [list ---- java -XX:ParallelGCThreads=1 -jar $srma R=$gatkrefseq NUM_THREADS=$threads {*}$opts {*}$optsio]
	catch_exec java -XX:ParallelGCThreads=1 -jar $srma R=$gatkrefseq NUM_THREADS=$threads {*}$opts {*}$optsio
	if {$resultfile ne "-"} {
		file rename -force -- $tempresult $resultfile
	}
}
