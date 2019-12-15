proc cg_realign_abra {args} {
	upvar job_logdir job_logdir
	set regionfile {}
	set threads 2
	set refseq {}
	set bamfile -
	set sourcefile -
	set resultfile -
	set inputformat -
	set outputformat -
	cg_options realign_abra args {
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
		realign around indels using abra
	}
	analysisinfo_write $sourcefile $resultfile realign abra realign_version [version abra]
	if {$inputformat eq "-"} {set inputformat [ext2format $sourcefile bam {bam sam}]}
	if {$outputformat eq "-"} {set outputformat [ext2format $resultfile bam {bam sam}]}
	set refseq [refseq $refseq]
	set abra [findjar abra]
	set sourcefile [tempbam $sourcefile $inputformat $refseq]
	if {$resultfile eq "-"} {
		set tempresult [tempfile]
	} else {
		set tempresult [filetemp $resultfile 0 1]
	}
	set indexext [indexext $resultfile]
	if {$regionfile eq "" || ![file exists $regionfile]} {
		putslog "making regionfile"
		set regionfile [cg_bam2reg -mincoverage 3 $sourcefile]
	}
	putslog "making $resultfile"
	if {![file exists $sourcefile.[indexext $sourcefile]]} {exec samtools index $sourcefile}
	if {$regionfile ne ""} {
		set regionfile [tempbed $regionfile $refseq]
	}
	set workingdir [scratchdir]/work
	file mkdir $workingdir
	# puts stderr [list java -Xmx4G -XX:ParallelGCThreads=1 -jar $abra --in $sourcefile --out $tempresult \
			--ref $refseq --targets $regionfile --threads $threads --working $workingdir]
	catch_exec java -Xmx4G -XX:ParallelGCThreads=1 -jar $abra --in $sourcefile --out $tempresult \
			--ref $refseq --targets $regionfile --threads $threads --working $workingdir
	if {$resultfile eq "-"} {
		file2stdout $tempresult
	} else {
		file rename -force -- $tempresult $resultfile
	}
}
