proc cg_realign_gatk {args} {
	set regionfile {}
	set threads 2
	set refseq {}
	set bamfile -
	set sourcefile -
	set resultfile -
	set inputformat -
	set outputformat -
	cg_options realign_gatk args {
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
		realign around indels using gatk
	}
	if {$inputformat eq "-"} {set inputformat [ext2format $sourcefile bam {bam cram sam}]}
	if {$outputformat eq "-"} {set outputformat [ext2format $resultfile bam {bam cram sam}]}
	set refseq [refseq $refseq]
	set gatkrefseq [gatk_refseq $refseq]
	set dict [file root $gatkrefseq].dict
	if {$sourcefile eq "-"} {
		set tempfile [tempfile].bam
		exec samtools view -b -u -o $tempfile <@ stdin
		exec samtools index $tempfile
		set sourcefile $tempfile
	} elseif {$inputformat ne "bam"} {
		set tempfile [tempfile].bam
		exec samtools view -b -u $sourcefile -o $tempfile
		analysisinfo_write $sourcefile $tempfile
		set sourcefile $tempfile
	}
	if {$resultfile eq "-"} {
		set tempresult [tempfile]
	} else {
		set tempresult [filetemp $resultfile 0 1]
	}
	if {$regionfile eq ""} {
		set regionfile [cg_bam2reg -mincoverage 3 $sourcefile]
	}
	putslog "making $resultfile"
	analysisinfo_write $sourcefile $resultfile realign gatk realign_version [version gatk3]
	if {![file exists $sourcefile.[indexext $sourcefile]]} {exec samtools index $sourcefile}
	if {$regionfile ne ""} {
		set bedfile [tempbed $regionfile $refseq]
		lappend realignopts -L $bedfile
	}
	gatk3exec {-XX:ParallelGCThreads=1 -Xms512m -Xmx8g} RealignerTargetCreator \
		-R $gatkrefseq -I $sourcefile -o $tempresult.intervals {*}$realignopts
	if {[loc_compare [version gatk3] 2.7] >= 0} {
		set extra {--filter_bases_not_stored}
	} else {
		set extra {}
	}
	lappend extra --filter_mismatching_base_and_quals
	if {$resultfile eq "-"} {
		set compressionlevel 0
	} else {
		set compressionlevel [defcompressionlevel 5]
	}
	gatk3exec {-XX:ParallelGCThreads=1 -Xms512m -Xmx8g} IndelRealigner -R $gatkrefseq \
		-targetIntervals $tempresult.intervals -I $sourcefile --bam_compression $compressionlevel \
		-o $tempresult {*}$extra
	catch {file delete -- $tempresult.intervals}
	if {$outputformat eq "cram"} {
		file rename $tempresult $tempresult.bam
		exec samtools view -C $tempresult.bam -T $::env(REFSEQ) -o $tempresult
		file delete $tempresult.bam
		file delete $tempfile
		if {$resultfile eq "-"} {
			file2stdout $tempresult
		} else {
			catch {file rename -force -- $tempresult.[indexext $tempresult] $tempresult.[indexext $tempresult]}
			file rename -force -- $tempresult $resultfile
		}
	} elseif {$resultfile eq "-"} {
		file2stdout $tempresult
	} else {
		catch {file rename -force -- $tempresult.[indexext $tempresult] $tempresult.[indexext $tempresult]}
		file rename -force -- $tempresult $resultfile
	}
}
