proc bam_sort {args} {
	set method biobambam
	set sort coordinate
	set inputformat bam
	cg_options bam_sort args {
		-method {
			if {$value ni {biobambam samtools}} {error "bamsort: unsupported -method $value"}
			set method $value
		}
		-sort {
			if {$value ni {coordinate name hash}} {error "bamsort: unsupported -sort $value"}
			set sort $value
		}
		-inputformat {
			set inputformat $value
		}
	} {sourcefile resultfile}
	if {$method eq "biobambam"} {
		set tempresult [filetemp $resultfile]
		catchstderr_exec bamsort I=$sourcefile inputformat=$inputformat SO=$sort tmpfile=[scratchfile] index=1 indexfilename=$tempresult.bai O=$tempresult
		file rename -force $tempresult.bai $resultfile.bai
		file rename -force $tempresult $resultfile
	} else {	
		set opts {}
		if {$method eq "name"} {
			lappend opts -n
		}
		if {[catch {version samtools 1}]} {
			if {[catch {exec samtools sort {*}$opts $sourcefile $resultfile.temp 2>@ stdout} msg]} {
				error $msg
			}
			if {[file exists $resultfile.temp.bam]} {
				file rename -force $resultfile.temp.bam $resultfile
			} else {
				file rename -force $resultfile.temp $resultfile
			}
		} else {
			if {[catch {exec samtools sort {*}$opts $sourcefile > $resultfile.temp 2>@ stdout} msg]} {
				error $msg
			}
			file rename -force $resultfile.temp $resultfile
		}
	}
}
