proc searchpath {envvar args} {
	set name [lindex $args 0]
	if {[info exists ::env($envvar)]} {
		if {![file exists $::env($envvar)]} {
			error "$name not found at $::env($envvar) (from env var $envvar)"
		}
		return $::env($envvar)
	} else {
		set dirlist [split [get ::env(PATH) ""] :]
		list_addnew dirlist $::externdir
		foreach pattern $args {
			foreach dir $dirlist {
				if {![catch {glob $dir/$pattern} dirs]} {
					return [lindex [lsort -dict $dirs] 0]
				}
			}
		}
		error "$name not found in PATH"
	}
}

proc picard {cmd args} {
	set picard [findpicard]
	if {[catch {
		exec java -jar $picard/$cmd.jar {*}$args
	} msg] && ![regexp "done. Elapsed time:" $msg]} {
		error $msg
	}
	return $msg
}

proc findpicard {} {
	global picard
	if {![info exists picard]} {
		set picard [searchpath PICARD picard picard*]
	}
	return $picard
}

proc gatk {} {
	global gatk
	if {![info exists gatk]} {
		set gatk [searchpath GATK gatk GenomeAnalysisTK*]/GenomeAnalysisTK.jar
	}
	return $gatk
}

proc gatkversion {{minversion {}}} {
	global gatkversion
	set gatk [gatk]
	if {![info exists gatkversion]} {
		set gatkversion [exec java -jar $gatk --version]
	}
	return $gatkversion
}

proc bcl2fastq {} {
	global bcl2fastq
	if {![info exists bcl2fastq]} {
		set path [searchpath bcl2fastq bcl2fastq*]
		if {[file isdir $path]} {
			set bcl2fastq $path/bin/bcl2fastq
		} else {
			set bcl2fastq $path
		}
	}
	return $bcl2fastq
}
