proc searchpath {envvar args} {
	set name [lindex $args 0]
	if {$envvar ne "" && [info exists ::env($envvar)]} {
		if {![file exists $::env($envvar)]} {
			error "$name not found at $::env($envvar) (from env var $envvar)"
		}
		return $::env($envvar)
	} else {
		set dirlist [split [get ::env(PATH) ""] [pathsep]]
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
	catchstderr_exec java -XX:ParallelGCThreads=1 -jar $picard/$cmd.jar {*}$args
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

proc findjar {program {envvar {}}} {
	global extprograms
	if {![info exists extprograms($program)]} {
		set extprograms($program) [searchpath $envvar $program/$program.jar $program/$program*.jar $program/*.jar $program.jar $program*.jar]
	}
	return $extprograms($program)
}

proc findR {} {
	global R
	if {![info exists R]} {
		set R [searchpath RCG Rcg R]
	}
	return $R
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
