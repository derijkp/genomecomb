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
	catchstderr_exec java -jar $picard/$cmd.jar {*}$args
}

proc findpicard {} {
	global picard
	if {![info exists picard]} {
		set picard [searchpath PICARD picard picard*]
	}
	return $picard
}

proc picardversion {{minversion {}}} {
	global picardversion
	if {![info exists picardversion]} {
		catch {picard MarkDuplicates --version} picardversion
		set picardversion [lindex [split $picardversion \n] 0]
	}
	if {$minversion ne ""} {
		if {[lindex [ssort -natural [list $minversion $picardversion]] 0] ne "$minversion"} {
			error "picard version ($picardversion) smaller than $minversion"
		}
	}
	return $picardversion
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
	if {$minversion ne ""} {
		if {[lindex [ssort -natural [list $minversion $gatkversion]] 0] ne "$minversion"} {
			error "gatk version ($gatkversion) smaller than $minversion"
		}
	}
	return $gatkversion
}

proc samtoolsversion {{minversion {}}} {
	global samtoolsversion
	if {![info exists samtoolsversion]} {
		catch {exec samtools} temp
		regexp {Version: ([^\n]+)} $temp temp samtoolsversion
	}
	if {$minversion ne ""} {
		if {[lindex [ssort -natural [list $minversion $samtoolsversion]] 0] ne "$minversion"} {
			error "samtools version ($samtoolsversion) smaller than $minversion"
		}
	}
	return $samtoolsversion
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
