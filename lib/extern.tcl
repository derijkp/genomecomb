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

proc gatkexec {args} {
	global gatk
	if {![info exists gatk(version)]} {
		set gatk(usecommand) 0
		if {![info exists ::env(GATK)] || [file isfile $::env(GATK)]} {
			if {[info exists ::env(GATK)]} {
				set gatk(command) $::env(GATK)
			} else {
				set gatk(command) gatk
			}
			catch {exec $gatk(command) HaplotypeCaller --version} msg
			if {[regexp {GATK.*Version:(.*)} $msg temp version]} {
				set gatk(usecommand) 1
				set gatk(version) $version
				set ::gatkjava java
			} elseif {[info exists ::env(GATK)]} {
				error "gatk command given in env var GATK gives following unexpected result when trying to get version:\n$msg"
			}
		} 
		if {![info exists gatk(version)]} {
			set gatk(command) 0
			set gatk(jar) [searchpath GATK gatk GenomeAnalysisTK*]/GenomeAnalysisTK.jar
			if {![catch {exec java -XX:ParallelGCThreads=1 -jar $gatk(jar) --version} msg]} {
				set version $msg
				set ::gatkjava java
			} elseif {![catch {exec java1.8 -XX:ParallelGCThreads=1 -jar $gatk(jar) --version} version]} {
				set ::gatkjava java1.8
			} elseif {![catch {exec java1.7 -XX:ParallelGCThreads=1 -jar $gatk(jar) --version} version]} {
				set ::gatkjava java1.7
			} else {
				error "Cannot determine gatk version:\n$msg"
			}
			set gatk(version) $version
		}
	}
	if {[lindex $args 0] eq "version"} {return $gatk(version)}
	set javaopts [list_pop args 0]
	if ($gatk(usecommand)) {
		catch_exec $gatk(command) --java-options $javaopts {*}$args
	} else {
		catch_exec $::gatkjava {*}$javaopts -jar $gatk(jar) -T {*}$args
	}
	return ""
}

proc gatk3exec {args} {
	global gatk3
	if {![info exists gatk3(version)]} {
		if {[catch {
			set gatk3(jar) [searchpath GATK3 gatk3 GenomeAnalysisTK*]/GenomeAnalysisTK.jar
		} msg]} {
			if {[catch {
				set gatk3(jar) [searchpath GATK gatk GenomeAnalysisTK*]/GenomeAnalysisTK.jar
			}]} {
				error $msg
			}
		}
		
		if {![catch {exec java -XX:ParallelGCThreads=1 -jar $gatk3(jar) --version} msg]} {
			set version $msg
			set ::gatkjava java
		} elseif {![catch {exec java1.8 -XX:ParallelGCThreads=1 -jar $gatk3(jar) --version} version]} {
			set ::gatkjava java1.8
		} elseif {![catch {exec java1.7 -XX:ParallelGCThreads=1 -jar $gatk3(jar) --version} version]} {
			set ::gatkjava java1.7
		} else {
			error "Cannot determine gatk3 version:\n$msg"
		}
		set gatk3(version) $version
	}
	if {[lindex $args 0] eq "version"} {return $gatk3(version)}
	set javaopts [list_pop args 0]
	catch_exec $::gatkjava {*}$javaopts -jar $gatk3(jar) -T {*}$args
}

#proc gatk {} {
#	global gatk
#	if {![info exists gatk]} {
#		set gatk [searchpath GATK gatk GenomeAnalysisTK*]/GenomeAnalysisTK.jar
#	}
#	return $gatk
#}

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
