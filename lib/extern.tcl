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
					return [lindex [bsort $dirs] end]
				}
			}
		}
		error "$name not found in PATH"
	}
}

proc picard {cmd args} {
	set maxmem 4G
	set picard [findpicard]
	if {[file exists $picard/$cmd.jar]} {
		catch_exec java -Xms1G -Xmx${maxmem} -Djava.io.tmpdir=[scratchdir]] -XX:ParallelGCThreads=1 -jar $picard/$cmd.jar {*}$args
	} else {
		catch_exec java -Xms1G -Xmx${maxmem} -Djava.io.tmpdir=[scratchdir]] -XX:ParallelGCThreads=1 -jar $picard/picard.jar $cmd {*}$args
	}
}

proc findpicard {} {
	global picard
	if {![info exists picard]} {
		set picard [searchpath PICARD picard2 picard*]
		if {$picard eq ""} {
			set picard [searchpath PICARD picard picard*]
		}
	}
	return $picard
}

proc findjava {jar {version {}}} {
	if {![info exists javavcmd($jar)]} {
		if {![catch {exec java$version -XX:ParallelGCThreads=1 -jar $jar} msg]
			|| (![regexp {Unsupported major.minor version} $msg] && ![regexp {couldn't execute} $msg])
		} {
			set javavcmd($jar) java$version
		} elseif {![catch {exec java -XX:ParallelGCThreads=1 -jar $jar} msg]
			|| (![regexp {Unsupported major.minor version} $msg] && ![regexp {couldn't execute} $msg])
		} {
			set javavcmd($jar) java
		} elseif {![catch {exec java1.8 -XX:ParallelGCThreads=1 -jar $jar} msg2]
			|| (![regexp {Unsupported major.minor version} $msg2] && ![regexp {couldn't execute} $msg2])
		} {
			set javavcmd($jar) java1.8
		} elseif {![catch {exec java1.7 -XX:ParallelGCThreads=1 -jar $jar} msg2]
			|| (![regexp {Unsupported major.minor version} $msg2] && ![regexp {couldn't execute} $msg2])
		} {
			set javavcmd($jar) java1.7
		} else {
			error "Cannot determine java version for $jar:\n$msg"
		}
	}
	return $javavcmd($jar)
}

proc execjar {args} {
	set java java
	set mem 1G
	set javaopts ""
	set finishedpattern {}
	set pos 0
	set javaversion {}
	foreach {key value} $args {
		if {$key eq "-mem"} {
			set mem $value
			incr pos 2
		} elseif {$key eq "-javaopts"} {
			set javaopts $value
			incr pos 2
		} elseif {$key eq "-javaversion"} {
			set javaversion $value
			incr pos 2
		} elseif {$key eq "-finishedpattern"} {
			set mem $value
			incr pos 2
		} else {
			break
		}
	}
	lappend javaopts -Xms$mem -Xmx$mem -XX:ParallelGCThreads=1
	set args [lrange $args $pos end]
	set cmd [list_pop args 0]
	set jar [findjar $cmd [string toupper $cmd]]
	set java [findjava $jar $javaversion]
	set error [catch {
		exec $java {*}$javaopts -jar $jar {*}$args
	} msg opt]
	if {$error} {
		if {$::errorCode ne "NONE"} {
			if {$finishedpattern eq "" || ![regexp $finishedpattern $msg]} {
				dict unset opt -level
				return -options $opt $msg
			}
		}
	}
	return ""
}

proc gatkexec {args} {
	if {[lindex $args 0] eq "-finishedpattern"} {
		# hack: If the given (regexp) pattern is present in the output, ignore errors
		# This hack was added because haplotype caller sometimes fails with an error after
		# completely processing the data and producing a compete gvcf (probably in cleanup)
		# Using this these errors are ignored (if the analysis complete pattern is encountered)
		# and there is no "ERROR" or "EXCEPTION" in the output
		set finishedpattern [lindex $args 1]
		set args [lrange $args 2 end]
	} else {
		set finishedpattern {}
	}
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
			if {
				[regexp {The Genome Analysis Toolkit \(GATK\) v([^\n]+)} $msg temp version]
				|| [regexp {GATK.*Version:(.+)} $msg temp version]
				|| [regexp {Version:([^\n]+).*GATK} $msg temp version]
			} {
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
		set error [catch {
			exec $gatk(command) --java-options $javaopts {*}$args
		} msg opt]
	} else {
		set error [catch {
			exec $::gatkjava {*}$javaopts -jar $gatk(jar) -T {*}$args
		} msg opt]
	}
	if {$error} {
		if {$::errorCode ne "NONE"} {
			if {$finishedpattern eq "" || ![regexp $finishedpattern $msg] || [regexp ERROR $msg] || [regexp EXCEPTION $msg]} {
				dict unset opt -level
				return -options $opt $msg
			}
		}
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
		
		if {![catch {exec java1.8 -XX:ParallelGCThreads=1 -jar $gatk3(jar) --version} version]} {
			set ::gatkjava java1.8
		} elseif {![catch {exec java1.7 -XX:ParallelGCThreads=1 -jar $gatk3(jar) --version} version]} {
			set ::gatkjava java1.7
		} elseif {![catch {exec java -XX:ParallelGCThreads=1 -jar $gatk3(jar) --version} msg]} {
			set version $msg
			set ::gatkjava java
		} else {
			error "Cannot determine gatk3 version:\n$msg"
		}
		set gatk3(version) $version
	}
	if {[lindex $args 0] eq "version"} {return $gatk3(version)}
	set javaopts [list_pop args 0]
	catch_exec $::gatkjava {*}$javaopts -jar $gatk3(jar) -T {*}$args
}

proc gatk_IndexFeatureFile_job {file {name {}}} {
	upvar job_logdir job_logdir
	if {[file extension $file] eq ".gz"} {set target $file.tbi} else {set target $file.idx}
	job gatkindex_$name[file tail $file] -deps {
		$file
	} -targets {
		$target
	} -vars {
		file
	} -code {
		cg_gatk_index $file
	}
	return $target
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
		set R [searchpath RCG dirR-4.2.1/R dirR-4.1.2/R dirR Rcg R]
	}
	return $R
}

proc R {args} {
	set vars {}
	set listvars {}
	set outRfile {}
	cg_options cat args {
		-vars {
			set vars $value
		}
		-listvars {
			set listvars $value
		}
		-outRfile {
			set outRfile $value
		}
	} cmd 1 1
	set pre "setwd(\"[pwd]\")\n"
	foreach var $vars {
		upvar 1 $var uvar
		if {![info exists uvar]} {
			error "error running R with -vars: variable $var does not exist"
		}
		if {$uvar ne "" && [string is double $uvar]} {
			append pre "$var=$uvar\n"
		} else {
			append pre "$var=\"$uvar\"\n"
		}
	}
	foreach var $listvars {
		upvar 1 $var uvar
		if {![info exists uvar]} {
			error "error running R with -vars: variable $var does not exist"
		}
		set type num
		foreach el $uvar {
			if {![string is double $el]} {set type string ; break}
		}
		if {$type eq "num"} {
			append pre "$var=c([join $uvar ,])\n"
		} else {
			append pre "$var=c(\"[join $uvar \",\"]\")\n"
		}
	}
	if {$outRfile eq ""} {
		set outRfile [tempfile].R
	}
	file_write $outRfile $pre$cmd
	puts "running R from $outRfile"
	exec [findR] --vanilla < $outRfile >@ stdout 2>@ stderr
	# ps: you can run R from script with parameters like this (but easier to include in runfile):
	# dirR --vanilla --slave --no-restore --file=rfile --args arg1 arg2
}

proc findpython3 {} {
	global python3
	if {![info exists python3]} {
		if {![catch {exec which python3} temp]} {
			set python3 $temp
		} elseif {![catch {exec which python} temp]} {
			set python3 $temp
		} else {
			error "python3 not found (not even python)"
		}
	}
	return $python3
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
