proc gatkjava {} {
	if {![info exists ::gatkjava]} {
		version gatk
	}
	return $::gatkjava
}

proc version {item {minversion {}}} {
	global _versions
	if {![info exists _versions($item)]} {
		set _versions($item) ?
		switch $item {
			genomecomb {
				set _versions($item) $::genomecomb::version.$::genomecomb::patchlevel
			}
			fileformat {
				set _versions($item) 0.99
			}
			dbdir {
				if {![catch {dbdir} dbdir] && [gzfile $dbdir/README*.txt] ne ""} {
					set readme [gzfile $dbdir/README*.txt]
					if {[file exists $readme]} {
						set temp [file_read $readme]
						regexp {version: ([^\n]+)} $temp temp _versions($item)
					} else {
						set _versions($item) ?
					}
				}
			}
			samtools {
				catch {exec samtools} temp
				regexp {Version: ([^\n]+)} $temp temp _versions($item)
			}
			gatk {
				set _versions($item) [gatkexec version]
			}
			gatk3 {
				set _versions($item) [gatk3exec version]
			}
			picard {
				catch {picard MarkDuplicates --version} version_picard
				set _versions($item) [lindex [split $version_picard \n] 0]
				regexp {Version:([^\n]+)} $_versions($item) temp _versions($item)
			}
			fastqc {
				catch {exec fastqc -v} temp
				regsub {^[^0-9]*} $temp {} _versions($item)
			}
			java {
				catch {exec java -XX:ParallelGCThreads=1 -version} temp
				regsub {java version "([^"]+)"} $temp {\1 } temp
				set _versions($item) [join [split $temp \n] {, }]
			}
			gatkjava {
				catch {exec [gatkjava] -XX:ParallelGCThreads=1 -version} temp
				regsub {java version "([^"]+)"} $temp {\1 } temp
				set _versions($item) [join [split $temp \n] {, }]
			}
			tcl {
				set _versions($item) $::tcl_patchLevel
			}
			R {
				catch {exec [findR] --version} temp
				set _versions($item) [lindex [split $temp \n] 0]
			}
			fastq-mcf {
				set temp [file_read $::externdir/ea-utils/README]
				regexp {\(version ([^)]+)\)} $temp temp temp
				set _versions($item) "$temp adapted"
			}
			fastq-stats {
				set temp [file_read $::externdir/ea-utils/README]
				regexp {\(version ([^)]+)\)} $temp temp temp
				set _versions($item) "$temp adapted"
			}
			plink {
				catch {exec plink --version --noweb} temp
				regsub {^[^v]*v} $temp {\1 } temp
				set temp [string trim [lindex [split $temp \n] 0]]
				regsub -all {[ |]+} $temp { } temp
				set _versions($item) $temp
			}
			primer3 {
				catch {exec primer3_core -version} temp
				regexp {release ([^)]+)} $temp temp _versions($item)
			}
			lz4 {
				catch {exec lz4 -v} temp
				regexp { v([0-9.]+)} $temp temp _versions($item)
			}
			zst - zstd {
				set _versions($item) [zstversion]
			}
			dot {
				catch {exec dot -V} temp
				regsub {^[^0-9]*} $temp {} temp
				set _versions($item) $temp
			}
			os {
				set os ? ; set release ? ; set kernel_version ? ; set machine ?
				catch {exec uname -s} os
				catch {exec uname -r} release
				catch {exec uname -v} kernel_version
				catch {exec uname -m} machine
				set _versions($item) "$release $os $machine $kernel_version"
			}
			biobambam {
				catch {biobambam bammarkduplicates2 -h} temp
				regexp {version ([0-9.]+)} $temp temp temp
				set _versions($item) [string trimright $temp .]
			}
			albacore {
				catch {exec read_fast5_basecaller.py -v} temp
				regexp {version ([0-9.]+)} $temp temp temp
				set _versions($item) [string trimright $temp .]
			}
			freebayes {
				catch {exec freebayes --version} temp
				regsub {version: *} $temp {} temp
				set _versions($item) [string trim $temp]
			}
			razip {
				catch {exec grep Version [exec which razip]_docs/razf.c} temp
				regsub {.*Version: *} $temp {} temp
				set _versions($item) [string trim $temp]
			}
			bgz {set _versions($item) [version bgzip]}
			gz {set _versions($item) [version gzip]}
			default {
				if {[info command version_$item] ne "" || [auto_load version_$item]} {
					set _versions($item) [version_$item]
				} else {
					if {![catch {exec $item --version} temp]} {
					} elseif {![catch {exec $item -version} temp]} {
					} elseif {![catch {exec $item -h} temp]} {
					} elseif {![catch {exec $item -V} temp]} {
					} elseif {![catch {exec $item -v} temp]} {
					} else {
						catch {exec $item} temp
					}
					set line1 [lindex [split $temp \n] 0]
					if {[regexp dir= $line1]} {set line1 [lindex [split $temp \n] 1]}
					if {[regexp {^couldn't execute} $line1 temp]} {
						set _versions($item) ?
					} elseif {[regexp {[0-9.]+-?[abr][0-9]+$} $line1 temp]} {
						set _versions($item) $temp
					} elseif {[regexp {[0-9.]+[0-9]$} $line1 temp]} {
						set _versions($item) $temp
					} elseif {[regexp {([0-9.]+\.[0-9]*)} $line1 temp]} {
						set _versions($item) $temp
					} elseif {[regexp {([0-9]+[0-9.a-zA-Z]*)} $line1 temp temp]} {
						set _versions($item) $temp
					} elseif {[regsub {^[^0-9\n]+} $temp {} temp]} {
						set temp [lindex [split [string trim $temp] \n] 0]
						set _versions($item) $temp
					} else {
						regsub {^.*[Vv]ersion:? } $temp {} temp
						set temp [lindex [split [string trim $temp] \n] 0]
						set _versions($item) $temp
					}
				}
			}
		}
		if {$_versions($item) eq ""} {
			set _versions($item) ?
		} else  {
			set _versions($item) [string_change $_versions($item) [list \t \\t \n \\n]]
		}
	}
	if {$minversion ne ""} {
		if {[lindex [bsort [list $minversion $_versions($item)]] 0] ne "$minversion"} {
			error "$item version ($_versions($item)) smaller than $minversion"
		}
	}
	return $_versions($item)
}

proc minversion {version minversion} {
	if {[lindex [bsort [list $minversion $version]] 0] ne "$minversion"} {
		return 0
	} else {
		return 1
	}
}

proc versions {args} {
	if {![llength $args]} {
		set args {genomecomb dbdir fastqc fastq-stats fastq-mcf bwa bowtie2 samtools gatk biobambam picard plink primer3 java tcl R gnusort8 tabix zst os}
	}
	set result {}
	foreach item $args {
		lappend result $item [version $item]
	}
	return $result
}

proc cg_version {args} {
	set item genomecomb
	cg_options versions args {
	} {item} 0 2 {
		returns the version of the given item/program
	}
	puts [version $item {*}$args]
}

proc cg_versions {args} {
	# extra commands used, but not in def list: {wget gzip gunzip zcat cat paste tail wc cp ln bigWigToBedGraph find grep chmod tar}
	# give no version: bgzip bzcat razip
	cg_options versions args {
	} {} 0 ... {
		returns the (current) versions of the given programs as a tsv file
	}
	if {![llength $args]} {
		set args {genomecomb dbdir fastqc fastq-stats fastq-mcf bwa bowtie2 samtools gatk biobambam picard plink primer3 java R gnusort8 tabix zst os}
	}
	puts "item\tversion"
	foreach item $args {
		puts $item\t[version $item]
	}
}
