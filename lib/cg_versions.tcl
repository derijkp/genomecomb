proc version {item {minversion {}}} {
	global _versions
	if {![info exists _versions($item)]} {
		set _versions($item) ?
		switch $item {
			genomecomb {
				set _versions($item) $::genomecomb::version.$::genomecomb::patchlevel
			}
			dbdir {
				if {![catch {dbdir} dbdir] && [file exists $dbdir/README.txt]} {
					set temp [file_read $dbdir/README.txt]
					regexp {version: ([^\n]+)} $temp temp _versions($item)
				}
			}
			samtools {
				catch {exec samtools} temp
				regexp {Version: ([^\n]+)} $temp temp _versions($item)
			}
			gatk {
				set gatk [gatk]
				set _versions($item) [exec java -jar $gatk --version]
			}
			picard {
				catch {picard MarkDuplicates --version} version_picard
				set _versions($item) [lindex [split $version_picard \n] 0]
			}
			fastqc {
				set temp [exec fastqc -v]
				regsub {^[^0-9]*} $temp {} _versions($item)
			}
			java {
				catch {exec java -version} temp
				regsub {java version "([^"]+)"} $temp {\1 } temp
				set _versions($item) [join [split $temp \n] {, }]
			}
			fastq-mcf {
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
				catch {exec lz4c -v} temp
				regexp { v([0-9.]+)} $temp temp _versions($item)
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
			default {
				if {![catch {exec $item --version} temp]} {
				} elseif {![catch {exec $item -version} temp]} {
				} elseif {![catch {exec $item -h} temp]} {
				} elseif {![catch {exec $item -V} temp]} {
				} elseif {![catch {exec $item -v} temp]} {
				} else {
					catch {exec $item} temp
				}
				set line1 [lindex [split $temp \n] 0]
				if {[regexp {[0-9.]+[0-9]$} $line1 temp]} {
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
		if {$_versions($item) eq ""} {set _versions($item) ?}
	}
	if {$minversion ne ""} {
		if {[lindex [ssort -natural [list $minversion $_versions($item)]] 0] ne "$minversion"} {
			error "$item version ($_versions($item)) smaller than $minversion"
		}
	}
	return $_versions($item)
}

proc cg_version {args} {
	set item genomecomb
	cg_options versions args {
	} {item} 0 2
	puts [version $item {*}$args]
}

proc cg_versions {args} {
	# extra commands used, but not in def list: {wget gzip gunzip zcat cat paste tail wc cp ln bigWigToBedGraph find grep chmod tar}
	# give no version: bgzip bzcat razip
	cg_options versions args {
	} {} 0
	if {![llength $args]} {
		set args {genomecomb dbdir fastq-mcf bwa bowtie2 samtools gatk picard fastqc plink primer3 java R gnusort8 tabix lz4 os}
	}
	puts "item\tversion"
	foreach item $args {
		puts $item\t[version $item]
	}
}
