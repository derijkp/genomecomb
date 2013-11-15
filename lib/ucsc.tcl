#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc ucsc2region {ucsc_file} {
	catch {close $f}
	set f [open $ucsc_file]
	set temp [string range [gets $f] 1 end]
	set header [split $temp \t]
	set poss [list_cor $header {chrom chromStart chromEnd name}]
	puts [join [list chromosome start end name] \t]
	while {![eof $f]} {
		set line [getnotempty $f]
		foreach {chrom chromStart chromEnd name} [list_sub $line $poss] break
		puts [join [list $chrom $chromStart $chromEnd $name] \t]
	}
	close $f
}

proc cg_ucsc2region {args} {
	global scriptname action
	if {[llength $args] != 1} {
		puts stderr "format is: $scriptname $action ucsc_file"
		puts stderr "convert ucsc tab file to regionfile for further use in filtering"
		exit 1
	}
	foreach {ucsc_file} $args break
	ucsc2region $ucsc_file
}

proc ucsc_wibfile {file dir} {
	if {[file exists $file]} {
		return $file
	}
	set testfile $dir/[file tail $file]
	if {[file exists $testfile]} {
		return $testfile
	}
	if {[file exists [file tail $file]]} {
		return [file tail $file]
	}
	puts "Downloading $file"
	wgetfile ftp://hgdownload.cse.ucsc.edu/$file $testfile
	return $testfile
}

proc ucscwiggle2reg {ucsc_file resultfile {precision 1} {formula {}} {addnum {}}} {
	if {$formula eq ""} {
		proc formula {value} {return $value}
	} else {
		proc formula {value} "return \[expr \{$formula\}\]"
	}
	if {$addnum ne ""} {
		proc putsresult {o args} {
			puts $o [join $args \t]
		}
		set useaddnum 1
	} else {
		proc putsresult {o args} {
			puts $o [join [lrange $args 0 3] \t]
		}
		set useaddnum 0
	}
	puts "Making $resultfile"
	catch {close $f}; catch {close $b}; catch {close $o}; 
	set dir [file dir [file normalize $ucsc_file]]
	set f [gzopen $ucsc_file]
	set o [open $resultfile.temp w]
	putsresult $o chromosome begin end score num
	set num {}
	set ucsc_wib {}
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 3]
	lappend poss {*}[list_cor $header {count offset span lowerLimit dataRange file}]
	if {[lsearch $poss -1] != -1} {
		error "necessary field not found while converting ucsc wiggle track"
	}
	set progress 10000
	set format "%.${precision}f"
	set pvalue NaN
	set pbegin NaN
	set pend NaN
	while 1 {
		if {[eof $f]} break
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		foreach {chrom begin end count offset span lowerLimit dataRange file} [list_sub $line $poss] break
		if {$begin != $pend} {
			if {$pvalue ne NaN} {
				putsresult $o $chrom $pbegin $cur $pvalue $pnum
			}
			set pvalue NaN
			set pbegin NaN
		}
		set pend $end
		if {$progress >= 10000} {
			puts $chrom:$begin
			set progress 0
		}
		incr progress
		if {$file ne $ucsc_wib} {
			catch {close $b}
			set ucsc_wib $file
			set file [ucsc_wibfile $file $dir]
			puts "Opening $file"
			set b [open $file]
			fconfigure $b -translation binary
		}
		seek $b $offset
		set bin [read $b $count]
		if {![binary scan $bin cu$count list]} {
			error "error scanning wib file"
		}
		set cur $begin
		foreach value $list {
			if {$value != 128} {
				set value [expr {$lowerLimit + $dataRange * $value / 127.0}]
				if {$useaddnum} {
					if {$value >= $addnum} {set num 1} else {set num 0}
				}
				set value [formula $value]
				set value [format $format $value]
			} else {
				set value NaN
				set num 0
			}
			if {$value ne $pvalue} {
				if {$pvalue ne NaN} {
					putsresult $o $chrom $pbegin $cur $pvalue $pnum
				}
				set pvalue $value
				set pbegin $cur
				set pnum $num
			}
			incr cur $span
		}
	}
	if {$pvalue ne NaN} {
		putsresult $o $chrom $pbegin $cur $pvalue $pnum
	}
	catch {close $f}; catch {close $b}; catch {close $o}
	file rename -force $resultfile.temp $resultfile
}

proc cg_ucscwiggle2reg {args} {
	set len [llength $args]
	set precision 1
	set formula {$value}
	set addnum {}
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-p {
				set precision $value
			}
			-f {
				set formula $value
			}
			-n {
				set addnum $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {($len < 1) && ($len > 2)} {
		puts "Wrong number of arguments"
		errorformat ucscwiggle2reg
		exit 1
	}
	foreach {file resultfile} $args break
	if {$resultfile eq ""} {
		set resultfile reg_$file
	}
	if {[file exists $resultfile]} {
		puts "Skipping $resultfile: file exists"
		return
	}
	ucscwiggle2reg $file $resultfile $precision $formula $addnum
}

proc ucscwb2reg {file resultfile {precision 1} {formula {}} {addnum {}}} {
	if {$formula eq ""} {
		proc formula {value} {return $value}
	} else {
		proc formula {value} "return \[expr \{$formula\}\]"
	}
	if {$addnum ne ""} {
		proc putsresult {o args} {
			puts $o [join $args \t]
		}
		set useaddnum 1
	} else {
		proc putsresult {o args} {
			puts $o [join [lrange $args 0 3] \t]
		}
		set useaddnum 0
	}
	set dir [file dir [file normalize $file]]
	set f [open $file]
	gets $f
	set wbfile [gets $f]
	close $f
	puts "Making $resultfile"
	set wbfile [ucsc_wibfile $wbfile $dir]
	if {![file exists $wbfile.bedgraph]} {
		exec bigWigToBedGraph $wbfile $wbfile.bedgraph.temp
		file rename -force $wbfile.bedgraph.temp $wbfile.bedgraph
	}
	set f [gzopen $wbfile.bedgraph]
	set o [open $resultfile.temp w]
	putsresult $o chromosome begin end score num
	set num {}
	set format "%.${precision}f"
	set pvalue NaN
	set pbegin NaN
	set pend NaN
	while 1 {
		if {[eof $f]} break
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		foreach {chrom begin end value} $line break
		if {$value >= $addnum} {set num 1} else {set num 0}
		set value [formula $value]
		set value [format $format $value]
		if {$begin != $pend || $value ne $pvalue} {
			if {$pvalue ne NaN} {
				putsresult $o $chrom $pbegin $pend $pvalue $pnum
			}
			set pbegin $begin
			set pvalue $value
			set pnum $num
		}
		set pend $end
	}
	if {$pvalue ne NaN} {
		putsresult $o $chrom $pbegin $pend $pvalue $pnum
	}
	catch {close $f}; catch {close $o}
	file rename -force $resultfile.temp $resultfile
}

proc cg_ucscwb2reg {args} {
	set len [llength $args]
	set precision 1
	set formula {$value}
	set addnum {}
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-p {
				set precision $value
			}
			-f {
				set formula $value
			}
			-n {
				set addnum $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {($len < 1) && ($len > 2)} {
		puts "Wrong number of arguments"
		errorformat ucscwiggle2reg
		exit 1
	}
	foreach {file resultfile} $args break
	if {$resultfile eq ""} {
		set resultfile [file join [file dir $file] reg_[file tail $file]]
	}
	if {[file exists $resultfile]} {
		puts "Skipping $resultfile: file exists"
		return
	}
	ucscwb2reg $file $resultfile $precision $formula $addnum
}

