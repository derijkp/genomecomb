#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg {args} {
	# puts "cg $args"
	eval exec cg $args 2>@ stderr
}

proc opencgifile {file headerVar {numlinesVar {}}} {
	if {$numlinesVar ne ""} {
		upvar $numlinesVar numlines
	}
	set numlines 0
	global cache
	upvar $headerVar header
	if {[file extension $file] eq ".gz"} {
		set f [open "|zcat $file"]
	} else {
		set f [open $file]
	}
	while {![eof $f]} {
		set line [gets $f]
		incr numlines
		if {[string length $line] && [string index $line 0] ne "#"} break
	}
	if {[string index $line 0] eq ">"} {
		set header [string range $line 1 end]
	} else {
		set header $line
	}
	while {![eof $f]} {
		set line [split [gets $f] \t]
		incr numlines
		if {[llength $line]} break
	}
	set cache($f) $line
	incr numlines -1
	return $f
}

proc opencgifile {file headerVar {numlinesVar {}}} {
	if {$numlinesVar ne ""} {
		upvar $numlinesVar numlines
	}
	set numlines 0
	global cache
	upvar $headerVar header
	if {[file extension $file] eq ".gz"} {
		set f [open "|zcat $file"]
	} else {
		set f [open $file]
	}
	while {![eof $f]} {
		set line [gets $f]
		incr numlines
		if {[string length $line] && [string index $line 0] ne "#"} break
	}
	if {[string index $line 0] eq ">"} {
		set header [string range $line 1 end]
	} else {
		set header $line
	}
	while {![eof $f]} {
		set line [split [gets $f] \t]
		incr numlines
		if {[llength $line]} break
	}
	set cache($f) $line
	incr numlines -1
	return $f
}

proc cggets {f} {
	global cache
	set line {}
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {[llength $line]} break
	}
	set result $cache($f)
	set cache($f) $line
	return $result
}

proc getnotempty {f} {
	set line {}
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {[llength $line]} break
	}
	return $line
}

proc assert {check message} {
	if {![uplevel expr $check]} {
		error $message
	}
}

proc readcgimap {f} {
	global cache
	set line $cache($f)
	while {![eof $f]} {
		foreach {flags chromosome offsetInChr gap1 gap2 gap3 weight mateRec} $line break
		set end [expr {$offsetInChr + 35 + $gap1 + $gap2 + $gap3}]
		binary scan $weight c weight
		incr weight -33
		set lastdnbrecord [expr {$flags & 0x01}]
		if {[expr {$flags & 0x02}]} {set side r} else {set side l}
		if {[expr {$flags & 0x04}]} {set strand -} else {set strand +}
		lappend list [list $flags $chromosome $offsetInChr $end $strand $side $gap1 $gap2 $gap3 $weight $mateRec]
		set line [split [gets $f] \t]
		if {$lastdnbrecord} {
			set cache($f) $line
			return $list
		}
	}
}

proc region2bins {start end {level 0}} {
	if {$level == 0} {
		set bins 0
		set stop 0
	} else {
		set bins {}
		array set levels {1 1 2 9 3 73 4 585 5 4681}
		set stop $levels($level)
	}
	set pre 4681
	set start [expr {$start >> 14}]
	set end [expr {$end >> 14}]
	while {$pre > $stop} {
		for {set i $start} {$i <= $end} {incr i} {
			lappend bins [expr {$pre+$i}]
		}
		set pre [expr {$pre >> 3}]
		set start [expr {$start >> 3}]
		set end [expr {$end >> 3}]
	}
	return [lsort -integer $bins]
}

proc histogram {list aVar} {
	upvar $aVar a
	unset -nocomplain a
	foreach el $list {
		if {![info exists a($el)]} {
			set a($el) 1
		} else {
			incr a($el)
		}
	}
}

proc max {args} {
	lmath_max [list_remove $args {} - ?]
}

proc min {args} {
	lmath_min [list_remove $args {} - ?]
}

proc opensqlite3 {dbfile query} {
	set f [open "| sqlite3 -separator \"\t\" $dbfile \"$query\""]
}

proc cg_bgzip args {
	set pos 0
	set keep 0
	foreach {key value} $args {
		switch -- $key {
			-k {
				set keep 1
			}
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	foreach file $args {
		set ext [file extension $file]
		switch $ext {
			.gz {
				set error [catch {exec tabix $file} errormsg]
				if {[regexp {was bgzip used to compress} $errormsg]} {
					putslog "bgzip $file"
					exec gunzip -d -c $file > $file.temp2
					exec bgzip -c $file > $file.temp
					file delete $file.temp2
					if {$keep} {file rename -force $file $file.old}
					file rename -force $file.temp $file
				}
			}
			.rz {
				putslog "bgzip $file"
				set result [file root $file].gz
				exec razip -d -c $file > $result.temp2
				exec bgzip -c $result.temp2 > $result.temp
				file delete $result.temp2
				file rename -force $result.temp $result
				if {!$keep} {file delete $file}
			}
			.bz2 {
				putslog "bgzip $file"
				set result [file root $file].gz
				exec bzcat $file > $result.temp2
				exec bgzip -c $result.temp2 > $result.temp
				file delete $result.temp2
				file rename -force $result.temp $result
				if {!$keep} {file delete $file}
			}
			default {
				putslog "bgzip $file"
				exec bgzip -c $file > $file.gz.temp
				file rename -force $file.gz.temp $file.gz
				if {!$keep} {file delete $file}
			}
		}
	}
}

proc cg_razip args {
	set pos 0
	set keep 0
	foreach {key value} $args {
		switch -- $key {
			-k {
				set keep 1
			}
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	foreach file $args {
		set ext [file extension $file]
		switch $ext {
			.gz {
				putslog "razip $file"
				set result [file root $file].rz
				exec gunzip -d -c $file > $result.temp2
				exec razip -c $result.temp2 > $result.temp
				file delete $result.temp2
				file rename -force $result.temp $result
				if {!$keep} {file delete $file}
			}
			.rz {
				putslog "$file already razip"
			}
			.bz2 {
				putslog "razip $file"
				set result [file root $file].rz
				exec bzcat $file > $result.temp2
				exec razip -c $result.temp2 > $result.temp
				file delete $result.temp2
				file rename -force $result.temp $result
				if {!$keep} {file delete $file}
			}
			default {
				putslog "razip $file"
				exec razip -c $file > $file.rz.temp
				file rename -force $file.rz.temp $file.rz
				if {!$keep} {file delete $file}
			}
		}
	}
}

proc cg_unzip args {
	foreach file $args {
		putslog "Uncompressing $file"
		gunzip $file
	}
}

proc gzopen {file {pos -1}} {
	if {![file exists $file]} {
		exiterror "Error: couldn't open \"$file\": no such file or directory"
	}
	if {[inlist {.rz} [file extension $file]]} {
		if {$pos == -1} {
			set f [open "| razip -d -c $file"]
		} else {
			set f [open "| razip -d -c -b $pos $file"]
		}
	} elseif {[inlist {.bgz .gz} [file extension $file]]} {
		if {$pos == -1} {
			set f [open "| zcat $file"]
		} else {
			error "positioning not supported in (b)gz files"
		}
	} elseif {[inlist {.bz2} [file extension $file]]} {
		if {$pos == -1} {
			set f [open "| bzcat $file"]
		} else {
			error "positioning not supported in bz2 files"
		}
	} else {
		set f [open $file]
		if {$pos != -1} {
			seek $f $pos
		}
	}
	return $f
}

proc gzclose {f} {
	if {[catch {close $f} error]} {
		if {$error eq "child killed: write on pipe with no readers"} return
		error $error
	}
}

proc gzroot filename {
	if {[inlist {.rz .gz .bgz .bz2} [file extension $filename]]} {
		return [file root $filename]
	} else {
		return $filename
	}
}

proc gziscompressed filename {
	if {[inlist {.rz .gz .bgz .bz2} [file extension $filename]]} {
		return 1
	} else {
		return 0
	}
}

proc gzfile {args} {
	foreach filename $args {
		if {![catch {glob $filename $filename.rz $filename.bgz $filename.gz $filename.bz2} list]} {
			return [lindex $list 0]
		}
	}
	return [lindex $args 0]
}

proc gzfiles {args} {
	set result {}
	foreach filename $args {
		if {![catch {glob $filename $filename.rz $filename.bgz $filename.gz $filename.bz2} list]} {
			lappend result {*}$list
		}
	}
	return $result
}

proc gzcat {filename} {
	switch [file extension $filename] {
		.rz {set cat "razip -d -c"}
		.gz - .bgz {set cat zcat}
		.bz2 {set cat bzcat}
		default {set cat cat}
	}
	return $cat
}

proc gztemp {filename} {
	set ext [file extension $filename]
	switch $ext {
		.gz {
			set tempfile [tempfile]
			exec gunzip -d -c $filename > $tempfile
			set ::gztemp_files($tempfile) 1
			return $tempfile
		}
		.rz {
			set tempfile [tempfile]
			exec razip -d -c $filename > $tempfile
			set ::gztemp_files($tempfile) 1
			return $tempfile
		}
		.bz2 {
			set tempfile [tempfile]
			exec bzcat $filename > $tempfile
			set ::gztemp_files($tempfile) 1
			return $tempfile
		}
		default {
			set ::gztemp_files($filename) 0
			return $filename
		}
	}
}

proc gzrmtemp {filename} {
	if {$::gztemp_files($filename)} {
		file delete $filename
	}
}

proc overlap {start1 end1 start2 end2} {
	if {$start2 >= $end1} {return [expr {$end1-$start2}]}
	if {$end2 < $start1} {return [expr {$end2-$start1}]}
	if {$start2 > $start1} {set start1 $start2}
	if {$end2 < $end1} {set end1 $end2}
	expr {$end1-$start1}
}

proc timestamp {} {
	clock format [clock seconds] -format "%Y-%m-%d %H:%M:%S"
}

proc chrindexseek {file f chr} {
	set indexfile [gzroot $file].chrindex
	if {![file exists $indexfile]} {
		set tf [gzopen $file]
		set header [gets $tf]
		set chrpos [lsearch $header chromosome]
		set prevchr {}
		set list {}
		set o [open $indexfile w]
		while {![eof $tf]} {
			set pos [tell $tf]
			set line [gets $tf]
			set chr [lindex $line $chrpos]
			if {$chr ne $prevchr} {
				puts $o $chr\t$pos
				set prevchr $chr
			}
		}
		close $tf
		close $o
	}
	set trfchrpos [split [string trim [file_read $indexfile]] \n\t]
	if {[dict exists $trfchrpos chr$chr]} {
		set chr chr$chr
	}
	if {[catch {set fpos [dict get $trfchrpos $chr]}]} {
		seek $f 0 end
	} else {
		seek $f $fpos start
	}
}

proc binsearch {table index value} {
	set begin 0
	set end [llength $table]
	while 1 {
		set mid [expr {($begin+$end)/2}]
		if {$mid == $begin} break
		set v [lindex $table $mid $index]
		if {$value < $v} {
			set end $mid
		} elseif {$value > $v} {
			set begin $mid
		} else {
			break
		}
	}
	return $mid
}

proc ifcatch {command varName args} {
	upvar $varName result
	set error [uplevel [list catch $command $varName]]
	if {!$error} {
		return $result
	}
	set switchlist [list_pop args]
	if {[lindex $switchlist end-1] ne "default"} {
		lappend switchlist default "error [list $result]"
	}
	uplevel switch $args [list $result $switchlist]
}

proc gunzip {file args} {
	if {[llength $args]} {
		set resultfile [lindex $args 0]
	} else {
		set resultfile [file root $file]
	}
	set error [catch {exec gunzip -S [file ext $file] -c $file > $resultfile.temp} result]
	if $error {
		if {![regexp "decompression OK, trailing garbage ignored" $result]} {
			error $result
		}
	}
	file rename -force $resultfile.temp $resultfile
	if {![llength $args]} {
		file delete $file
	}
	return $result
}

#doc scratchdir title {
#scratchdir
#} shortdescr {
# returns a directory in which temporary files can be stored. This directory is specific to one proces:
# no other processes will (should) write in this directory. Subsequent calls to the function within one process
# will allways be the same directory, The program has to take care not to overwrite its own files
# Temporary files returned by tempfile are also in this directory (named like _Extral_temp_1.tmp)
# The program should also not overwrite these
# The temporary directory is deleted when the program exits by an atexit handler
# 
#}

proc scratchdir {} {
	global env
	if {![info exists env(SCRATCHDIR)]} {
		putslog "Could not find SCRATCHDIR, using tempdir"
		return [tempdir]
	}
	if {![info exists ::Extral::scratchdir]} {
		for {set i 0} {$i < 20} {incr i} {
			set testdir [file join $env(SCRATCHDIR) scratchExtral.[pid]-[Extral::randstring 20]]
			if {[file exists $testdir]} continue
			if {[catch {
				file mkdir $testdir
				if {$::tcl_platform(platform) eq "unix"} {
					file attributes $testdir -permissions 0700
				}
				set files [glob -nocomplain $testdir/*]
				if {[llength $files]} {
					error "Very fishy: there are files in the temporary directory I just created"
				}
				set ::Extral::scratchdir $testdir
				set ::Extral::scratchnum 0
			}]} continue
			break
		}
	}
	if {![info exists ::Extral::scratchdir]} {
		error "couldn't create temporary directory"
	}
	return $::Extral::scratchdir
}

proc scratchfile {{action {get}} {type file}} {
	switch $action {
		get {
			set scratchdir [scratchdir]
			return [file join $scratchdir _Extral_scratch_[incr ::Extral::scratchnum].tmp]
		}
		clean {
			set scratchdir [scratchdir]
			catch {file delete -force $scratchdir}
			unset ::Extral::scratchdir
		}
		cleanall {
			catch {eval file delete -force [glob [file join $scratch_dir scratchExtral*]]}
		}
		default {
			return -code error "bad option \"$action\": must be get, clean or cleanall"
		}
	}
}

proc chanexec {in out pipe} {
	set o [open "|\ $pipe\ >@\ $out 2>@\ stderr" w]
	if {[info exists ::filebuffer($in)]} {
		foreach line $::filebuffer($in) {
			puts $o $line
		}
		unset ::filebuffer($in)
	}
	fcopy $in $o
	if {$in ne "stdin"} {catch {close $in}}
	if {$out ne "stdout"} {catch {close $out}}
	close $o
}

proc wgetfile {url {resultfile {}}} {
	if {$resultfile eq ""} {
		set resultfile [file tail $url]
	}
	catch {exec wget --tries=45 -O $resultfile.temp $url} errmsg
	if {![file exists $resultfile.temp]} {
		return {}
	}
	if {[regexp "No such file" $errmsg]} {
		file delete $resultfile.temp
		return {}
	}
	file rename -force $resultfile.temp $resultfile
	return $resultfile
}

proc cg_checktsv {file} {
	set f [gzopen $file]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 4 0]
	set line [split [gets $f] \t]
	# foreach {pchr pbegin pend ptype} [list_sub $line $poss] break
	set prev [list_sub $line $poss]
	set len [llength $header]
	set linenr 0
	while 1 {
		if {[llength $line] != $len} {
			puts "line $linenr is of wrong length: [llength $line] iso $len\t$line"
		}
		if {[eof $f]} break
		set line [split [gets $f] \t]
		set llen [llength $line]
		if {!$llen && [eof $f]} break
		incr linenr
		set cur [list_sub $line $poss]
		if {[list $prev $cur] ne [lsort -dict [list $prev $cur]]} {
			puts "line $linenr is sorted wrong:\t$line"
		}
		set prev $cur
	}
	close $f
}

proc progress {cmd args} {
	if {![llength [info commands winfo]]} {
		global progresslevel
		if {![info exists progresslevel]} {
			set progresslevel -1
		}
		# no Tk, so we are running commandline: ignore a lot of the progress commands
		switch $cmd {
			start {
				incr progresslevel
			}
			stop {
				if {$progresslevel == -1} return
				incr progresslevel -1
			}
			protect {
				foreach {code on_error} $args break
				set error [catch {uplevel 1 $code} result]
				if {$error} {
					if {$on_error eq ""} {
						set errorInfo $::errorInfo
						return -code error -errorinfo $errorInfo $result
					} else {
						set ::errorResult $result
						uplevel 1 $on_error
					}
				} else {
					return $result
				}
			}
		}
	} else {
		switch $cmd {
			onerror {
				Classy::Progress on_error [subst {
					[lindex $args 0]
				}]
			}
			startdisplay {
				if {([Classy::Progress level] == -1)} {
				}
				uplevel 1 [list Classy::Progress start] $args
			}
			start {
				if {([Classy::Progress level] == -1)} {
				}
				uplevel 1 [list Classy::Progress $cmd] $args
			}
			stop {
				Classy::Progress stop
				if {([Classy::Progress level] == -1)} {
					Classy::Progress on_error {}
				}
			}
			default {
				uplevel 1 [list Classy::Progress $cmd] $args
			}
		}
	}
}

proc exiterror errormessage {
	puts stderr $errormessage
	exit 1
}

proc catprog file {
	set ext [file extension $file]
	switch $ext {
		.rz - .gz - .bgz {
			return zcat
		}
		.bz2 {
			return bzcat
		}
		default {
			return cat
		}
	}
}

proc getline f {
	set line [split [gets $f] \t]
	while {![llength $line] && ![eof $f]} {
		set line [split [gets $f] \t]
	}
	return $line
}

proc mklink {src dest} {
	set src [file normalize $src]
	set dest [file normalize $dest]
	set pos 0
	set ssrc [file split $src]
	set sdest [file split $dest]
	# puts $ssrc\n$sdest
	foreach s $ssrc d $sdest {
		if {$s ne $d} break
		incr pos
	}
	if {$pos > 1} {
		set prelen [expr {[llength $sdest]-$pos -1}]
		set src [file join {*}[list_fill $prelen ..] {*}[lrange $ssrc $pos end]]
	}
	file link -symbolic $dest $src
}

proc gzmklink {src dest} {
	set src [gzfile $src]
	set ext_s [file extension $src]
	set ext_d [file extension $dest]
	if {$ext_s ne $ext_d && [inlist {.gz .bgz .rz .bz2} $ext_s]} {
		mklink $src $dest$ext_s
	} else {
		mklink $src $dest
	}
}

if 0 {

	ifcatch {error test} result {
		test1 {puts ERRORtest1}
		test2 {puts ERRORtest2}
		default {puts ERRORdefault}
	}

	ifcatch {error test2} result {
		test1 {puts ERRORtest1}
		test2 {puts ERRORtest2}
	}

	ifcatch {error test} result {
		test1 {puts ERRORtest1}
		test2 {puts ERRORtest2}
	}

	ifcatch {set a 1} result {
		test1 {puts ERRORtest1}
		test2 {puts ERRORtest2}
	}

}

proc reload {} {
	global appdir
	foreach file [glob $appdir/lib/*.tcl $appdir/lib-exp/*.tcl] {
		puts "sourcing $file"
		if {[catch {source $file} e]} {
			puts "error: $e"
		}
	}
}
