#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg {cmd args} {
	# puts "cg $args"
	if {[info exists ::stderr_redirect] && [lsearch -regexp $args 2>.?] == -1} {
		set tempfile $::stderr_redirect
		set redir [list 2> $::stderr_redirect]
	} else {
		set tempfile [tempfile]
		set redir ""
	}
	if {[string length $args] < 2000} {
		set error [catch {exec cg $cmd {*}$args {*}$redir} result]
	} else {
		set temprunfile [tempfile]
		set poss [list_concat [list_find -glob $args ">*"] [list_find -glob $args "<*"]]
		if {[llength $poss]} {
			set pos [min $poss]
			set redirect [lrange $args $pos end]
			set code [lrange $args 0 [expr {$pos-1}]]
		} else {
			set code $args
			set redirect {}
		}
		file_write $temprunfile [list cg_$cmd {*}$code]\n
		set error [catch {exec cg source $temprunfile {*}$redirect {*}$redir} result]
	}
	if {$error} {
		if {[file exists $tempfile]} {
			set errmessage [file_read $tempfile]\n
		} else {
			set errmessage {}
		}
		if {[info exists ::stderr_redirect]} {file delete $tempfile}
		return -code error $errmessage$result
	} else {
		if {[info exists ::stderr_redirect]} {file delete $tempfile}
		return $result
	}
}

proc bgcg_progress {bgexechandleVar args} {
	upvar #0 $bgexechandleVar bgexechandle
	if {![isint $args]} {
		append ::bgerror [lindex $args 0]\n
		return
	}
	if {[catch {
		progress set $args
	}]} {
		puts error
		Extral::bgexec_cancel $bgexechandle
	}
}

proc bgcg {progresscommand channelvar cmd args} {
	# puts "progresscommand cg $args"
	if {[info exists ::stderr_redirect]} {
		set tempfile $::stderr_redirect
	} else {
		set tempfile [tempfile]
	}
	if {[string length $args] < 2000} {
		set ::bgerror {}
		Extral::bgexec -progresscommand [list $progresscommand $channelvar] -no_error_redir -channelvar $channelvar \
				cg $cmd {*}$args 2>@1
		if {$::bgerror ne ""} {error $::bgerror}
	} else {
		set poss [list_concat [list_find -glob $args ">*"] [list_find -glob $args "<*"]]
		if {[llength $poss]} {
			set pos [min $poss]
			set redirect [lrange $args $pos end]
			set code [lrange $args 0 [expr {$pos-1}]]
		} else {
			set code $args
			set redirect {}
		}
		set temprunfile [tempfile]
		file_write $temprunfile [list cg_$cmd {*}$code]\n
		set ::bgerror {}
		Extral::bgexec -progresscommand [list $progresscommand $channelvar] -no_error_redir -channelvar $channelvar \
				cg source $temprunfile {*}$redirect 2>@1
		if {$::bgerror ne ""} {error $::bgerror}
	}
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
	if {[llength $args] == 1} {set args [lindex $args 0]}
	lmath_max [list_remove $args {} - ?]
}

proc min {args} {
	if {[llength $args] == 1} {set args [lindex $args 0]}
	lmath_min [list_remove $args {} - ?]
}

proc opensqlite3 {dbfile query} {
	set f [open "| sqlite3 -separator \"\t\" $dbfile \"$query\""]
}

proc compress {file {ext .lz4}} {
	set file [file normalize $file]
	if {[file exists $file$ext]} {file delete $file$ext}
	if {[inlist {.rz} $ext]} {
		exec razip $file
	} elseif {[inlist {.lz4} $ext]} {
		exec lz4c -9 -c $file > $file.lz4.temp
		file rename -force $file.lz4.temp $file.lz4
	} elseif {[inlist {.gz} $ext]} {
		exec gzip $file
	} elseif {[inlist {.bgz} $ext]} {
		exec bgzip $file
	} elseif {[inlist {.bz2} $ext]} {
		exec bzip2 $file
	} else {
		error "Unknown extension $ext"
	}
}

proc gzopen {file {pos -1}} {
	if {![file exists $file]} {
		exiterror "Error: couldn't open \"$file\": no such file or directory"
	}
	set file [file normalize $file]
	set ext [file extension $file]
	if {[inlist {.rz} $ext]} {
		if {$pos == -1} {
			set f [open "| razip -d -c $file"]
		} else {
			set f [open "| razip -d -c -b $pos $file"]
		}
	} elseif {[inlist {.lz4} $ext]} {
		if {$pos == -1} {
			set f [open "| lz4c -d -c $file"]
		} else {
			error "positioning not supported in lz4 files"
		}
	} elseif {[inlist {.bgz .gz} $ext]} {
		if {$pos == -1} {
			set f [open "| zcat $file"]
		} else {
			error "positioning not supported in (b)gz files"
		}
	} elseif {[inlist {.bz2} $ext]} {
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
		if {[regexp {Successfully decoded [0-9]+ bytes} $error]} return
		error $error
	}
}

proc gzcatch {cmd} {
	if {[catch {uplevel $cmd} error]} {
		if {$error eq "child killed: write on pipe with no readers"} return
		if {[regexp {Successfully decoded [0-9]+ bytes} $error]} return
		error $error
	}
}

proc gzroot filename {
	if {[inlist {.rz .lz4 .gz .bgz .bz2} [file extension $filename]]} {
		return [file root $filename]
	} else {
		return $filename
	}
}

proc gziscompressed filename {
	if {[inlist {.rz .lz4 .gz .bgz .bz2} [file extension $filename]]} {
		return 1
	} else {
		return 0
	}
}

proc gzfile {args} {
	foreach filename $args {
		if {![catch {glob $filename $filename.rz $filename.lz4 $filename.bgz $filename.gz $filename.bz2} list]} {
			return [lindex $list 0]
		}
	}
	return [lindex $args 0]
}

proc gzexists {filename {checkcompressed 1}} {
	if {$checkcompressed} {
		expr {[file exists $filename] || [file exists $filename.rz] || [file exists $filename.lz4] || [file exists $filename.gz] ||[file exists $filename.bgz] || [file exists $filename.bz2]}
	} else {
		file exists $filename
	}
}

proc checkfile {args} {
	foreach filename $args {
		if {![catch {glob $filename} list]} {
			return [lindex $list 0]
		}
	}
	return [lindex $args 0]
}

proc gzfiles {args} {
	set result {}
	foreach filename $args {
		if {![catch {glob $filename $filename.rz $filename.lz4 $filename.bgz $filename.gz $filename.bz2} list]} {
			lappend result {*}$list
		}
	}
	return $result
}

proc gzarraynames {aVar pattern} {
	upvar $aVar a
	set result [lsort [list_remdup [list_concat [array names a $pattern] [array names a $pattern.rz] [array names a $pattern.lz4] [array names a $pattern.gz] [array names a $pattern.bgz] [array names a $pattern.bz2]]]]
	return $result
}

proc checkfiles {args} {
	set result {}
	foreach filename $args {
		if {![catch {glob $filename} list]} {
			lappend result {*}$list
		}
	}
	return $result
}

proc gzcat {filename} {
	switch [file extension $filename] {
		.rz {set cat "razip -d -c"}
		.lz4 {set cat "lz4c -d -c"}
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
			set tempfile [scratchfile get]
			exec gunzip -d -c $filename > $tempfile
			set ::gztemp_files($tempfile) 1
			return $tempfile
		}
		.rz {
			set tempfile [scratchfile get]
			exec razip -d -c $filename > $tempfile
			set ::gztemp_files($tempfile) 1
			return $tempfile
		}
		.lz4 {
			set tempfile [scratchfile get]
			exec lz4c -d -c $filename $tempfile
			set ::gztemp_files($tempfile) 1
			return $tempfile
		}
		.bz2 {
			set tempfile [scratchfile get]
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

proc reg_compare {loc1 loc2} {
	if {![llength $loc1]} {return 1}
	if {![llength $loc2]} {return -1}
	foreach {chr1 start1 end1} $loc1 break
	foreach {chr2 start2 end2} $loc2 break
	set chrcomp [loc_compare $chr1 $chr2]
	if {$chrcomp != 0} {return $chrcomp}
	if {$start1 == $start2} {return 0}
	if {$start2 >= $end1} {return [expr {$end1 - $start2 -1}]}
	if {$end2 <= $start1} {return [expr {$start1 - $end2 + 1}]}
	return 0
}

proc timestamp {} {
	clock format [clock seconds] -format "%Y-%m-%d %H:%M:%S"
}

proc chrindexseek {file f chr} {
	set root [gzroot $file]
	file mkdir $root.index
	set indexfile [indexdir_file $file chrindex ok]
	if {!$ok} {
		set tf [gzopen $file]
		set header [gets $tf]
		set chrpos [tsv_basicfields $header 1 0]
		set prevchr {}
		set list {}
		set o [open $indexfile w]
		while {![eof $tf]} {
			set pos [tell $tf]
			set line [gets $tf]
			set chr [chr_clip [lindex $line $chrpos]]
			if {$chr ne $prevchr} {
				puts $o $chr\t$pos
				set prevchr $chr
			}
		}
		close $tf
		close $o
	}
	set trfchrpos [split [string trim [file_read $indexfile]] \n\t]
	set chr [chr_clip $chr]
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

proc decompress {file args} {
	if {[llength $args]} {
		set resultfile [lindex $args 0]
	} else {
		set resultfile [file root $file]
	}
	if {[file extension $file] eq ".lz4"} {
		set error [catch {exec lz4c -d $file > $resultfile.temp} result]
	} else {
		set error [catch {exec zcat $file > $resultfile.temp} result]
	}
	if $error {
		if {![regexp "decompression OK, trailing garbage ignored" $result] && ![regexp {Successfully decoded} $errmessage]} {
			error $result
		}
	}
	file rename -force $resultfile.temp $resultfile
	if {![llength $args]} {
		file delete $file
	}
	return $result
}

proc gunzip {file args} {
	if {[llength $args]} {
		set resultfile [lindex $args 0]
	} else {
		set resultfile [file root $file]
	}
	set error [catch {exec zcat $file > $resultfile.temp} result]
	if $error {
		if {![regexp "decompression OK, trailing garbage ignored" $result] && ![regexp {Successfully decoded} $result]} {
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
		# putslog "Could not find SCRATCHDIR, using tempdir"
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
		error "couldn't create scratch directory in $env(SCRATCHDIR) (defined by env variable SCRATCHDIR)"
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

proc file_add {file args} {
	set f [open $file a]
	foreach arg $args {
		puts $f $arg
	}
	close $f
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
	if {[catch {
		exec wget -c --tries=45 -O $resultfile.temp $url 2>@ stderr
	} errmsg]} {
		if {[file size $resultfile.temp] == 0} {
			file delete $resultfile.temp
		}
		return {}
	}
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
			next {
				uplevel 1 [list Classy::Progress $cmd] $args
				update idletasks
			}
			cancel {
				Classy::Progress cancel
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
		.lz4 {
			return {lz4c -d}
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

proc find_link file {
	set file [file_absolute $file]
	while 1 {
		if {[catch {
			set file [file join [file dir $file] [file readlink $file]]
			set file [file_absolute $file]
		}]} break
	}
	return $file
}

proc file_link {linkname linkdest} {
	set dir [file dir $linkname]
	if {[file pathtype $linkdest] eq "absolute"} {
		set dest $linkdest
	} else {
		set dest [file join $dir $linkdest]
	}
	if {![file exists $dest]} {
		file_write $dest temp
		file link $linkname $linkdest
		file delete $dest
	} else {
		file link $linkname $linkdest
	}
}

proc mklink {src dest} {
	set src [file_absolute $src]
	set dest [file_absolute $dest]
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
	set err [catch {file link $dest} link]
	if {!$err || $link ne "$src"} {
		file delete $dest
	}
	if {![file exists $dest]} {
		if {[file exists $src]} {
			file link -symbolic $dest $src
		} else {
			set keeppwd [pwd]
			cd [file dir $dest]
			exec ln -s $src [file tail $dest]
			cd $keeppwd
		}
	}
}

proc gzmklink {src dest} {
	set src [gzfile $src]
	set ext_s [file extension $src]
	set ext_d [file extension $dest]
	if {$ext_s ne $ext_d && [inlist {.gz .bgz .rz .lz4 .bz2} $ext_s]} {
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

proc dict_get_default {d key {default {}}} {
	if {[catch {dict get $d $key} result]} {
		return $default
	} else {
		return $result
	}
}

proc samples {header {pattern {}}} {
	set names {}
	foreach col $header {
		set pos [string first - $col]
		if {$pos != -1} {
			incr pos
			set name [string range $col $pos end]
			if {$pattern eq "" || [string match $pattern $name]} {
				lappend names $name
			}
		}
	}
	list_remdup $names
}

proc genomecombenv {} {
	global auto_path env appdir tcl_dirtcl genomecombdir externdir
	if {![info exists appdir]} {
		set appdir ${genomecomb::dir}
	}
	if {[file dir [file dir $appdir]] eq [get tcl_dirtcl ""]} {
		# we are being run from a dirtcl installation in apps/cg
		set genomecombdir $tcl_dirtcl
		set externdir $genomecombdir/bin
		set env(PATH) $appdir/bin:$externdir:$genomecombdir:$env(PATH)
	} elseif {[file tail $appdir] eq "cg_viz"} {
		# we are being run from dev cg_viz
		set genomecombdir [file dir $appdir]
		set externdir $genomecombdir/extern
		set env(PATH) $genomecombdir/bin:$externdir:$genomecombdir:$env(PATH)
	} else {
		# we are being run from dev 
		set genomecombdir $appdir
		set externdir $genomecombdir/extern
		set env(PATH) $genomecombdir/bin:$externdir:$genomecombdir:$env(PATH)
	}
	return $genomecombdir
}

proc trans {trans value} {
	if {[dict exists $trans $value]} {
		return [dict get $trans $value]
	} {
		return $value
	}
}

proc lforeach {args} {
	set result {}
	set pattern [list_pop args]
	set code "lappend result \"$pattern\""
	foreach {*}$args $code
	return $result
}

proc sourcename base {
	if {![regexp {^[^-]*-(.+)$} $base temp name]} {
		set name $base
	}
	return $name
}

proc razip_job {file args} {
	set deps [list $file {*}$args]
	uplevel [list job razip-$file -checkcompressed 0 -deps $deps -targets $file.rz -rmtargets $file -code {
		if {![file exists $dep]} {error "error compressing: file $dep does not exist"}
		cg_razip $dep
	}]
}

proc file_absolute {file} {
	if {[string range $file 0 1] eq "~/"} {
		set result [file join $::env(HOME) [string range $file 2 end]]
	} else {
		set result {}
		foreach el [file split [file join [pwd] $file]] {
			if {$el eq ".."} {
				if {[llength $result] <= 1} {error "file_absolute error: cannot .. past root"}
				list_pop result
			} elseif {$el ne "." && $el ne ""} {
				lappend result $el
			}
		}
	}
	file join {*}$result
}

proc file_tempwrite {file} {
	if {![file exists $file.temp]} {return $file.temp}
	set num 2
	while {[file exists $file.temp$num]} {incr num}
	return $file.temp$num
}

proc trimformat args {
	string trimright [string trimright [::format {*}$args] 0] .
}

proc oargserr {cmd def} {
	set result [list $cmd]
	foreach el $def {
		set field [lindex $el 0]
		if {[llength $el] > 1} {
			lappend result ?$field?
		} else {
			lappend result $field
		}
	}
	return [join $result " "]
}

proc oargs {cmd def arg} {
	set parseopts 1
	set pos 0
	set len [llength $def]
	if {[lindex $def end] eq "args"} {
		set useargs 1
		list_pop def
		incr len -1
	} else {
		set useargs 0
	}
	foreach el $def {
		if {[llength $el] > 1} {
			set defa([lindex $el 0]) [lindex $el 1]
		}
	}
	set todo {}
	set optargs {}
	foreach el $arg {
		if {[info exists field]} {
			set a($field) $el
			if {![info exists defa($field)]} {
				if {$useargs} {
					lappend optargs -$field $el
				} else {
					error "unknown option -$field, cmd should be: [oargserr $cmd $def]"
				}
			}
			unset field
		} elseif {$parseopts && [string index $el 0] eq "-"} {
			if {$el eq "--"} {
				# no more options
				set parseopts 0
			} else {
				set field [string range $el 1 end]
			}
		} else {
			lappend todo $el
		}
	}
	if {[info exists field]} {
		error "option -$field without value, should be: [oargserr $cmd $def]"
	}
	set apos 0
	set len [llength $todo]
	foreach el $def {
		set field [lindex $el 0]
		if {[info exists a($field)]} {
			uplevel [list set $field $a($field)]
		} elseif {$apos < $len} {
			uplevel [list set $field [lindex $todo $apos]]
			incr apos
		} elseif {[info exists defa($field)]} {
			uplevel [list set $field $defa($field)]
		} else {
			error "missing arg(s): $field, should be: [oargserr $cmd $def]"
		}
	}
	set todo [lrange $todo $apos end]
	lappend todo {*}$optargs
	if {[llength $todo] && !$useargs} {
		error "too many args: [join $todo ,], should be: [oargserr $cmd $def]"
	}
	uplevel [list set args $todo]
}

proc file_resolve {file {lastlinkVar {}}} {
	if {$lastlinkVar ne ""} {upvar $lastlinkVar lastlink}
	if {$::tcl_platform(platform) eq "unix"} {
		set file [file normalize $file]
		while 1 {
			if {[catch {set link [file readlink $file]}]} break
			if {[file pathtype $link] ne "absolute"} {set link [file normalize [file join [file dir $file] $link]]}
			set lastlink $file
			set file [file normalize $link]
		}
	}
	return $file
}

proc multimatch {patterns value} {
	set final 1
	foreach temp $patterns {
		if {![string match $temp $value]} {
			set final 0; break
		}
	}
	return $final
}
