proc bcol_indexlines {file indexfile} {
	set time [file mtime $file]
	set ext [file extension $file]
	if {[inlist {.rz .bgz .gz} $ext]} {set compressed 1} else {set compressed 0}
	if {![file exists $indexfile] || [file mtime $indexfile] < $time} {
		file mkdir [file dir $indexfile]
		if {$compressed} {
			progress start 1 "uncompressing $file for indexing, please be patient"
			progress message "uncompressing $file for indexing, please be patient (no progress shown)"
			if {$ext eq ".rz"} {
				set indexdir [gzroot $file].index
				set tempfile $indexdir/[file root [file tail $file]]
				if {![file exists $tempfile]} {
					gunzip $file $tempfile
				}
			} else {
				gunzip $file
				set file [file root $file]
				set tempfile $file
			}
			progress stop
		} else {
			set tempfile $file
		}
		set f [open $tempfile]
		set header [tsv_open $f comment]
		close $f
		if {[catch {tsv_basicfields $header 4} poss]} {
			if {![catch {tsv_basicfields $header 3} poss]} {
				lappend poss -1
			} else {
				set poss {-1 -1 -1 -1}
			}
		}
		if {![file exists $indexfile] || [file mtime $indexfile] < $time} {
			progress start [file size $tempfile] "Indexing $file, please be patient"
			progress message "Indexing $file, please be patient"
			exec bcol_indexfile $tempfile $indexfile.temp $indexfile.bin.temp {*}$poss
			file rename -force $indexfile.bin.temp $indexfile.bin
			file rename -force $indexfile.temp $indexfile
			progress stop
		}
		if {$compressed} {
			if {$ext eq ".rz"} {
				file delete $tempfile
			} else {
				cg razip $file
				set file $file.rz
			}
		}
	}
}

proc bcol_open indexfile {
	set result [dict create objtype bcol file $indexfile binfile $indexfile.bin]
	set f [open $indexfile]
	set header [tsv_open $f comment]
	set comment [split $comment \n]
	if {[lindex $comment 0] ne "# binary column"} {error "file \"$indexfile\" is not a binary column file"}
	if {$header ne "begin type offset"} {error "file \"$indexfile\" has an incorrect table header"}
	foreach line $comment {
		set line [string range $line 1 end]
		dict set result [lindex $line 0] [lindex $line 1]
	}
	set table {}
	while {![eof $f]} {
		set line [gets $f]
		if {![llength $line]} continue
		lappend table [split $line \t]
	}
	close $f
	set tabled [dict create]
	list_foreach {num type offset} $table {
		if {$type eq "end"} break
		dict set tabled $num $offset
	}
	dict set result table $table
	dict set result tabled $tabled
	dict set result max [lindex $table end 0]
	set fi [open $indexfile.bin]
	fconfigure $fi -encoding binary -translation binary
	dict set result fi $fi
	return $result
}

proc bcol_size bcol {
	expr {[dict get $bcol max]+1}
}

proc bcol_close bcol {
	catch {close [dict get $bcol fi]}
}

proc bcol_get {bcol start end} {
	if {[dict get $bcol objtype] ne "bcol"} {error "This is not a bcol object: $bcol"}
	set type [dict get $bcol type]
	set max [dict get $bcol max]
 	set table [dict get $bcol table]
	set tabled [dict get $bcol tabled]
	if {$start > $max} {
		return {}
	}
	if {$end > $max} {set end $max}
	set len [expr {$end-$start+1}]
	if {$len == 0} {return {}}
	set offset 0
	if {[llength $table] > 2} {
		set offset 0
		list_foreach {num type noffset} $table {
			if {$start < $num} {
				break
			}
			set offset $noffset
		}
	}
	set pos [expr {4*$start}]
	set binfile [dict get $bcol binfile]
	set f [open $binfile]
	fconfigure $f -encoding binary -translation binary
	seek $f $pos
	set b [read $f [expr {4*$len}]]
	binary scan $b i$len sresult
	close $f
	set result {}
	foreach v $sresult {
		catch {set offset [dict get $tabled $start]}
		incr start
		lappend result [expr {($v & 0xffffffff) + $offset}]
	}
	return $result
}

proc cg_bcol {cmd args} {
	switch $cmd {
		get {
			foreach indexfile $args break
			set bcol [bcol_open $indexfile]
			set result [bcol_get $bcol {*}[lrange $args 1 end]]
			bcol_close $bcol
			puts $result
			return $result
		}
		size {
			foreach indexfile $args break
			set bcol [bcol_open $indexfile]
			set size [bcol_size $bcol]
			bcol_close $bcol
			puts $size
			return $size
		}
		default {
			error "unknown subcommand $cmd of bcol, must be one of: get, size"
		}
	}
}

proc cg_index file {
	set time [file mtime $file]
	set indexdir [gzroot $file].index
	set ext [file extension $file]
	if {[inlist {.rz .bgz .gz} $ext]} {set compressed 1} else {set compressed 0}
	file mkdir $indexdir
	set indexfile $indexdir/lines.bcol
	bcol_indexlines $file $indexfile
	if {![file exists $indexdir/info.tsv] || [file mtime $indexdir/info.tsv] < $time} {
		set f [gzopen $file]
		set header [tsv_open $f]
		catch {close $f}
		set bcol [bcol_open $indexfile]
		set size [bcol_size $bcol]
		bcol_close $bcol
		set f [open $indexdir/info.tsv w]
		puts $f key\tvalue
		puts $f file\t$file
		puts $f lineindexfile\t[file tail $indexfile]
		puts $f header\t$header
		puts $f size\t$size
		close $f
	}
	return $indexfile
}

proc cg_size file {
	set indexfile [cg_index $file]
	set bcol [bcol_open $indexfile]
	set size [bcol_size $bcol]
	bcol_close $bcol
	puts $size
	return $size
}
