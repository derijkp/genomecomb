#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_downloadgenome {build result} {
	set keepdir [pwd]
	set result [file normalize $result]
	set files {}
	file mkdir $result.temp
	cd $result.temp
	set files {}
	foreach chr {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 M X Y} {
		lappend files chr$chr.fa.gz
		if {[file exists chr$chr.fa.gz]} {
			puts "Skipping chromosome $chr: chr$chr.fa.gz already there"
			continue
		}
		putslog "Downloading chromosome chr$chr.fa.gz"
		catch {exec wget --tries=45 ftp://hgdownload.cse.ucsc.edu/goldenPath/$build/chromosomes/chr$chr.fa.gz >@ stdout  2>@ stderr}
	}
	if {![file exists $result]} {
		putslog "Converting and indexing"
		exec zcat {*}$files | cg genome_indexfasta [file tail $result]
		file rename {*}[glob [file tail $result]*] ..
	}
	set rfile [file dir $result]/reg_[file root [file tail $result]].tsv
	putslog "Making $rfile"
	cd ..
	set data [file_read $result.index]
	set o [open $rfile w]
	puts $o chromosome\tbegin\tend
	list_foreach {chromosome begin len} [split [string trim $data] \n] {
		puts $o $chromosome\t0\t$len
	}
	close $o
	cd $keepdir
}

proc cg_genome_indexfasta {resultfile} {
	set f stdin
	set o [open $resultfile w]
	set oi [open $resultfile.index w]
	set line [gets $f]
	set pos 0
	while {![eof $f]} {
		set len [string length $line]
		set name [string range $line 1 end]
		if {![regexp {chromosome ([^ ,]+)[ ,]} $name temp chr]} {
			if {![regexp {chr([^ ,]+)} $name temp chr]} {
				set chr [lindex $name end]
			}
		}
		puts $name
		set line ">chr$chr $name"
		puts $o $line
		incr pos $len
		incr pos 1
		set seqlen 0
		while {![eof $f]} {
			set line [gets $f]
			if {[string index $line 0] eq ">"} {
				break
			}
			set len [string length $line]
			incr seqlen $len
			puts -nonewline $o $line
		}
		puts $o ""
		set a($name) [list $pos $seqlen]
		puts $oi "$name\t$pos $seqlen"
		flush $oi
		flush $o
		incr pos $seqlen
		incr pos 1
	}
	close $oi
	close $o
	close $f
}

proc genome_makefastaindex {fastafile} {
	putslog "Making index for $fastafile"
	set f [open $fastafile]
	set oi [open $fastafile.index.temp w]
	while {![eof $f]} {
		set line [gets $f]
		set name [string range $line 1 end]
		if {![regexp {chromosome ([^ ,]+)[ ,]} $name temp chr]} {
			if {![regexp {chr([^ ,]+)} $name temp chr]} {
				set chr [lindex $name end]
			}
		}
		set start [tell $f]
		set seq [gets $f]
		set seqlen [string length $seq]
		putslog $name
		puts $oi "$chr\t$start $seqlen"
	}
	close $oi
	close $f
	file rename -force $fastafile.index.temp $fastafile.index
}

proc genome_open {file} {
	global genomefasta
	set f [open $file]
	set fastaindex {}
	if {![file exists $file.index]} {
		genome_makefastaindex $file
	}
	foreach {name data} [split [string trim [file_read $file.index]] \t\n] {
		if {![regexp {chr([^ ,]+)} $name temp chr]} {set chr [lindex $name end]}
		dict set fastaindex $chr $data
	}
	set genomefasta($f) $fastaindex
	return $f
}

proc genome_close {f} {
	global genomefasta
	close $f
	unset genomefasta($f)
}

proc genome_get {f chr start end} {
	global genomefasta
	set fastaindex $genomefasta($f)
	if {[catch {
		set temp [dict get $fastaindex $chr]
	}]} {
		regsub ^chr $chr {} chr
		set temp [dict get $fastaindex $chr]
	}
	foreach {gstart glen} $temp break
	set pos [expr {$gstart+$start}]
	seek $f $pos
	read $f [expr {$end-$start}]
}

proc cg_make_genomecindex {ifasfile} {
	package require cindex
	set ifasfile [file normalize $ifasfile]
	set f [open $ifasfile]
	set base [file root $ifasfile]
	file mkdir $base.ssa
	while {![eof $f]} {
		set name [gets $f]
		puts $name
		if {![string length $name]} break
		if {![regexp {chromosome ([^ ,]+)[ ,]} $name temp chr]} {
			if {![regexp {chr([^ ,]+)} $name temp chr]} {
				set chr [lindex $name end]
			}
		}
		time {
			set seq [gets $f]
			set seq [string toupper $seq]
		}
		set result $base.ssa/[file tail $base]-$chr
		if {![file exists $result.ssa]} {
			puts "creating index $chr"
			set o [cindex create $seq]
			puts "saving index $result"
			cindex save $o $result.temp
			file rename $result.temp.ssa $result.ssa
		}
	}
	close $f
}

proc cg_genome_get {args} {
	foreach {genome chr start end} $args break
	if {[llength $args] == 4} {
		set f [genome_open $genome]
		puts [genome_get $f $chr $start $end]
		close $f
	} elseif {[llength $args] == 2} {
		global genomefasta
		set f [genome_open $genome]
		set fastaindex $genomefasta($f)
		foreach {gstart glen} [dict get $fastaindex $chr] break
		seek $f $gstart
		fcopy $f stdout -size $glen
		close $f
	} else {
		error "format is: genome_get chromosome ?start end?"
	}
}

