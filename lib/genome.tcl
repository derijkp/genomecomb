#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_fas2ifas {srcfile destfile} {
	# make new one
	set f [gzopen $srcfile]
	set o [open $destfile.temp w]
	set line [gets $f]
	set ids {}
	puts $o $line
	set curid [string range $line 1 end]
	lappend ids $curid
	set a(s,$curid) [tell $o]
	while {![eof $f]} {
		set line [gets $f]
		if {$line eq ""} continue
		if {[string index $line 0] eq ">"} {
			putslog "chromosome $curid finished"
			set a(e,$curid) [tell $o]
			puts $o \n$line
			set curid [string range $line 1 end]
			lappend ids $curid
			set a(s,$curid) [tell $o]
		} else {
			puts -nonewline $o $line
		}
	}
	puts $o ""
	putslog "chromosome $curid finished"
	set a(e,$curid) [tell $o]
	close $o
	gzclose $f
	# sort
	set sids [ssort -natural $ids]
	if {$sids ne $ids} {
		set f [open $destfile.temp]
		set o [open $destfile.temp2 w]
		set oi [open $destfile.index.temp2 w]
		foreach name $sids {
			if {![regexp {chromosome ([^ ,]+)[ ,]} $name temp chr]} {
				if {![regexp {chr([^ ,]+)} $name temp chr]} {
					set chr [lindex $name end]
				}
			}
			set len [expr {$a(e,$name) - $a(s,$name)}]
			puts $o "\>$name"
			puts $oi "$chr\t[tell $o] $len"
			seek $f $a(s,$name)
			fcopy $f $o -size [expr {$len+1}]
		}
 		close $oi
		close $o
		close $f
		file rename -force $destfile.temp2 $destfile
		file rename -force $destfile.index.temp2 $destfile.index
		file delete $destfile.temp
	} else {
		set oi [open $destfile.index.temp w]
		foreach name $sids {
			if {![regexp {chromosome ([^ ,]+)[ ,]} $name temp chr]} {
				if {![regexp {chr([^ ,]+)} $name temp chr]} {
					set chr [lindex $name end]
				}
			}
			set len [expr {$a(e,$name) - $a(s,$name)}]
			puts $oi "$chr\t$a(s,$name) $len"
		}
 		close $oi
		file rename $destfile.temp $destfile
		file rename $destfile.index.temp $destfile.index
	}
}

# index from stdin stream
proc cg_genome_indexfasta {resultfile} {
	set f stdin
	set o [open $resultfile w]
	set oi [open $resultfile.index w]
	set line [gets $f]
	set pos 0
	while {![eof $f]} {
		set len [string length $line]
		if {$len == 0} continue
		set name [string range $line 1 end]
		if {![regexp {chromosome ([^ ,]+)[ ,]} $name temp chr]} {
			if {![regexp {chr([^ ,]+)} $name temp chr]} {
				set chr [lindex $name end]
			}
		}
		puts $name
		set line ">chr$chr $name"
		puts $o $line
		incr pos [string length $line]
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
		puts $oi "$chr\t$pos $seqlen"
		flush $oi
		flush $o
		incr pos $seqlen
		incr pos 1
	}
	puts $oi "\t$pos 0"
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
	if {[file isdir $file]} {
		set file [lindex [glob $file/genome_*.ifas] 0]
	}
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

proc genome_get {f chr start end {correctend 0}} {
	global genomefasta
	if {$end < $start} {error "end ($end) is smaller than start ($start)"}
	set fastaindex $genomefasta($f)
	if {[catch {
		set temp [dict get $fastaindex $chr]
	}]} {
		set chr [chr_clip $chr]
		set temp [dict get $fastaindex $chr]
	}
	foreach {gstart glen} $temp break
	if {$end > $glen} {
		if {!$correctend} {
			error "trying to get sequence beyond end of chromosome ($end > $glen)"
		} else {
			set end $glen
		}
	}
	set pos [expr {$gstart+$start}]
	seek $f $pos
	read $f [expr {$end-$start}]
}

proc genome_mask {dbdir seq chr estart eend {freql 0} {freqN 0.2} {delsize 5} {repeats s} {snpdbpatterns snp}} {
	global dbsnpdata
	# init dbsnpdata
	if {![info exists dbsnpdata($dbdir)]} {
		set dbsnpfiles {}
		foreach snpdbpattern $snpdbpatterns {
			lappend dbsnpfiles {*}[gzfiles $dbdir/var_*$snpdbpattern*.tsv.gz]
		}
		set dbsnpposs {}
		foreach dbsnp $dbsnpfiles {
			set dbsnpheader [cg select -h $dbsnp]
			set temp [tsv_basicfields $dbsnpheader 4]
			lappend temp [lsearch $dbsnpheader freq]
			lappend dbsnpposs $temp
		}
		set dbsnpdata($dbdir) [dict create dbsnpfiles $dbsnpfiles dbsnpposs $dbsnpposs]
	} else {
		set dbsnpfiles [dict get $dbsnpdata($dbdir) dbsnpfiles]
		set dbsnpposs [dict get $dbsnpdata($dbdir) dbsnpposs]
	}
	# repeats are already masked, change lowercase to N if hardmasking is required
	if {$repeats eq "0"} {
		set seq [string toupper $seq]
	} elseif {$repeats eq "N"} {
		regsub -all {[a-z]} $seq N seq
	}
	# mask snps
	set list {}
	foreach snpposs $dbsnpposs dbsnp $dbsnpfiles {
		set temp [tabix $dbsnp chr$chr $estart $eend]
		lappend list {*}[list_subindex $temp $snpposs]
	}
	set list [ssort -natural -decreasing $list]
	list_foreach {c s e type ffreq} $list {
		if {$e <= $estart || $s > $eend} continue
		set freq 0
		foreach el [split $ffreq ,] {
			if {[isdouble $el] && $el > $freq} {set freq $el}
		}
		if {$freq <= $freql} continue
		set start [expr {$s-$estart}]
		if {$start < 0} {set start 0}
		if {$type eq "ins"} {
			set end $start
		} else {
			set end [expr {$e-$estart-1}]
			if {$end < $start} {set end $start}
		}
		if {$type eq "del" && ($delsize != -1) && ([expr {$end-$start}] > $delsize)} continue
		set base [string range $seq $start $end]
		if {$freq > $freqN} {
			regsub -all . $base N base
		} else {
			set base [string tolower $base]
		}
		set seq [string_replace $seq $start $end $base]
	}
	return $seq
}

proc cg_make_genomecindex {ifasfile} {
	package require cindex
	set ifasfile [file_absolute $ifasfile]
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
		set result $base.ssa/[file tail $base]-$chr
		if {![file exists $result.ssa]} {
			puts "creating index $chr"
			time {
				set seq [gets $f]
				set seq [string toupper $seq]
			}
			set o [cindex create $seq]
			puts "saving index $result"
			cindex save $o $result.temp
			file rename -force $result.temp.ssa $result.ssa
		} else {
			puts "skipping $chr: already made"
			time {gets $f}
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
		error "format is: genome_get genome chromosome ?start end?"
	}
}

