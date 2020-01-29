#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc t {} {
	global threads
	set num 0
	foreach thread $threads {
		puts "*$num* [list_subindex $thread {0 1 2 4 5}]"
		incr num
	}
}

proc oformat {value {num 2}} {
	if {[catch {format %0.${num}f $value} result]} {
		return $value
	}
	return $result
}

proc getmategap {libfile} {
	set f [opencgifile $libfile header]
	while {![eof $f]} {
		set line [gets $f]
		foreach {id type armId indArm objArm min max} $line break
		if {$type eq "mategap"} {return $min}
	}
	return {}
}

proc map2sv {files prefix} {
	global appdir
	if {[llength $files] == 1 && [file isdir [lindex $files 0]]} {
		set cgdir [lindex $files 0]
		set files [glob -nocomplain $cgdir/MAP/*/mapping_*.tsv.bz2]
		if {![llength $files]} {
			set files [glob -nocomplain $cgdir/mapping_*.tsv.bz2]
		}
		if {![llength $files]} {
			error "No files found"
		}
	}

	set rdir [file dir $prefix]
	file mkdir $rdir
	if {[file exists $rdir/working_in_scratch.txt]} {
		set scratchdir [file_read $rdir/working_in_scratch.txt]
		file mkdir $scratchdir
	} else {
		set scratchdir [scratchdir]
		file_write $rdir/working_in_scratch.txt $scratchdir
	}
#	set rem [glob -nocomplain $scratchdir/*]
#	if {[llength $rem]} {
#		file delete {*}$rem
#	}
	if {![file exists ${prefix}_map2sv_sort_FINISHED]} {
		file mkdir $scratchdir/tmp
		set tail [file tail $prefix]
		set scratchprefix $scratchdir/tmp/$tail-
		set num 0
		foreach file $files {
			set root [file root [file root [file tail $file]]]
			if {[llength [glob -nocomplain $scratchdir/$tail-$root-*]]} {
				puts "already done: $file"
			} else {
				puts "doing: $file"
				set cat [gzcat $file]
				exec {*}$cat $file | $appdir/bin/map2sv $num | $appdir/bin/distr2chr $scratchprefix$root-
				file rename -force -- {*}[glob $scratchprefix$root-*] $scratchdir
			}
			incr num
		}
		set files [glob $scratchdir/$tail-*]
		unset -nocomplain a
		foreach file $files {
			set tail [file tail $file]
			set chr [lindex [split $tail -] end]
			lappend a($chr) $file
		}
		foreach chr [bsort [array names a]] {
			set rfile $prefix-$chr-paired.tsv
			if {[file exists $rfile]} {
				puts "$rfile exists"
				continue
			} else {
				puts "Making $rfile"
			}
#			file mkdir [file dir $file]/tmp
#			if {[file extension $file] eq ".tsv"} continue
#			if {[file extension $file] eq ".rz"} continue
#			if {[file extension $file] eq ".bgz"} continue
			set f [open $rfile.temp w]
			puts $f [join {chromosome bin strand1 start1 end1 weight1 numl type chr2 strand2 start2 end2 weight2 numr dist num fnum side} \t]
			close $f
			exec cat {*}$a($chr) | gnusort8 -T $scratchdir/tmp -t \t -n -s -k5 >> $rfile.temp
			file rename -force -- $rfile.temp $rfile
		}
		file delete -force $scratchdir/tmp {*}$files
	}
	exec touch ${prefix}_map2sv_sort_FINISHED
}

proc bam2sv {bamfile prefix} {
	global appdir
	set bamfile [file_absolute $bamfile]
	set prefix [file_absolute $prefix]
	file mkdir [file dir $prefix]
	if {![file exists ${prefix}_map2sv_sort_FINISHED]} {

		foreach chr [array names o] {
			close $o($chr)
		}
		catch {close $f}
		set f [open "| [list samtools view $bamfile]"]
		# numl,numr,num,bin not directly useful here
		set numl ""; set numr ""; set num "" ; set bin ""; set side ""
		unset -nocomplain todo
		unset -nocomplain o
		set curchr ""
		set fnum -1
		while {![eof $f]} {
			set line [split [gets $f] \t]
			incr fnum
			if {![llength $line]} continue
			foreach {qname flags} $line break
			set chr2 [lindex $line 6]
			if {[sam_unmapped $flags]} continue
			if {[info exists todo($qname)] || $chr2 eq "*" || [sam_unmapped_mate $flags]} {
				if {[info exists todo($qname)]} {
					set line2 $line
					set line1 $todo($qname)
					unset todo($qname)
				} else {
					unset -nocomplain line2
					set line1 $line
				}
				foreach {qname1 flags1 chr1 pos1 mapq1 cigar1 chr2 pos2 tlen1 seq1 qual1} $line1 break
				if {![info exists o($chr1)]} {
					set o($chr1) [open $prefix-$chr1-paired.tsv.temp w]
					puts $o($chr1) [join {chromosome bin strand1 start1 end1 weight1 numl type chr2 strand2 start2 end2 weight2 numr dist num fnum side} \t]
				}
				incr pos1 -1
				# this should take into account cigar, and is thus not correct
				# but need to get it to work good enough fast ...
				set end1 [expr {$pos1+[string length $seq1]}]
				set weight1 $mapq1
				if {[sam_reverse $flags1]} {set strand1 -} else {set strand1 +}
				if {[info exists line2]} {
					foreach {qname2 flags2 chr2 pos2 mapq2 cigar2 temp temp temp seq2 qual2} $line2 break
					if {$qname1 ne $qname2} {error "names of pair do not match:\n$line1\n$line2"}
					incr pos2 -1
					# this should take into account cigar, and is thus not correct
					# but need to get it to work good enough fast ...
					set end2 [expr {$pos2+[string length $seq2]}]
					set weight2 $mapq2
					if {[sam_reverse_mate $flags1]} {set strand2 -} else {set strand2 +}
					set dist [expr {$pos2-$end1}]
					if {$chr1 ne $chr2} {
						set type c
						set dist -1
					} elseif {$strand1 eq $strand2} {
						set type r
					} else {
						set type u
					}
				} else {
					set dist ""
					set type s
					foreach {chr2 strand2 pos2 end2 weight2 numr dist} {{} {} {} {} {} {} {}} break
				}
				puts $o($chr1) [join [list $chr1 $bin $strand1 $pos1 $end1 $weight1 $numl $type $chr2 $strand2 $pos2 $end2 $weight2 $numr $dist $num $fnum $side] \t]
			} else {
				set todo($qname) $line
			}
		}

		foreach chr [array names o] {
			close $o($chr)
			cg select -s {chromosome start1 end1 chr2 start2 end2} $prefix-$chr1-paired.tsv.temp $prefix-$chr1-paired.tsv.temp2
			file delete $prefix-$chr1-paired.tsv.temp
			file rename -force -- $prefix-$chr1-paired.tsv.temp2 $prefix-$chr1-paired.tsv
		}
		catch {close $f}
	}
	exec touch ${prefix}_map2sv_sort_FINISHED
}

proc sv2db {files} {
	foreach file $files {
		puts $file
		set dbfile $file.sqlite
		if {[file exists $dbfile]} {
			puts "skipping $file: exists"
			continue
		}
		exec sqlite3 $dbfile {
			create table hits (
				chr1 integer,
				bin integer,
				strand1 char,
				start1 integer,
				end1 integer,
				weight1 integer,
				numl integer,
				type char,
				chr2 integer,
				strand2 integer,
				start2 integer,
				end2 integer,
				weight2 integer,
				numr integer,
				dist integer,
				num integer,
				fnum integer,
				side integer)
		}
		exec sqlite3 -separator \t $dbfile ".import \"[file_absolute $file]\" hits"
		exec sqlite3 $dbfile "create index hits_type on hits(type)"
		exec sqlite3 $dbfile "create index hits_bin on hits(bin)"
		svdbinfo $dbfile
	}
}

proc svinfo {pairfile} {
	puts "svinfo on file $pairfile"
	svhisto $pairfile
	set list [list_remove [split [file_read [file root [gzroot $pairfile]].disthisto] \n] {}]
	list_shift list
	list_pop list
	set total [lmath_sum [list_subindex $list 1]]
	set slist [lsort -integer -index 1 $list]
	set mode [lindex $slist end 0]
	if {$mode eq ""} return
	if {$mode == -1} {
		set mode [lindex $slist end-1 0]
	}
	set minnum [expr {round(0.05*$total)}]
	set min10num [expr {round(0.10*$total)}]
	set min25num [expr {round(0.25*$total)}]
	set curnum 0
	set min 0
	set min10 0
	list_foreach {gapsize num} $list {
		if {!$min && $curnum > $minnum} {
			set min $gapsize
		}
		if {!$min10 && $curnum > $minnum} {
			set min10 $gapsize
		}
		if {$curnum > $min25num} {
			set min25 $gapsize
			break
		}
		incr curnum $num
	}
	set curnum 0
	set max 0
	set max10 0
	list_foreach {gapsize num} [list_reverse $list] {
		if {!$max && $curnum > $minnum} {
			set max $gapsize
		}
		if {!$max10 && $curnum > $minnum} {
			set max10 $gapsize
		}
		if {$curnum > $min25num} {
			set max25 $gapsize
			break
		}
		incr curnum $num
	}
	set o [open [gzroot $pairfile].numinfo.temp w]
	puts $o key\tvalue
	puts $o mode\t$mode
	puts $o min\t$min
	puts $o max\t$max
	puts $o min10\t$min10
	puts $o max10\t$max10
	puts $o min25\t$min25
	puts $o max25\t$max25
	close $o
	file rename -force -- [gzroot $pairfile].numinfo.temp [gzroot $pairfile].numinfo
	putslog "Made svinfo $pairfile.numinfo"
}

proc svhisto {pairfile} {
	set out [file root [gzroot $pairfile]].disthisto
	set f [gzopen $pairfile]
	set header [gets $f]
	set distpos [lsearch $header dist]
	set num 1
	while {![eof $f]} {
		incr num
		if {![expr $num%1000000]} {putsprogress $num}
		set dist [lindex [split [gets $f] \t] $distpos]
		if {![info exists a($dist)]} {
			set a($dist) 1
		} else {
 			incr a($dist)
		}
	}
	unset -nocomplain a()
	unset -nocomplain a(-)
	set f [open $out.temp w]
	puts $f "dist\tnumber"
	set total 0
	foreach dist [lsort -integer [array names a]] {
		puts $f $dist\t$a($dist)
		incr total $a($dist)
	}
	puts $f $total
	close $f
	catch {file rename -force -- $out $out.old}
	file rename -force -- $out.temp $out
	# draw histo
	set tempfile [tempfile]
	file_write $tempfile [subst -nocommands {
		#set terminal postscript eps enhanced color
		set terminal png
		set output "$out.png"
		set ylabel "Number of reads"
		set xlabel "distance between paired ends"
		plot [0:600] "$out" using 1:2
	}]
	catch {exec gnuplot $tempfile}
}

proc table2tsv {table} {
	set result {}
	foreach line $table {
		append result [join $line \t]\n
	}
	return $result
}

proc sv_maxima {list {min 0}} {
	set maxima {}
	set pos -1
	set maxnum 0
	set maxpos 0
	set cutoff 0
	set rise 1
	foreach num $list {
		if {$rise} {
			if {$num > $maxnum} {
				set maxnum $num
				set maxpos $pos
				set cutoff [expr {$maxnum - 0.01}]
			} elseif {$num < $cutoff} {
				if {$maxnum > $min} {
					lappend maxima [list $maxpos $maxnum]
				}
				set rise 0
				set minpos $pos
				set minnum $num
			}
		} else {
			if {$num < $minnum} {
				set minnum $num
				set minpos $pos
				set cutoff [expr {$minnum + 0.01}]
			} elseif {$num > $cutoff} {
				lappend minima [list $minpos $minnum]
				set rise 1
				set maxpos $pos
				set maxnum $num
			}
		}
		incr pos
	}
	if {[llength $maxima] == 1} {
		set line [lindex $maxima 0]
		lappend line [lindex $line 0]
		return [list $line]
	}
	# join out close together peaks of similar height (artefacts)
	set result {}
	set pos -1
	set prev [lindex $maxima 0]
	foreach line [lrange $maxima 1 end] {
		incr pos
		if {![llength $prev]} {
			set prev $line
			continue
		}
		foreach {p1 h1} $prev break
		foreach {p2 h2} $line break
		set diff [expr {$p2-$p1}]
		set hdiff [expr {abs($h2-$h1)}]
		if {($diff < 30) && ($hdiff < 5)} {
			set mh 0
			list_foreach {pos h} $minima {
				if {($pos >= $p1) && ($pos <= $p2)} {
					set mh $h
					break
				}
			}
			set mdiff [expr {max($h2,$h1)-$mh}]
			if {$mdiff < 5} {
				lappend prev $p2
				lappend result $prev
				set prev {}
			} else {
				lappend prev $p1
				lappend result $prev
				set prev $line
			}
		} else {
			lappend prev $p1
			lappend result $prev
			set prev $line
		}
	}
	if {[llength $prev]} {
		lappend prev [lindex $prev 0]
		lappend result $prev
	}
	return $result
}

proc svtools_loadindex {file xfield} {
	global svi
	set f [gzopen $file]
	set svi(header) [tsv_open $f]
	catch {close $f}
	set svi(xfield) $xfield
	set svi(xfieldpos) [lsearch $svi(header) $xfield]
	set indexname [gzroot $file].${xfield}_index
	set o [open $indexname]
	set svi(step) [gets $o]
	set svi(findex) [gets $o]
	set xmin [gets $o]
	set svi(xmin) $xmin
	set svi(xmax) [gets $o]
	set svi(fx) [expr {$svi(xmin)-$xmin%10000}]
	set temp [split [string trim [read $o]] \n]
	set svi(index) $temp
	close $o
}

proc svtools_aprgoto {pairfile start} {
	global svi
	svtools_loadindex $pairfile end1
	set xfieldpos $svi(xfieldpos)
	if {$start < $svi(findex)} {set start $svi(findex)}
	set index $svi(index)
	set fpos [expr {round([lindex $index [expr {($start-$svi(findex))/10000}]])}]
	set f [gzopen $pairfile $fpos]
	while {![eof $f]} {
		set fpos [tell $f]
		set line [gets $f]
		set pos [lindex [split $line \t] $xfieldpos]
		if {$pos > $start} break
	}
	return $f
}

proc e {table} {
	[edit].editor set [table2tsv $table]
}

proc pf {list} {
	set temp {}
	foreach el $list {
		lappend temp [format %0.1f $el]
	}
	return $temp
}

proc maxpos {vector} {
	set maxval [lmath_max $vector]
	set maxpos [lsearch $vector $maxval]
}

proc svdiffscore {diff} {
	if {$diff < 10} {
		return 2
	} elseif {$diff < 15} {
		return 1
	} elseif {$diff > 100} {
		return -1000
	} else {
		return [expr {-($diff-10)/10}]
	}
}

proc cluster_fts {data {checkpos {}} {maxdist 20}} {
	global infoa
	set checknum 20
	set cutoffdist 150
	set match 3
	set mismatch -5
	set mismatchlong -1
	set maxdist2 [expr {2*$maxdist}]
	set result {}
	if {$checkpos eq ""} {
		set checkpos $infoa(end1pos)
	}
	while 1 {
		if {[llength $data] <= 7} break
		if {[llength $data] <= $checknum} {
			set checknum 5
		}
		set xs [list_subindex $data $checkpos]
		set ds [lmath_calc [lrange $xs 1 end] - [lrange $xs 0 end-1]]
		# search for region with highest number of smalldists
		set pos 0
		set numsmall 0
		foreach d [lrange $ds 0 [expr {$checknum-1}]] {
			if {$d < $maxdist} {incr numsmall}
		}
		set maxnumsmall $numsmall
		set rpos 0
		set cpos $checknum
		set startpos 0
		foreach d [lrange $ds $checknum end] {
			if {$d < $maxdist} {incr numsmall}
			if {[lindex $ds $rpos] < $maxdist} {incr numsmall -1}
			if {$numsmall > $maxnumsmall} {
				set maxnumsmall $numsmall
				set startpos $cpos
			}
			incr rpos
			incr cpos
		}
		set percent [expr {double($maxnumsmall)/$checknum}]
		if {$percent < 0.5} break
		# do not start at bad one
		while 1 {
			set d [lindex $ds $startpos]
			if {![isint $d]} break
			if {$d < $maxdist} break
			incr startpos
		}
		if {![isint $d]} break
		# find optimum before match
		set score [expr {$maxnumsmall * $match}]
		set limit $mismatch
		set maxscore $score
		set maxpos $startpos
		set pos $startpos
		while {$pos >= 0} {
			set diff [lindex $ds $pos]
			if {$diff > $cutoffdist} break
			incr score [svdiffscore $diff]
			if {$score >= $maxscore} {
				set maxscore $score
				set maxpos $pos
				if {$score > 20} {
					set limit [expr {round($maxscore/4)}]
				}
			}
			if {$score <= $limit} break
			incr pos -1
		}
		set rstart $maxpos
		# find optimum after match
		set score [expr {$maxnumsmall * $match}]
		set limit $mismatch
		set maxscore $score
		set pos $startpos
		set range [lrange $ds $pos end]
		foreach diff $range {
			if {$diff > $cutoffdist} break
			incr score [svdiffscore $diff]
			if {$score >= $maxscore} {
				set maxscore $score
				set maxpos $pos
				if {$score > 20} {
					set limit [expr {round($maxscore/4)}]
				}
			}
			if {$score <= $limit} break
			incr pos
		}
		incr maxpos
		set rend $maxpos
		# region found
		set list [lrange $data $rstart $rend]
		lappend result $list
		incr rstart -1
		incr rend
		set data [list_concat [lrange $data 0 $rstart] [lrange $data $rend end]]
	}
	return $result
}

proc rkde {list} {
	# set list [list_subindex $table 11]
	set tempfile [tempfile]
	set outfile [tempfile]
	file_write $tempfile [join $list \n]
	set r [open "|R --slave" w]
	puts $r [subst {x = read.table("$tempfile")}]
	puts $r {
		p = x$V1
		min = min(p)
		max = max(p)
		d = density(p,bw=10,kernel="g",from=min,to=max,n=(max-min+1))
		md = max(d$y)
		t = data.frame(d$x,d$y)
	}
	puts $r [subst {write.table(t,file="$outfile")}]
	close $r
	set r [open $outfile]
	set c [csv_file $r " "]
	close $r
	list_shift c
	return [list_subindex $c {1 2}]
}

proc kde {list {kernel {}}} {
	set list [lsort -real $list]
	set s [lindex $list 0 0]
	set e [lindex $list end 0]
	if {$kernel eq ""} {
		# gausian kernel, bandwidth = 10
		# set kernel {0.0060 0.0079 0.0104 0.0136 0.0175 0.0224 0.0283 0.0355 0.0440 0.0540 0.0656 0.0790 0.0940 0.1109 0.1295 0.1497 0.1714 0.1942 0.2179 0.2420 0.2661 0.2897 0.3123 0.3332 0.3521 0.3683 0.3814 0.3910 0.3970 0.3989 0.3970 0.3910 0.3814 0.3683 0.3521 0.3332 0.3123 0.2897 0.2661 0.2420 0.2179 0.1942 0.1714 0.1497 0.1295 0.1109 0.0940 0.0790 0.0656 0.0540 0.0440 0.0355 0.0283 0.0224 0.0175 0.0136 0.0104 0.0079 0.0060}
		set kernel $::infoa(kernel)
	}
	set ke [llength $kernel]
	# join [list_subindex $list {0 1 2}] \n
	set plen [expr {$e-$s}]
	set len [expr {$plen+2*$ke}]
	set data [list_fill $len 0]
	set prev [lindex $list 0]
	set num 0
	list_foreach p $list {
		if {$p != $prev} {
			set prev $p
			set p [expr {$p-$s+$ke}]
			lset data $p $num
			set num 0
		}
		incr num
	}
	# convolving with the filter over this is the same as doing the discrete kde
	set data [lmath_filter $data $kernel]
	# draw $data $s
	# data starts at [expr {$s-$ke}]
	return [list [expr {$s-$ke}] $data]
}

proc listmaxima {data {min 0} {grace 0.01}} {
	set maxima {}
	set pos 0
	set maxnum 0
	set maxpos 0
	set cutoff 0
	set rise 1
	foreach num $data {
		if {$rise} {
			if {$num > $maxnum} {
				set maxnum $num
				set maxpos $pos
				set cutoff [expr {$maxnum - $grace}]
			} elseif {$num < $cutoff} {
				if {$maxnum > $min} {
					lappend maxima [list $maxpos $maxnum]
				}
				set rise 0
				set minpos $pos
				set minnum $num
			}
		} else {
			if {$num < $minnum} {
				set minnum $num
				set minpos $pos
				set cutoff [expr {$minnum + $grace}]
			} elseif {$num > $cutoff} {
				lappend minima [list $minpos $minnum]
				set rise 1
				set maxpos $pos
				set maxnum $num
			}
		}
		incr pos
	}
	return $maxima
}

proc simple_group {list {dist 50} {minnum 3}} {
	set prev [lindex $list 0]
	set result {}
	set curresult {}
	foreach v $list {
		set d [expr {$v-$prev}]
		if {$d <= $dist} {
			lappend curresult $v
		} else {
			if {[llength $curresult] >= $minnum} {
				lappend result $curresult
			}
			set curresult {}
		}
		set prev $v
	}
	if {[llength $curresult] >= $minnum} {
		lappend result $curresult
	}
	return $result
}

proc kde_maxima {dists {maxsize 50} {min 1.5}} {
	# returns begin and end of "peak"
	if {![llength $dists]} {return {}}
	set dists [lsort -real $dists]
	set result {}
	set groups [simple_group $dists]
	foreach group $groups {
		if {[llength $group] < 5} continue
		set s [lindex $group 0]
		set e [lindex $group end]
		set size [expr {$e-$s}]
		# simply return first and last point if range < maxsize
		if {$size < $maxsize} {
			lappend result [list $s $e]
			continue
		}
		foreach {s data} [kde $group] break
set ::kdedata $data
set ::kdes $s
#		draw $data $s
		set maxima [lsort -real -index 1 -decreasing [listmaxima $data $min]]
		if {![llength $maxima]} continue
		set stopat [expr {max(1.5,0.2*[lindex $maxima 0 1])}]
		set minsize [expr {0.6*$maxsize}]
		list_foreach {pos h} $maxima {
			set h [lindex $data $pos]
			if {$h < $stopat} break
			if {$h == 0} continue
			set pos1 [expr {$pos-1}]
			set pos2 [expr {$pos+1}]
			set h1 [lindex $data $pos1]
			set h2 [lindex $data $pos2]
			for {set i 0} {$i < $maxsize} {incr i} {
				if {$h1 < $h2} {
					incr pos2
					set h2 [lindex $data $pos2]
					if {$h2 < 0.1} break
				} else {
					incr pos1 -1
					set h1 [lindex $data $pos1]
					if {$h1 < 0.1} break
				}
			}
			set size [expr {$pos2-$pos1}]
			if {$size < $minsize} continue
			# if {[expr {$pos2-$pos1}] < 10} continue
			set s1 [expr {$s+$pos1+1}]
			set s2 [expr {$s+$pos2-1}]
			# set mp [expr {($s1+$s2)/2}]
			lappend result [list $s1 $s2]
			set ph1 $h1
			while {$pos1 > 0} {
				incr pos1 -1
				set h1 [lindex $data $pos1]
				if {$h1 > $ph1} break
				set ph1 $h1
			}
			set ph2 $h2
			while 1 {
				incr pos2
				set h2 [lindex $data $pos2]
				if {![isdouble $h2]} break
				if {$h2 > $ph2} break
				set ph2 $h2
			}
			set data [lreplace $data $pos1 $pos2 {*}[list_fill [expr {$pos2-$pos1+1}] 0]]
		}
	}
	return $result
}

proc sv_threadadd {threads pos item} {
	set t [lindex $threads $pos]
	lappend t $item
	lset threads $pos $t
	return $threads
}

proc sv_addtothreads {threads start maxima table} {
	global infoa
	set distpos $infoa(distpos)
	set mode $infoa(mode)
	# if empty region, reset main thread
	set numthreads [llength $threads]
	set mpos [lindex $threads 0 end 0]
	set lpos $mpos
	foreach thread $threads {
		set cpos [lindex $thread end 0]
		if {$cpos > $lpos} {set lpos $cpos}
	}
	if {[expr {$start-$lpos}] >= 200} {
		set prevmains [expr {$mode-50}]
		set prevmaine [expr {$mode+50}]
		set temp [list $start $prevmains $prevmaine {} $numthreads]
		set threads [sv_threadadd $threads 0 $temp]
	}
	# put possibly matching threads and maxima in list
	set matchlist {}
	set mpos -1
	set maxlist {}
	set totnum [llength $table]
	set usednum 0
	foreach l $maxima {
		incr mpos
		foreach {ds de} $l break
		set  part {}
		foreach line $table {
			set dist [lindex $line $distpos]
			if {($dist >= $ds) && ($dist <= $de)} {
				incr usednum
				lappend part $line
			}
		}
		lappend maxlist [list $start $ds $de $part]
		set pos -1
		foreach thread $threads {
			incr pos
			set last [lindex $thread end]
			foreach {temp tds tde} $last break
			set overlap [overlap $ds $de $tds $tde]
			set minoverlap [expr {min($de-$ds,$tde-$tds)/3}]
			if {$overlap >= $minoverlap} {
				lappend matchlist [list $mpos $pos $overlap]
			}
		}
	}
	# sort matchlist, and use best matches
	set matchlist [lsort -index 2 -integer -decreasing $matchlist]
	unset -nocomplain a
	set exnum [expr {$totnum-$usednum}]
	list_foreach {mpos pos} $matchlist {
		if {[info exists a(m,$mpos)] || [info exists a(p,$pos)]} continue
		set line [lindex $maxlist $mpos]
		lappend line $numthreads $exnum
		set threads [sv_threadadd $threads $pos $line]
		set a(m,$mpos) 1
		set a(p,$pos) 1
	}
	set mpos -1
	foreach line $maxlist {
		incr mpos
		if {[info exists a(m,$mpos)]} continue
		lappend line $numthreads $exnum
		lappend threads [list $line]
	}
	return $threads
}

proc sv_threaddist {thread} {
	set temp {}
	list_foreach {p s e} $thread {
		lappend temp [expr {($e+$s)/2}]
	}
	set m [lmath_average $temp]
}

proc sv_checkthreads {threads start mode resultVar} {
	upvar $resultVar result
	global infoa
	set clustmin 6
	set trfpos $infoa(trfpos)
	set distpos $infoa(distpos)
	set weight1pos $infoa(weight1pos)
	set weight2pos $infoa(weight2pos)
	# prune
	set pos 0
	foreach thread $threads {
		if {[llength $thread] >= 16} {
			set num [lindex $thread 0]
			set thread [lrange $thread end-15 end]
			lset threads $pos $thread
		}
		incr pos
	}
	set mainthread [list_shift threads]
	# list_subindex $mainthread {0 1 2}
	# check for new main
	set mpos [lindex $mainthread end 0]
	if {[expr {$start - $mpos}] > 1000} {
		set bestpos -1
		set bestd $mode
		set pos -1
		foreach thread $threads {
			incr pos
			set tpos [lindex $thread end 0]
			if {[expr {$start - $tpos}] < 200} {
				set cdist [lmath_average [lrange [lindex $thread end] 1 2]]
				set d [expr {abs($cdist-$mode)}]
				if {$d < $bestd} {set bestd $d; set bestpos $pos}
			}
		}
		if {$bestd < 50} {
			set mainthread [list_pop threads $bestpos]
		}
	}
	set lmode [expr {round([sv_threaddist $mainthread])}]
	set pos -1
	set resultthreads [list $mainthread]
	foreach thread $threads {
		incr pos
		# list_subindex $thread {0 1 2}
		set x [lindex $thread end 0]
		set diff [expr {$start - $x}]
		if {$diff > 50} {
			# thread ends
			if {[llength $thread] < 2} continue
			set parts [list_remdup [list_concat [list_subindex $thread 3]]]
			if {[llength $parts] < 7} continue
			set parts [lsort -index 2 -integer $parts]
			set chr [lindex $parts 0 0]
			set cstart [lmath_min [list_subindex $parts 1]]
			set cend [lmath_max [list_subindex $parts 2]]
			set chr2 [lindex $parts 0 5]
			set start2 [lmath_min [list_subindex $parts 6]]
			set end2 [lmath_max [list_subindex $parts 7]]
			set patchsize [expr {$cend-$cstart+1}]
			set gapsize [expr {round( [lmath_average [list_subindex $parts $distpos]] )}]
			set totnum [llength $parts]
			set exnum [lmath_average [list_subindex $thread 5]]
			set tnum [lmath_average [list_subindex $thread 4]]
			set problems {}
			# test heterozygosity (and check if len exists)
			# ---------------------------------------------
			set tposs [list_subindex $thread 0]
			set mposs [list_subindex $mainthread 0]
			set common [list_common $mposs $tposs]
			if {[llength $common] <= 1} {
				set zyg hom
			} else {
				set commonp [expr {double([llength $common])/[llength $tposs]}]
				if {$commonp < 0.35} {set zyg hom} else {set zyg het}
			}
			# linreg
			# ------
			foreach {b1 sd1 b2 sd2} {0 0 0 0} break
			if {[llength $thread] > 2} {
				set xs [list_subindex $thread 0]
				set ys1 [list_subindex $thread 1]
				set ys2 [list_subindex $thread 2]
				set ys [lmath_calc [lmath_calc $ys1 + $ys2] / 2]
				set mid [expr {[llength $xs]/2}]
				set xs1 [lrange $xs 0 $mid]
				set ys1 [lrange $ys 0 $mid]
				if {![expr {[llength $xs]%2}]} {incr mid}
				foreach {a b1 sd1} [linreg $xs1 $ys1] break
				set b1 [oformat $b1]; set sd1 [oformat $sd1]
				set xs2 [lrange $xs $mid end]
				set ys2 [lrange $ys $mid end]
				foreach {a b2 sd2} [linreg $xs2 $ys2] break
				set b2 [oformat $b2]; set sd2 [oformat $sd2]
			}
			# quality
			# -------
			set diff [expr {$gapsize - $lmode}]
			if {$diff > 0} {
				set type del
			} else {
				set type ins
				set diff [expr {-$diff}]
			}
			set num [llength $parts]
			set weight1 [lmath_average [list_subindex $parts $weight1pos]]
			set weight2 [lmath_average [list_subindex $parts $weight2pos]]
			set weight [expr {round ([min $weight1 $weight2])}]
			if {$type eq "del"} {
				set opsdiff [expr {$patchsize-$mode}]
			} else {
				set opsdiff [expr {$patchsize+$diff-$mode}]
			}
			set psdiff [expr {abs($opsdiff)}]
			# test trf
			# --------
			set numnontrf [llength [list_find -exact [list_subindex $parts $trfpos] 0]]
			if {$numnontrf < $clustmin} {
				lappend problems trfartefact
			} elseif {$numnontrf < [expr {0.2*$num}]} {
				lappend problems trfbad
			} elseif {$numnontrf < [expr {0.4*$num}]} {
				lappend problems trfpoor
			}
			# smallins
			# -------
			if {$diff <= $infoa(maxsize)} {
				if {$zyg eq "hom"} {
					lappend problems msmall
				} else {
					lappend problems hsmall
				}
			} elseif {($type eq "ins") && ($diff <= 115)} {
				if {($b1 < -0.1) && ($b2 > 0.1)} {
					lappend problems pdip
				}
			}
			# score
			set score [svscore $mode $type $diff $zyg $problems $gapsize $num $numnontrf $weight $patchsize $b1 $sd1 $b2 $sd2 $totnum $opsdiff $tnum $exnum]
			# add to result
			# -------------
			if {$cend > $start2} {
				set b $start2
				set e $cend
			} else {
				set b $cend
				set e $start2
			}
			lappend result [list $chr $b $e $type [expr {round($diff)}] $zyg $score $problems $cstart $cend $chr2 $start2 $end2 $gapsize $num $numnontrf $weight $patchsize $b1 $sd1 $b2 $sd2 $totnum [oformat $opsdiff 0] $tnum $exnum]
		} else {
			lappend resultthreads $thread
		}
	}
	return $resultthreads
}

proc kde_distcluster {dtable} {
	global infoa

	set distpos $infoa(distpos)
	# r [list_subindex $dtable $distpos]
	set prev [lindex $dtable 0 $distpos]
	set num 0
	set list {}
	set part {}
	set pairs {}
	set end 0
	foreach line $dtable {
		set v [lindex $line $distpos]
		set d [expr {$v-$prev}]
		if {$d == 0} {
			incr num
			lappend pairs $line
		} elseif {$d < 50} {
			lappend part [list $prev $num $pairs]
			set pairs [list $line]
			set num 1
			set prev $v
		} else {
			lappend part [list $prev $num $pairs]
			set pairs [list $line]
			set tot [lmath_sum [list_subindex $part 1]]
			set s [lindex $part 0 0]
			set e [lindex $part end 0]
			if {$tot > 7} {
				lappend list [list $s $e $tot $part]
			}
			set part {}
			set num 1
			set prev $v
		}
	}
	lappend part [list $prev $num]
	set tot [lmath_sum [list_subindex $part 1]]
	set s [lindex $part 0 0]
	set e [lindex $part end 0]
	if {$tot > 7} {
		lappend list [list $s $e $tot $part]
	}
	# find clusters
	set clusters {}
	# gausian kernel, bandwidth = 7
	set kernel {0.0100 0.0146 0.0209 0.0293 0.0402 0.0540 0.0711 0.0918 0.1161 0.1438 0.1746 0.2076 0.2420 0.2763 0.3091 0.3388 0.3639 0.3830 0.3949 0.3989 0.3949 0.3830 0.3639 0.3388 0.3091 0.2763 0.2420 0.2076 0.1746 0.1438 0.1161 0.0918 0.0711 0.0540 0.0402 0.0293 0.0209 0.0146 0.0100}
	# gausian kernel, bandwidth = 8
	# set kernel {0.0238 0.0317 0.0417 0.0540 0.0688 0.0863 0.1065 0.1295 0.1550 0.1826 0.2119 0.2420 0.2721 0.3011 0.3282 0.3521 0.3719 0.3867 0.3958 0.3989 0.3958 0.3867 0.3719 0.3521 0.3282 0.3011 0.2721 0.2420 0.2119 0.1826 0.1550 0.1295 0.1065 0.0863 0.0688 0.0540 0.0417 0.0317 0.0238}
	# gausian kernel, bandwidth = 9
	# set kernel {0.0262 0.0338 0.0430 0.0540 0.0670 0.0822 0.0995 0.1190 0.1406 0.1640 0.1890 0.2152 0.2420 0.2687 0.2948 0.3194 0.3419 0.3614 0.3774 0.3892 0.3965 0.3989 0.3965 0.3892 0.3774 0.3614 0.3419 0.3194 0.2948 0.2687 0.2420 0.2152 0.1890 0.1640 0.1406 0.1190 0.0995 0.0822 0.0670 0.0540 0.0430 0.0338 0.0262}
	# gausian kernel, bandwidth = 10
	# set kernel {0.0060 0.0079 0.0104 0.0136 0.0175 0.0224 0.0283 0.0355 0.0440 0.0540 0.0656 0.0790 0.0940 0.1109 0.1295 0.1497 0.1714 0.1942 0.2179 0.2420 0.2661 0.2897 0.3123 0.3332 0.3521 0.3683 0.3814 0.3910 0.3970 0.3989 0.3970 0.3910 0.3814 0.3683 0.3521 0.3332 0.3123 0.2897 0.2661 0.2420 0.2179 0.1942 0.1714 0.1497 0.1295 0.1109 0.0940 0.0790 0.0656 0.0540 0.0440 0.0355 0.0283 0.0224 0.0175 0.0136 0.0104 0.0079 0.0060}
	set ke [llength $kernel]
	# join [list_subindex $list {0 1 2}] \n
	list_foreach {s e num part} $list {
		if {$num < 8} continue
		set plen [expr {$e-$s}]
		if {(($s > 1000) && ($plen < 150))
			|| (($s > 600) && ($plen < 100))
			|| (($num < 20) && ($plen < 100))
			|| ($plen < 50)
		} {
			set elements [list_subindex $part 2]
			if {[llength $elements] <= 1} {
				set elements [lindex $elements 0]
			} else {
				set elements [list_concat {*}$elements]
			}
			if {[llength $elements] < 8} continue
			lappend clusters [list [lmath_average [list_subindex $elements $distpos]] [llength $elements] $elements]
			continue
		}
		set len [expr {$plen+2*$ke}]
		set data [list_fill $len 0]
		set elementdata [list_fill $len {}]
		unset -nocomplain a
		list_foreach {p v l} $part {
			set p [expr {$p-$s+$ke}]
			lset data $p $v
			lset elementdata $p $l
		}
		# convolving with the filter over this is the same as doing the discrete kde
		set data [lmath_filter $data $kernel]
		set maxima [sv_maxima $data 3]
		set end 0
		list_foreach {maxpos1 ph maxpos2} $maxima {
			set peak $ph
			set pos $maxpos1
			set minpos [expr {$pos-1}]
			set minh $peak
			set stop [max [expr {$pos-50}] 0 $end]
			while {$pos >= $stop} {
				incr pos -1
				if {$pos < 0} break
				set h [lindex $data $pos]
				if {$h < $minh} {
					set minpos $pos
					set minh $h
				} elseif {$h > [expr {$minh+0.1}]} {
					break
				}
			}
			set start $minpos
			set pos $maxpos2
			set minpos [expr {$pos-1}]
			set minh $peak
			set len [llength $data]
			set stop [min [expr {$pos+50}] $len]
			while {$pos < $stop} {
				incr pos
				set h [lindex $data $pos]
				if {$h < $minh} {
					set minpos $pos
					set minh $h
				} elseif {$h > [expr {$minh+0.1}]} {
					break
				}
			}
			set end $minpos
			set elements [lrange $elementdata $start $end]
			if {[llength $elements] <= 1} {
				set elements [lindex $elements 0]
			} else {
				set elements [list_concat {*}$elements]
			}
			if {[llength $elements] < 8} continue
			lappend clusters [list [lmath_average [list_subindex $elements $distpos]] [llength $elements] $elements]
		}
	}

	# join [list_subindex $clusters {0 1}] \n
	return $clusters
}

proc gauss u {expr {(1/sqrt(2*acos(-1)))*exp(-0.5*$u*$u)}}
proc makekernel {bw} {
	set size [expr {$bw*2}]
	set hkernel {}
	for {set i 0} {$i < $size} {incr i} {
		set u [expr {$i/double($bw)}]
		lappend hkernel [format %.4f [gauss $u]]
	}
	set kernel [list_concat [list_reverse [lrange $hkernel 1 end]] $hkernel]
}

if 0 {
	# make kernels
	# epanechnikov kernel
	set scale 50
	set h 9.5
	set hkernel {} 
	set max [expr {(3/4.0)*(1)}]
	foreach u {0 1 2 3 4 5 6 7 8 9} {
		set u [expr {$u/double($h)}]
		lappend hkernel [expr {round($scale*(3/4.0)*(1-$u*$u)/$max)}]
	}
	puts [list set kernel [list_concat [list_reverse [lrange $hkernel 1 end]] $hkernel]]
	# gausian kernel
	proc gauss u {expr {(1/sqrt(2*acos(-1)))*exp(-0.5*$u*$u)}}
	set size 22
	set bw 9
	set hkernel {}
	for {set i 0} {$i < $size} {incr i} {
		set u [expr {$i/double($bw)}]
		lappend hkernel [format %.4f [gauss $u]]
	}
	puts [list set kernel [list_concat [list_reverse [lrange $hkernel 1 end]] $hkernel]]

}

proc r {data} {
	if {[isint [lindex $data 0]]} {
		return "c([join $data ,])"
	}
	set temp {}
	foreach el $data {
		lappend temp [format %.3f $el]
	}
	return "c([join $temp ,])"
}

proc linreg {xs ys} {
	if {[llength $xs] < 2} {
		return [list Nan Nan Nan]
	}
	set xm [lmath_average $xs]
	set ym [lmath_average $ys]
	set dx [lmath_calc $xs - $xm]
	set dy [lmath_calc $ys - $ym]
	set xsum [lmath_sum [lmath_calc $dx * $dy]]
	set ysum [lmath_sum [lmath_calc $dx * $dx]]
	if {$ysum != 0} {
		set b [expr {$xsum/$ysum}]
	} else {
		return [list Nan Nan Nan]
	}
	set a [expr {$ym - $b * $xm}]
	set eys [lmath_calc $a + [lmath_calc $xs * $b]]
	set res [lmath_calc $eys - $ys]
	set s 0.0
	foreach r $res {
		set s [expr {$s + $r*$r}]
	}
	if {$s != 0} {
		set sd [expr {sqrt($s/([llength $xs]-2))}]
	} else {
		set sd 0.0
	}
	return [list $a $b $sd]
}

proc svscore {mode type size zyg problems gapsize num numnontrf weight patchsize b1 sd1 b2 sd2 totnum opsdiff tnum exnum} {
	set psdiff [expr {abs($opsdiff)}]
	set maxpatch [expr {1.8*$mode}]
	set hmode [expr {$mode/2}]
	if {$patchsize > $maxpatch} {
		set score 0
	} elseif {$type eq "inv"} {
		set score 0
		if {$patchsize > 1000} {
			set score 0
		} elseif {$weight == 0} {
			if {$num < 10} {set score 0} else {set score 2}
		} else {
			set score 0
			if {$patchsize < $mode} {
				if {$psdiff < 100} {
					incr score 3
				} elseif {$psdiff < 200} {
					incr score 2
				} elseif {$psdiff < 300} {
					incr score 1
				}
			} else {
				if {$psdiff < 40} {
					incr score 3
				} elseif {$psdiff < 100} {
					incr score 2
				}
			}
			if {$num >= 50} {
				incr score 3
			} elseif {$num > 20} {
				incr score 1
			}
			if {$weight >= 30} {
				incr score 4
			} elseif {$weight > 20} {
				incr score 3
			} elseif {$weight > 5} {
				incr score 1
			}
		}
		return $score
	} elseif {$type eq "trans"} {
		set score 7
		if {$weight > 20} {
			set score 5
		} elseif {$weight > 5} {
			set score 3
		} elseif {$weight > 0} {
			set score 2
		} else {
			set score 1
		}
		if {$num > 100} {
			incr score 4
		} elseif {$num > 50} {
			incr score 2
		}
		if {[inlist $problems many]} {set score 0}
		return $score
	} elseif {$type eq "del"} {
		set score 0
		if {$size <= 10} {
			set score 0
		} elseif {$size > 400} {
			if {($num >= 20) && ($weight >= 10)} {
				set score 10
			} elseif {($num >= 10) && ($weight > 0)} {
				set score 9
			} else {
				set score 7
			}
		} else {
			if {$num >= 50} {
				incr score 4
			} elseif {$num > 20} {
				incr score 1
			}
			if {$weight >= 30} {
				incr score 3
			} elseif {$weight > 20} {
				incr score 2
			} elseif {$weight > 5} {
				incr score 1
			} elseif {$weight == 0} {
				incr score -1
			}
			if {$patchsize < $mode} {
				if {$psdiff < 100} {
					incr score 3
				} elseif {$psdiff < 150} {
					incr score 2
				} elseif {$psdiff < 250} {
					incr score 1
				}
			} else {
				if {$psdiff < 40} {
					incr score 3
				} elseif {$psdiff < 100} {
					incr score 2
				}
			}
		}
		if {$size <= 100} {
			if {$patchsize < $hmode} {set score [min $score 2]}
		}
	} elseif {$type eq "ins"} {
		if {$size <= 10} {
			set score 0
		} elseif {($size <= 115) && ($b1 < -0.2) && ($b2 > 0.2)} {
			set score 1
		} else {
			set score 0
			if {$num >= 50} {
				incr score 4
			} elseif {$num > 20} {
				incr score 1
			}
			if {$weight >= 30} {
				incr score 3
			} elseif {$weight > 20} {
				incr score 2
			} elseif {$weight > 5} {
				incr score 1
			} elseif {$weight == 0} {
				incr score -1
			}
			if {$patchsize < $mode} {
				if {$psdiff < 100} {
					incr score 3
				} elseif {$psdiff < 150} {
					incr score 2
				} elseif {$psdiff < 250} {
					incr score 1
				}
			} else {
				if {$psdiff < 40} {
					incr score 3
				} elseif {$psdiff < 100} {
					incr score 2
				}
			}
		}
		if {$size <= 100} {
			if {$patchsize < $hmode} {set score [min $score 2]}
		}
	} else {
		set score 0
	}
	if {$score > 4} {
		if {$patchsize > 500} {
			set score 4
		}
		if {$exnum > 10} {set score 4}
	}
	if {$score > 3} {
		if {$size <= 30} {
			set score 3
		}
		if {$tnum > 3} {set score 3}
	}
	if {$score > 2} {
		if {$tnum > 4} {set score 2}
		if {$exnum > 20} {set score 2}
	}
	if {[inlist $problems trfartefact]} {
		set score 0
	} elseif {[inlist $problems trfbad] && ($score > 1)}  {
		set score 1
	}
#	if {$score > 5} {
#		if {($sd1 ne "Nan") && (abs($sd1) > 0.15)} {
#			set score 5
#		} elseif {($sd2 ne "Nan") && (abs($sd2) > 0.15)} {
#			set score 5
#		}
#	}
	return $score
}

proc svwindow_checkinv {table rtable} {
	global infoa
	set result {}
	set distpos $infoa(distpos)
	set start1pos $infoa(start1pos)
	set end1pos $infoa(end1pos)
	set start2pos $infoa(start2pos)
	set end2pos $infoa(end2pos)
	set mode $infoa(mode)
	set weight1pos $infoa(weight1pos)
	set weight2pos $infoa(weight2pos)
	set chr1pos $infoa(chr1pos)
	set trfpos $infoa(trfpos)
	set clustmin 6
	set step 5
	set chr [lindex [lindex $rtable 0] $chr1pos]
	set rtable [lsort -integer -index $end1pos $rtable]
	set lists [cluster_fts $rtable $end1pos 30]
	foreach list $lists {
		set list [lsort -integer -index $end1pos $list]
		set groups {}
		foreach line $list {
			set found 0
			set tpos [lindex $line $end1pos]
			set tdist [lindex $line $distpos]
			set pos -1
			foreach group $groups {
				incr pos
				set gpos [lindex $group end $end1pos]
				set gdist [lindex $group end $distpos]
				set maxdist [expr {min(400,max(100,round(0.3*$gdist)))}]
				set diffexpected [expr {abs($tdist - ($gdist-($gpos-$tpos)))}]
				if {$diffexpected < $maxdist} {
					lappend group $line
					lset groups $pos $group
					set found 1
					break
				}
			}
			if {!$found} {
				lappend groups [list $line]
			}
		}
		# join [lindex $groups 0] \n
		foreach group $groups {
			if {[llength $group] < $clustmin} continue
			set size [expr {[lindex $group 0 $distpos]-$mode}]
			set cstart [lmath_min [list_subindex $group $start1pos]]
			set cend [lmath_max [list_subindex $group $end1pos]]
			set cstart [lindex $group 0 $end1pos]
			set cend [lindex $group end $end1pos]
			set patchsize [expr {$cend-$cstart+1}]
			if {$patchsize <= 1} continue
			set start2 [lmath_min [list_subindex $group $start2pos]]
			set end2 [lmath_max [list_subindex $group $end2pos]]
			set num [llength $group]
			# test heterozygosity
			set list {}
			foreach line $table {
				set pos [lindex $line $end1pos]
				if {($pos >= $cstart) && ($pos < $cend)} {
					lappend list $line
				}
			}
			set totnum [llength $list]
			if {[llength $list] > [expr {$num*0.4}]} {set zyg het} else {set zyg hom}
			# add to result
			set weight1 [lmath_average [list_subindex $group $weight1pos]]
			set weight2 [lmath_average [list_subindex $group $weight2pos]]
			set weight [expr {round ([min $weight1 $weight2])}]
			set opsdiff [expr {$patchsize-$mode}]
			set psdiff [expr {abs($opsdiff)}]
			# linreg
			# ------
			foreach {a b sd} [linreg [list_subindex $group $end1pos] [list_subindex $group $distpos]] break
			set b [oformat $b]; set sd [oformat $sd]
			set score [svscore $mode inv $size $zyg {} {} $num $num $weight $patchsize $b $sd {} {} $totnum $opsdiff {} {}]
			lappend result [list $chr $cend $start2 inv $size $zyg $score {} $cstart $cend $chr $start2 $end2 {} $num {} $weight $patchsize $b $sd {} {} $totnum [oformat $opsdiff 0] {} {}]
		}
		
	}
	return $result
}

proc svwindow_checktrans {table ctable mode} {
	global infoa
	set distpos $infoa(distpos)
	set start1pos $infoa(start1pos)
	set end1pos $infoa(end1pos)
	set start2pos $infoa(start2pos)
	set end2pos $infoa(end2pos)
	set mode $infoa(mode)
	set weight1pos $infoa(weight1pos)
	set weight2pos $infoa(weight2pos)
	set chr1pos $infoa(chr1pos)
	set clustmin 10
	set step 5
	set chr [lindex [lindex $ctable 0] $chr1pos]
	set chr2pos $infoa(chr2pos)
	set start2pos $infoa(start2pos)
	unset -nocomplain a
	foreach line $ctable {
		set chr2 [lindex $line $chr2pos]
		set dist [expr {[lindex $line $start2pos] - [lindex $line $end1pos]}]
		lset line $distpos $dist
		lappend a($chr2) $line
	}
	set result {}
	foreach chr2 [array names a] {
		set rtable $a($chr2)
		set rtable [lsort -integer -index $distpos $rtable]
		set lists [cluster_fts $rtable $distpos 40]
		foreach rtable $lists {
			if {[llength $rtable] < 5} continue
			set rtable [lsort -integer -index $end1pos $rtable]
			# join [lsort -integer -index 2 $rtable] \n
			set srtables [cluster_fts $rtable]
			foreach rtable $srtables {
				if {[llength $rtable] < $clustmin} continue
				set size [expr {[lindex $rtable 0 $distpos]-$mode}]
				set cstart [lmath_min [list_subindex $rtable $start1pos]]
				set cend [lmath_max [list_subindex $rtable $end1pos]]
				set patchsize [expr {$cend-$cstart+1}]
				if {$patchsize <= 1} continue
				set start2 [lmath_min [list_subindex $rtable $start2pos]]
				set end2 [lmath_max [list_subindex $rtable $end2pos]]
				set opsdiff [expr {$patchsize-$mode}]
				set psdiff [expr {abs($opsdiff)}]
				set num [llength $rtable]
				# test heterozygosity
				set list {}
				foreach line $table {
					set pos [lindex $line $end1pos]
					if {($pos >= $cstart) && ($pos < $cend)} {
						lappend list $line
					}
				}
				set totnum [llength $list]
				if {[llength $list] > [expr {$num*0.4}]} {set zyg het} else {set zyg hom}
				# add to result
				set weight1 [lmath_average [list_subindex $rtable $weight1pos]]
				set weight2 [lmath_average [list_subindex $rtable $weight2pos]]
				set weight [expr {round ([min $weight1 $weight2])}]
				set pos2 [lindex $rtable 0 $start2pos]
				# linreg
				# ------
				set xs [list_subindex $rtable $end1pos]
				set ys [lmath_calc [list_subindex $rtable $start2pos] - [list_subindex $rtable $end1pos]]
				foreach {temp b sd} [linreg $xs $ys] break
				set b [oformat $b]; set sd [oformat $sd]
				set score [svscore $mode trans $pos2 $zyg {} {} $num $num $weight $patchsize $b $sd {} {} $totnum $opsdiff {} {}]
				lappend result [list $chr $cend [expr {$cend+$mode}] trans $pos2 $zyg $score {} $cstart $cend $chr2 $start2 $end2 {} $num {} $weight $patchsize $b $sd {} {} $totnum [oformat $opsdiff] {} {}]
			}
		}
	}
	if {[llength [list_remdup [list_subindex $result 10]]] > 1} {
		set pos 0
		foreach line $result {
			lset result $pos 7 many
			lset result $pos 6 0
			incr pos
		}
	}
	return $result
}

proc svloadtrf {trf chr pstart start end trfposs trflistVar trflineVar} {
	upvar $trflistVar trflist
	upvar $trflineVar trfline
	set temp {}
	list_foreach {trfchr trfstart trfend} $trflist {
		set chrcompar [chr_compare $trfchr $chr]
		if {($trfchr == $chr) && ($trfend >= $start) && ($trfstart < $end)} {
			lappend temp [list $trfchr $trfstart $trfend]
		}
	}
	set trflist $temp
	while {![eof $trf]} {
		foreach {trfchr trfstart trfend} $trfline break
		set chrcompar [chr_compare $trfchr $chr]
		if {$chrcompar > 1} break
		if {$trfstart > $end} break
		set trflen [expr {$trfend-$trfstart}]
		if {($trfchr == $chr) && ($trfend >= $pstart) && ($trflen > 30)} {
			incr trfstart -6
			incr trfend +6
			lappend trflist [list $trfchr $trfstart $trfend]
		}
		set trfline [get_region $trf $trfposs]
	}
}

if 0 {

	set header {match check sample chr patchstart pos type size zyg problems gapsize/chr2 quality numreads numnontrf weight patchsize slope1 sd1 slope2 sd2 totnum psdiff threads exnum chr-2 patchstart-2 pos-2 type-2 size-2 zyg-2 problems-2 gapsize/chr2-2 quality-2 numreads-2 numnontrf-2 weight-2 patchsize-2 slope1-2 sd1-2 slope2-2 sd2-2 totnum-2 psdiff-2 threads-2 exnum-2}
	set cor [list_cor $header {quality type size zyg problems gapsize/chr2 numreads numnontrf weight patchsize slope1 sd1 slope2 sd2 totnum psdiff threads exnum}]
	set mode $infoa(mode)

	set line [split $line \t]

	foreach {score type size zyg problems gapsize num numnontrf weight patchsize b1 sd1 b2 sd2 totnum psdiff tnum exnum} [list_sub $line $cor] break
	set nscore [svscore $mode $type $size $zyg $problems $gapsize $num $numnontrf $weight $patchsize $b1 $sd1 $b2 $sd2 $totnum $psdiff $tnum $exnum]

	set file /complgen/sv/workGS103-9-paired-sv.tsv

}

proc svrescore {file} {

	set ofile [file dir $file]/rs_[file tail $file]
	file delete -force $ofile
	catch {close $f} ; catch {close $o}
	set f [open $file]
	set o [open $ofile w]
	set header [gets $f]
	puts $o $header
	set cor [list_cor $header {quality type size zyg problems gapsize/chr2 numreads numnontrf weight patchsize slope1 sd1 slope2 sd2 totnum psdiff threads exnum}]
	set scorepos [lsearch $header quality]
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		foreach {score type diff zyg problems gapsize num numnontrf weight patchsize b1 sd1 b2 sd2 totnum opsdiff tnum exnum} [list_sub $line $cor] break
		if {$type eq "del"} {
			set mode [expr {$gapsize-$diff}]
			break
		}
	}
	seek $f 0
	gets $f
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		foreach {score type diff zyg problems gapsize num numnontrf weight patchsize b1 sd1 b2 sd2 totnum opsdiff tnum exnum} [list_sub $line $cor] break
		set nscore [svscore $mode $type $diff $zyg $problems $gapsize $num $numnontrf $weight $patchsize $b1 $sd1 $b2 $sd2 $totnum $opsdiff $tnum $exnum]
		lset line $scorepos $nscore
		puts $o [join $line \t]
	}
	close $f
	close $o
	putslog "made svrescore $ofile"
#sv_evaluate $ofile

}

proc sv_evaluate {file {field check}} {

#	set file /complgen/sv/workGS103-9-paired-sv.tsv
#	set file /complgen/sv/rworkGS103-9-paired-sv.tsv
#	set file /complgen/sv/refGS103-20-paired-sv.tsv
#	set file /complgen/sv/rs_refGS103-20-paired-sv.tsv
	set f [open $file]
	set header [gets $f]
	set checkpos [lsearch $header $field]
	set scorepos [lsearch $header quality]
	set problemspos [lsearch $header problems]
	unset -nocomplain a
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set problems [lindex $line $problemspos]
		if {[inlist $problems msmall]} continue
		set check [lindex $line $checkpos]
		if {$check eq ""} continue
		set score [lindex $line $scorepos]
		set error [catch {dict get $a($check) $score} num]
		if {$error} {
			if {![info exists a($check)]} {
				set a($check) [dict create $score 1]
			} else {
				dict set a($check) $score 1
			}
		} else {
			incr num
			dict set a($check) $score $num
		}
	}
	close $f
	set header {check 0 1 2 3 4 5 6 7 8 9 10}
	set table [join $header \t]
	foreach check [array names a] {
		set line [list $check]
		foreach num [lrange $header 1 end] {
			if {[catch {dict get $a($check) $num} value]} {set value 0}
			lappend line $value
		}
		append table \n[join $line \t]
	}
	set w [edit]
	$w.editor set $table

}

if 0 {

package require Tclx
signal -restart error SIGINT
lappend auto_path ~/dev/genomecomb/lib ~/dev/genomecomb/lib-exp
package require Extral
#cd /complgen/sv

set trffile /complgen/refseq/hg18/reg_hg18_simpleRepeat.tsv
set pairfile /complgen/projects/cmt71/cmt71_02_a/sv/cmt71_02_a-22-paired.tsv.rz
set pairfile /media/solid/kr1270/d816_3/sv/d816_3-22-paired.tsv.gz

foreach {score type diff zyg problems gapsize num numnontrf weight patchsize b1 sd1 b2 sd2 totnum opsdiff} [list_sub $line $cor] break
set nscore [svscore $mode $type $diff $zyg $problems $gapsize $num $numnontrf $weight $patchsize $b1 $sd1 $b2 $sd2 $totnum $opsdiff]

set pairfile GS103/GS103-20-paired.tsv.rz
set pairfile sv78-20-pairs.tsv
set pairfile sv79-20-pairs.tsv
set pairfile GS102/GS102-9-paired.tsv.rz
set pairfile GS103/GS103-9-paired.tsv.rz
set pairfile GS102/GS102-20-paired.tsv.rz
set pairfile GS103/GS103-20-paired.tsv.rz

}

proc svfind {pairfile trffile} {
#putslog [list svfind $pairfile $trffile]
	set bpairfile [gzroot $pairfile]
	set outfile [file root $bpairfile]-sv.tsv
	set windowsize 50
	set maxpairs 10000
	# set bandwidth 50

	catch {close $trf}
	catch {close $o}
	catch {close $f}
#	if {[file exists $outfile]} {
#		putslog "$outfile exists: skipping"
#		return
#	}
	set lognum [expr {1000000 - 100000%$windowsize}]
	global infoa
	set f [open $bpairfile.numinfo]
	tsv_open $f
	array set infoa [split [string trim [read $f]] \n\t]
	close $f
	set bandwidth [expr {($infoa(max10)-$infoa(min10)/2)}]
	set infoa(step) 5
	set min $infoa(min)
	set max $infoa(max)
	set mode $infoa(mode)
	set infoa(hmode) [expr {$infoa(mode)-$infoa(mode)%$infoa(step)}]
	set infoa(hremove) [list_fill 5 [expr {$infoa(hmode)-2*$infoa(step)}] $infoa(step)]
	set infoa(kernel) [makekernel $bandwidth]
	set infoa(maxsize) $bandwidth
	set f [gzopen $pairfile]
	set header [gets $f]
	set iheader {chromosome start1 end1 weight1 numl chr2 start2 end2 weight2 numr type dist}
	set poss [list_cor $header $iheader]
	if {[lindex $poss 0] == -1} {lset poss 0 [lsearch $header chr1]}
	array set infoa [list \
		start1pos [lsearch $iheader start1] \
		end1pos [lsearch $iheader end1] \
		start2pos [lsearch $iheader start2] \
		end2pos [lsearch $iheader end2] \
		typepos [lsearch $iheader type] \
		chr1pos [lsearch $iheader chromosome] \
		distpos [lsearch $iheader dist] \
		weight1pos [lsearch $iheader weight1] \
		weight2pos [lsearch $iheader weight2] \
		chr2pos [lsearch $iheader chr2] \
		start2pos [lsearch $iheader start2] \
		trfpos [llength $iheader] \
	]
	set infoa(pairposs) [list $infoa(typepos) [lsearch $iheader start1] $infoa(end1pos) [lsearch $iheader start2] [lsearch $iheader end2]]
	set pairposs $infoa(pairposs)
	set end1pos $infoa(end1pos)
	set distpos $infoa(distpos)
	set typepos $infoa(typepos)
#check
#lassign {9790000 9790550} dbgstart dbgstop
#catch {close $f}
#set f [svtools_aprgoto $pairfile $dbgstart]
#set outfile test-sv.tsv
#check
	set dir [file dir [file_absolute $outfile]]
	set o [open $outfile.temp w]
	puts $o [join {check chromosome begin end type size zyg quality problems start1 end1 chr2 start2 end2 gapsize numreads numnontrf weight patchsize slope1 sd1 slope2 sd2 totnum psdiff threads exnum} \t]
	set list {}
	set mainrtable {}
	set rtable {}
	set rlist {}
	set prlist {}
	set mainctable {}
	set ctable {}
	set clist {}
	set pclist {}
	unset -nocomplain plist
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set line [list_sub $line $poss]
		set start [lindex $line $end1pos]
		if {[isint $start]} break
	}
	set start [expr {$start-$start%$windowsize}]
	set pstart [expr {$start-$windowsize}]
	set end [expr {$start + $windowsize -1}]
	set chr [lindex $line $infoa(chr1pos)]
	# prepare trf data
	set trffile [gztemp $trffile]
	set trf [gzopen $trffile]
	set trfposs [open_region $trf]
	set trflist {}
	chrindexseek $trffile $trf $chr
	while {![eof $trf]} {
		set trfline [get_region $trf $trfposs]
		foreach {trfchr trfstart trfend} $trfline break
		set chrcompar [chr_compare $trfchr $chr]
		if {$chrcompar >= 0} break
	}
	# go over file and find sv
	set prevmains [expr {$mode-50}]
	set prevmaine [expr {$mode+50}]
	set temp {}
	for {set i -400} {$i < 0} {incr i 50} {lappend temp [list $i $prevmains $prevmaine {}]}
	set threads [list $temp]
	while {![eof $f]} {
		if {![expr $start % $lognum]} {putslog $start}
		# load trf regions
		svloadtrf $trf $chr $pstart $start $end $trfposs trflist trfline
		# load pairs
		while {![eof $f]} {
			if {[llength $line]} {
				foreach {type start1 refpos start2 end2} [list_sub $line $pairposs] break
				if {$refpos >= $end} break
				set intrf 0
				list_foreach {trfchr trfstart trfend} $trflist {
					if {($start1 >= $trfstart) && ($refpos <= $trfend)} {
						set intrf 1
					}
					if {($start2 >= $trfstart) && ($end2 <= $trfend)} {
						set intrf 1
					}
				}
				lappend line $intrf
				if {$type eq "c"} {
					lappend clist $line
				} elseif {$type eq "r"} {
					lappend rlist $line
				} elseif {$type eq "s"} {
				} elseif {$type eq "i"} {
					lappend list $line
				} else {
					lappend list $line
				}
			}
			set line [getnotempty $f]
			set line [list_sub $line $poss]
		}
		if {[info exists plist] && ([llength $plist] < $maxpairs)} {
			# check for sv
			set result {}
			set minadd [expr {$start - $windowsize/2}]
			set maxadd [expr {$start + $windowsize/2}]
			# check for insertions and deletions
			set table [list_concat $plist $list]
			set dists [list_subindex $table $distpos]
			set maxima [kde_maxima $dists $infoa(maxsize)]
#if {$start >= $dbgstop} {error STOPPED}
# draw $::kdedata $::kdes
			set threads [sv_addtothreads $threads $start $maxima $table]
#			set onum [llength $maxima]
#			if {$onum < 4} {
#				set threads [sv_addtothreads $threads $start $maxima $table]
#			}
#check
			set temp {}
			set threads [sv_checkthreads $threads $start $mode temp]
#putsvars start
#puts [t]
#if {[llength $temp]} {error STOP}
			if {[llength $temp]} {lappend result {*}$temp}
			# check for inversions
			if {[llength $prlist] > 0} {
				lappend rtable {*}$prlist
				lappend mainrtable {*}$plist
				set rlastpos $start
			} elseif {[llength $rtable] && ([expr {$start-$rlastpos}] > 50)} {
				set temp [svwindow_checkinv $mainrtable $rtable]
				if {[llength $temp]} {
					lappend result {*}$temp
				}
				set mainrtable {}
				set rtable {}
			}
			# check for translocations
			if {[llength $pclist] > 0} {
				lappend ctable {*}$pclist
				lappend mainctable {*}$plist
				set clastpos $start
			} elseif {[llength $ctable] && ([expr {$start-$clastpos}] > 50)} {
				set temp [svwindow_checktrans $table $ctable $mode]
				lappend result {*}$temp
				set mainctable {}
				set ctable {}
			}
			set result [lsort -integer -index $end1pos $result]
#putsvars result
			foreach l $result {
				puts $o \t[join $l \t]
			}
			flush $o
		}
		set pstart $start
		set plist $list
		set pclist $clist
		set prlist $rlist
		set list {}
		set clist {}
		set rlist {}
		incr start $windowsize
		incr end $windowsize
	}

	gzclose $trf
	gzrmtemp $trffile
	close $o
	catch {close $f}
	catch {file delete $outfile.old}
	catch {file rename -force -- $outfile $outfile.old}
	file rename -force -- $outfile.temp $outfile
	putslog "Made svfind $outfile"

}

proc cg_bam2sv {args} {
	global scriptname action
	if {[llength $args] < 2} {
		error "format is: $scriptname $action bamfile out_prefix"
	}
	bam2sv {*}$args
}

proc cg_sv2db {args} {
	global scriptname action
	if {[llength $args] < 1} {
		error "format is: $scriptname $action file ..."
	}
	sv2db $args
}

proc cg_svcompare {args} {
	global scriptname action
	if {[llength $args] != 2} {
		error "format is: $scriptname $action svfile1 svfile2"
	}
	foreach {svfile1 svfile2} $args break
	svcompare $svfile1 $svfile2
}

proc cg_svrescore {args} {
	global scriptname action
	if {[llength $args] != 1} {
		error "format is: $scriptname $action svfile"
	}
	foreach {svfile} $args break
	svrescore $svfile
}

proc cg_map2sv {args} {
	global scriptname action
	if {[llength $args] < 2} {
		error "format is: $scriptname $action file ... out_prefix"
	}
	set prefix [list_pop args]
	map2sv $args $prefix
}

proc cg_svinfo {args} {
	global scriptname action
	if {[llength $args] < 1} {
		error "format is: $scriptname $action file ..."
	}
	foreach pairfile $args {
		putslog $pairfile
		svinfo $pairfile
	}
}

proc cg_svfind {args} {
	global scriptname action
	if {[llength $args] < 2} {
		error "format is: $scriptname $action pairedfile ... trffile"
	}
	set trffile [list_pop args]
	foreach {pairfile} $args {
		putslog $pairfile
		svfind $pairfile $trffile
	}
}

proc cg_sv_cg {args} {
	if {[llength $args] < 1} {
		error "format is: cg sv_cg subcmd ..."
  subcmds are: map2sv, find, info"
  find more extensive description using \"cg process_sv -h\""
	}
	set cmd [lindex $args 0]
	set args [lrange $args 1 end]
	switch $cmd {
		map2sv {
			cg_map2sv {*}$args
		}
		info {
			cg_svinfo {*}$args
		}
		find {
			cg_svfind {*}$args
		}
		default {
			error "Unkown subcommand $cmd to sv, must be one of map2sv, info, find"
		}
	}
}
