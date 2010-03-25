
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
	set num 0
	foreach file $files {
		puts $file
		exec zcat $file | $appdir/bin/map2sv $num | $appdir/bin/distr2chr $prefix
		incr num
	}
	exec touch ${prefix}_map2sv_FINISHED
	set files [glob $prefix-*]
	foreach file $files {
		puts $file
		set f [open ${file}-paired.tsv w]
		puts $f [join {chr1 bin strand1 start1 end1 weight1 numl type chr2 strand2 start2 end2 weight2 numr dist num fnum side} \t]
		close $f
		exec sort -nk 5 $file >> ${file}-paired.tsv
		file delete $file
	}
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
		exec sqlite3 -separator \t $dbfile ".import \"[file normalize $file]\" hits"
		exec sqlite3 $dbfile "create index hits_type on hits(type)"
		exec sqlite3 $dbfile "create index hits_bin on hits(bin)"
		svdbinfo $dbfile
	}
}

proc svinfo {pairfile} {
	svhisto $pairfile
	set list [list_remove [split [file_read [file root $pairfile].disthisto] \n] {}]
	list_shift list
	list_pop list
	set total [lmath_sum [list_subindex $list 1]]
	set list [lsort -integer -index 1 $list]
	set mode [lindex $list end 0]
	if {$mode == -1} {
		set mode [lindex $list end-1 0]
	}
	array set a [list_concat $list]
	set min $mode
	set max $mode
	set num $a($mode)
	set minnum [expr {round(0.95*$total)}]
	while {$num < $minnum} {
		incr min -1
		incr max
		incr num [get a($min) 0]
		incr num [get a($max) 0]
	}
	set o [open [rzroot $pairfile].numinfo w]
	puts $o key\tvalue
	puts $o mode\t$mode
	puts $o min\t$min
	puts $o max\t$max
	close $o
	putslog "finished $pairfile.numinfo"
}

proc svhisto {pairfile} {
	set out [file root [rzroot $pairfile]].disthisto
	set f [rzopen $pairfile]
	set header [gets $f]
	set distpos [lsearch $header dist]
	set num 1
	while {![eof $f]} {
		incr num
		if {![expr $num%1000000]} {putslog $num}
		set dist [lindex [split [gets $f] \t] $distpos]
		if {![info exists a($dist)]} {
			set a($dist) 1
		} else {
 			incr a($dist)
		}
	}
	unset a()
	set f [open $out w]
	puts $f "dist\tnumber"
	set total 0
	foreach dist [lsort -integer [array names a]] {
		puts $f $dist\t$a($dist)
		incr total $a($dist)
	}
	puts $f $total
	close $f
	# draw histo
	set tempfile [tempfile]
	file_write $tempfile [subst -nocommands {
		#set terminal postscript eps enhanced color
		set terminal png
		set output "$out.png"
		set ylabel "Number of reads"
		set xlabel "distance between paired ends"
		plot [200:500] "$out" using 1:2
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
	set f [rzopen $file]
	set svi(header) [tsv_open $f]
	catch {close $f}
	set svi(xfield) $xfield
	set svi(xfieldpos) [lsearch $svi(header) $xfield]
	if {[inlist {.rz} [file extension $file]]} {
		set indexname [file root $file].${xfield}_index
	} else {
		set indexname $file.${xfield}_index
	}
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
	set f [rzopen $pairfile $fpos]
	while {![eof $f]} {
		set fpos [tell $f]
		set line [gets $f]
		set pos [lindex [split $line \t] $xfieldpos]
		if {$pos > $start} break
	}
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

proc kde {list} {
	set list [lsort -real $list]
	set s [lindex $list 0 0]
	set e [lindex $list end 0]
	# gausian kernel, bandwidth = 9
	set kernel {0.0262 0.0338 0.0430 0.0540 0.0670 0.0822 0.0995 0.1190 0.1406 0.1640 0.1890 0.2152 0.2420 0.2687 0.2948 0.3194 0.3419 0.3614 0.3774 0.3892 0.3965 0.3989 0.3965 0.3892 0.3774 0.3614 0.3419 0.3194 0.2948 0.2687 0.2420 0.2152 0.1890 0.1640 0.1406 0.1190 0.0995 0.0822 0.0670 0.0540 0.0430 0.0338 0.0262}
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
	draw $data $s
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

proc svscore {mode type diff zyg problems gapsize num numnontrf weight patchsize b1 sd1 b2 sd2 totnum opsdiff} {
	set psdiff [expr {abs($opsdiff)}]
	if {$patchsize > 800} {
		set score 0
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
	return $score
}

proc svwindow_checkinv {table rtable minadd maxadd} {
	global infoa
	set result {}
	set distpos $infoa(distpos)
	set end1pos $infoa(end1pos)
	set mode $infoa(mode)
	set weight1pos $infoa(weight1pos)
	set weight2pos $infoa(weight2pos)
	set chr1pos $infoa(chr1pos)
	set trfpos $infoa(trfpos)
	set clustmin 6
	set step 5
	set chr [lindex [lindex $table 0] $chr1pos]
	set rtable [lsort -integer -index $distpos $rtable]
	set lists [cluster_fts $rtable $distpos 40]
	# llength [lindex $lists 0]
	foreach list $lists {
		# join [lsort -integer -index 2 $list] \n
		set list [lsort -integer -index $end1pos $list]
		set srtables [cluster_fts $list $end1pos 30]
		# llength [lindex $srtables 0]
		foreach rtable $srtables {
			if {[llength $rtable] < $clustmin} continue
			set size [expr {[lindex $rtable 0 $distpos]-$mode}]
			set cstart [lindex $rtable 0 $end1pos]
			set cend [lindex $rtable end $end1pos]
			set patchsize [expr {$cend-$cstart+1}]
			if {$patchsize <= 1} continue
			if {($cend >= $minadd) && ($cend < $maxadd)} continue
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
			set opsdiff [expr {$patchsize-$mode}]
			set psdiff [expr {abs($opsdiff)}]
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
			# linreg
			# ------
			foreach {a b sd} [linreg [list_subindex $rtable $end1pos] [list_subindex $rtable $distpos]] break
			set b [oformat $b]; set sd [oformat $sd]
			lappend result [list $chr $cstart $cend inv $size $zyg {} {} $score $num {} $weight $patchsize $b $sd {} {} $totnum [oformat $opsdiff 0]]
		}
	}
	return $result
}

proc svwindow_checktrans {table ctable minadd maxadd} {
	global infoa
	set distpos $infoa(distpos)
	set end1pos $infoa(end1pos)
	set mode $infoa(mode)
	set weight1pos $infoa(weight1pos)
	set weight2pos $infoa(weight2pos)
	set chr1pos $infoa(chr1pos)
	set clustmin 10
	set step 5
	set chr [lindex [lindex $table 0] $chr1pos]
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
				set cstart [lindex $rtable 0 $end1pos]
				set cend [lindex $rtable end $end1pos]
				set patchsize [expr {$cend-$cstart+1}]
				if {$patchsize <= 1} continue
				set opsdiff [expr {$patchsize-$mode}]
				set psdiff [expr {abs($opsdiff)}]
				if {($cend >= $minadd) && ($cend < $maxadd)} continue
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
				set quality 7
				set pos2 [lindex $rtable 0 $start2pos]
				if {$weight > 20} {
					set quality 5
				} elseif {$weight > 5} {
					set quality 3
				} elseif {$weight > 0} {
					set quality 2
				} else {
					set quality 1
				}
				if {$num > 100} {
					incr quality 4
				} elseif {$num > 50} {
					incr quality 2
				}
				# linreg
				# ------
				set xs [list_subindex $rtable $end1pos]
				set ys [lmath_calc [list_subindex $rtable $start2pos] - [list_subindex $rtable $end1pos]]
				foreach {temp b sd} [linreg $xs $ys] break
				set b [oformat $b]; set sd [oformat $sd]
				lappend result [list $chr $cstart $cend trans $pos2 $zyg {} $chr2 $quality $num {} $weight $patchsize $b $sd {} {} $totnum [oformat $opsdiff]]
			}
		}
	}
	if {[llength [list_remdup [list_subindex $result 7]]] > 1} {
		set pos 0
		foreach line $result {
			lset result $pos 6 many
			lset result $pos 8 0
			incr pos
		}
	}
	return $result
}

proc svwindow_checkindels {table minadd maxadd} {
	global infoa
	set result {}
	set distpos $infoa(distpos)
	set end1pos $infoa(end1pos)
	set weight1pos $infoa(weight1pos)
	set weight2pos $infoa(weight2pos)
	set chr1pos $infoa(chr1pos)
	set trfpos $infoa(trfpos)
	set clustmin 10
	set step 5
	set chr [lindex [lindex $table 0] $chr1pos]
	# extract clusters of similar gapsize by kde
	set step $infoa(step)
	set mode $infoa(mode)
	set dtable [lsort -integer -index $distpos $table]
	set clusters [kde_distcluster $dtable]
	if {[llength $clusters] <= 1} return
	set closest -1
	set bestdist 1000000000
	set num 0
	# remove baseline (normal gapsize) cluster
	list_foreach {len} $clusters {
		set dist [expr {abs($len-$mode)}]
		if {$dist < $bestdist} {
			set closest $num
			set bestdist $dist
			set bestlen $len
		}
		incr num
	}
	if {$bestdist > 100} {
		set closest -1
		set bestlen -1
	}
	# join [list_subindex $clusters {0 1}] \n
	if {$closest != -1} {
		set lmode [lindex $clusters $closest 0]
		list_pop clusters $closest
	} else {
		set lmode $mode
	}
	# analyse each remaining cluster
	# join [list_subindex $clusters {0 1}] \n
	list_foreach {len num data} $clusters {
		if {[llength $data] < 10} continue
		set data [lsort -integer -index $end1pos $data]
		# cluster by position
		# llength [cluster_fts $data]
		foreach list [cluster_fts $data] {
			if {[llength $list] < $clustmin} continue
			set cstart [lindex $list 0 2]
			set cend [lindex $list end 2]
			set patchsize [expr {$cend-$cstart+1}]
			if {$patchsize <= 1} continue
			if {($cend < $minadd) || ($cend >= $maxadd)} continue
			set gapsize [expr {round( [lmath_average [list_subindex $list $distpos]] )}]
			set problems {}
			# test heterozygosity (and check if len exists)
			# ---------------------------------------------
			set hlist {}
			foreach line $table {
				set pos [lindex $line 2]
				if {($pos >= $cstart) && ($pos <= $cend)} {
					lappend hlist $line
				}
			}
			set totnum [llength $hlist]
			set hlist [lsort -integer -index $distpos $hlist]
			set maxima [kde_distcluster $hlist]
			# join [list_subindex $maxima {0 1}] \n
			set lenexists 0
			list_foreach el $maxima {
				set diff [expr {abs($el-$len)}]
				if {($diff < 20) || ([expr {$diff/$len}] < 0.5)} {
					set lenexists 1
					break
				}
			}
			if {!$lenexists} continue
			if {[llength $maxima] > 1} {set zyg het} else {set zyg hom}
			# quality
			# -------
			set diff [expr {$gapsize - $lmode}]
			if {$diff > 0} {
				set type del
			} else {
				set type ins
				set diff [expr {-$diff}]
			}
			set num [llength $list]
			set weight1 [lmath_average [list_subindex $list $weight1pos]]
			set weight2 [lmath_average [list_subindex $list $weight2pos]]
			set weight [expr {round ([min $weight1 $weight2])}]
			if {$type eq "del"} {
				set opsdiff [expr {$patchsize-$mode}]
			} else {
				set opsdiff [expr {$patchsize+$diff-$mode}]
			}
			set psdiff [expr {abs($opsdiff)}]
			# test trf
			# --------
			set numnontrf [llength [list_find -exact [list_subindex $list $trfpos] 0]]
			if {$numnontrf < $clustmin} {
				lappend problems trfartefact
			} elseif {$numnontrf < [expr {0.2*$num}]} {
				lappend problems trfbad
			} elseif {$numnontrf < [expr {0.4*$num}]} {
				lappend problems trfpoor
			}
			# smallins
			# -------
			if {$diff < 50} {
				if {$zyg eq "hom"} {
					lappend problems msmall
				} else {
					lappend problems hsmall
				}
			} elseif {($type eq "ins") && ($diff < 80)} {
				if {$zyg eq "hom"} {
					lappend problems pdip
				} else {
					lappend problems pdip
				}
			}
			# linreg
			# ------
			set xs [list_subindex $list $end1pos]
			set ys [list_subindex $list $distpos]
			set mid [expr {([lindex $xs end]+[lindex $xs 0])/2}]
			set midpos 0
			foreach x $xs {
				if {$x >= $mid} break
				incr midpos
			}
			set xs1 [lrange $xs 0 $midpos]
			set ys1 [lrange $ys 0 $midpos]
			foreach {a b1 sd1} [linreg $xs1 $ys1] break
			set b1 [oformat $b1]; set sd1 [oformat $sd1]
			set xs2 [lrange $xs $midpos end]
			set ys2 [lrange $ys $midpos end]
			foreach {a b2 sd2} [linreg $xs2 $ys2] break
			set b2 [oformat $b2]; set sd2 [oformat $sd2]
			set score [svscore $mode $type $diff $zyg $problems $gapsize $num $numnontrf $weight $patchsize $b1 $sd1 $b2 $sd2 $totnum $opsdiff]
			# add to result
			# -------------
			lappend result [list $chr $cstart $cend $type [expr {round($diff)}] $zyg $problems $gapsize $score $num $numnontrf $weight $patchsize $b1 $sd1 $b2 $sd2 $totnum [oformat $opsdiff 0]]
		}
	}
	return $result
}

proc svloadtrf {trf chr pstart start end trfposs trflistVar trflineVar} {
	upvar $trflistVar trflist
	upvar $trflineVar trfline
	set temp {}
	list_foreach {trfchr trfstart trfend} $trflist {
		set trfchr [chr2num $trfchr]
		if {($trfchr == $chr) && ($trfend >= $start) && ($trfstart < $end)} {
			lappend temp [list $trfchr $trfstart $trfend]
		}
	}
	set trflist $temp
	while {![eof $trf]} {
		foreach {trfchr trfstart trfend} $trfline break
		set trfchr [chr2num $trfchr]
		if {$trfchr > $chr} break
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

proc sv_evaluate {file} {

	set file /complgen/sv/workGS103-9-paired-sv.tsv
	set f [open $file]
	set header [gets $f]
	set checkpos [lsearch $header check]
	set scorepos [lsearch $header quality]
	unset -nocomplain a
	while {![eof $f]} {
		set line [split [gets $f] \t]
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
lappend auto_path ~/dev/completegenomics/lib
package require Extral
cd /complgen/sv

set pairfile sv70-20-pairs.tsv
set pairfile GS103/GS103-20-paired.tsv.rz
set pairfile sv78-20-pairs.tsv
set pairfile sv79-20-pairs.tsv
set pairfile GS102/GS102-9-paired.tsv.rz
set pairfile GS103/GS103-9-paired.tsv.rz
set trffile /data/db/regdb-simple_repeats.tsv

}

proc svfind {pairfile trffile} {

	catch {close $trf}
	catch {close $o}
	catch {close $f}
	set bpairfile [file_rmrz $pairfile]
	set outfile [file root $bpairfile]-sv.tsv
	set windowsize 1000
	set lognum [expr {1000000 - 100000%$windowsize}]
	global infoa
	set f [open $bpairfile.numinfo]
	tsv_open $f
	array set infoa [split [string trim [read $f]] \n\t]
	close $f
	set infoa(step) 5
	set min $infoa(min)
	set max $infoa(max)
	set mode $infoa(mode)
	set infoa(hmode) [expr {$infoa(mode)-$infoa(mode)%$infoa(step)}]
	set infoa(hremove) [list_fill 5 [expr {$infoa(hmode)-2*$infoa(step)}] $infoa(step)]
	set f [rzopen $pairfile]
	set header [gets $f]
	set iheader {chr1 start1 end1 weight1 numl chr2 start2 end2 weight2 numr type dist}
	set poss [list_cor $header $iheader]
	array set infoa [list \
		end1pos [lsearch $iheader end1] \
		typepos [lsearch $iheader type] \
		chr1pos [lsearch $iheader chr1] \
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
	set typepos $infoa(typepos)
#check
#set outfile test-sv.tsv
#lassign {416000	417000} dbgstart dbgstop
#close $f
#svtools_aprgoto $pairfile $dbgstart
	set dir [file dir [file normalize $outfile]]
	set o [open $outfile w]
	puts $o [join {chr patchstart pos type size zyg problems gapsize/chr2 quality numreads numnontrf weight patchsize slope1 sd1 slope2 sd2 totnum psdiff} \t]
	set list {}
	set clist {}
	set rlist {}
	set pclist {}
	set prlist {}
	unset -nocomplain plist
	while {![eof $f]} {
		set line [split [gets $f] \t]
		set line [list_sub $line $poss]
		set start [lindex $line $end1pos]
		if {[isint $start]} break
	}
	set start [expr {$start-$start%$windowsize}]
	set pstart [expr {$start-$windowsize}]
	set end [expr {$start + $windowsize -1}]
	set chr [lindex $line $infoa(chr1pos)]
	# prepare trf data
	set trf [open $trffile]
	set trfposs [open_region $trf]
	set trflist {}
	if {[file exists $trffile.chr_index]} {
		set trfchrpos [split [string trim [file_read $trffile.chr_index]] \n\t]
		if {[dict exists $trfchrpos chr$chr]} {
			set fpos [dict get $trfchrpos chr$chr]
			seek $trf $fpos start
		} else {
			set fpos [dict get $trfchrpos $chr]
			seek $trf $fpos start
		}
	}
	while {![eof $trf]} {
		set trfline [get_region $trf $trfposs]
		foreach {trfchr trfstart trfend} $trfline break
		set trfchr [chr2num $trfchr]
		if {$trfchr >= $chr} break
	}
	# go over file and find sv
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
			set line [split [gets $f] \t]
			set line [list_sub $line $poss]
		}
		if {[info exists plist]} {
			# check for sv
			set result {}
			set minadd [expr {$start - $windowsize/2}]
			set maxadd [expr {$start + $windowsize/2}]
			# check for insertions and deletions
			set table [list_concat $plist $list]
			if {[llength $table] >= 20} {
				set temp [svwindow_checkindels $table $minadd $maxadd]
				lappend result {*}$temp
			}
#check
#if {$start >= $dbgstop} {error STOPPED}
			# check for inversions
			if {[expr {[llength $prlist] + [llength $rlist]}] > 7} {
				set rtable [list_concat $prlist $rlist]
				set temp [svwindow_checkinv $table $rtable $minadd $maxadd]
				lappend result {*}$temp
			}
			# check for translocations
			if {[expr {[llength $pclist] + [llength $clist]}] > 7} {
				set ctable [list_concat $pclist $clist]
				set temp [svwindow_checktrans $table $ctable $minadd $maxadd]
				lappend result {*}$temp
			}
			set result [lsort -integer -index $end1pos $result]
			foreach l $result {
				puts $o [join $l \t]
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

	close $trf
	close $o
	close $f
	putslog "finished $outfile"

}

if 0 {

package require Tclx
signal -restart error SIGINT
lappend auto_path ~/dev/completegenomics/lib
package require Extral
cd /complgen/sv
	cd /complgen/sv
	set svfile1 GS103/GS103-9-paired-sv.tsv
	set svfile2 GS102/GS102-9-paired-sv.tsv
	set outfile svcompar_test.tsv
	set outfile svcompar_GS102_GS103/svcompar_GS102_GS103-20.tsv
	set svfile1 sv79-20-pairs-sv.tsv
	set svfile2 sv78-20-pairs-sv.tsv
	set outfile svcompar_sv78_sv79-20.tsv
	set svfile1 GS103/GS103-9-paired-sv.tsv
	set svfile2 oldGS103-9-paired-sv.tsv

}

proc svcompare {svfile1 svfile2} {

	catch {close $f1}
	catch {close $f2}
	catch {close $o}
	set tempfile [tempfile]
	set o [open $tempfile w]
	set name1 [lindex [split [file tail $svfile1] -] 0]
	set name2 [lindex [split [file tail $svfile2] -] 0]
	set f1 [open $svfile1]
	set f2 [open $svfile2]
	set header1 [split [gets $f1] \t]
	set len [llength $header1]
	set end1pos [lsearch $header1 pos]
	set header2 [split [gets $f2] \t]
	set header [list_concat match sample $header1]
	foreach f $header2 {lappend header ${f}-2}
	puts $o [join $header \t]
	set poss1 [list_cor $header1 {patchstart pos type size zyg}]
	set poss2 [list_cor $header2 {patchstart pos type size zyg}]
	set list1 {}
	set list2 {}
	set line1 [split [gets $f1] \t]
	foreach {liststart1 listend1} [list_sub $line1 $poss1] break
	set line2 [split [gets $f2] \t]
	while {![eof $f1]} {
		while {![eof $f1]} {
			if {[llength $line1]} {
				foreach {start1 pos1 type1 size1 zyg1} [list_sub $line1 $poss1] break
#if {$pos1 == 135614960} {error STOPPED}
				if {$start1 > $listend1} break
				lappend list1 $line1
				if {$start1 < $liststart1} {set liststart1 $start1}
				set listend1 $pos1
			}
			set line1 [split [gets $f1] \t]
			set missing [max [expr {$len - [llength $line1]}] 0]
			if {$missing} {
				while {$missing} {lappend line1 {}; incr missing -1}
			}
		}
		if {![llength $list1]} break
		foreach tline2 $list2 {
			foreach {start2 pos2 type2 size2 zyg2} [list_sub $tline2 $poss2] break
			if {$pos2 < $liststart1} {
				puts $o df\t$name2\t[join $tline2 \t]
				set list2 [list_remove $list2 $tline2]
				flush $o
			}
		}
		while {![eof $f2]} {
			if {[llength $line2]} {
				foreach {start2 pos2} [list_sub $line2 $poss2] break
				if {$start2 >= $listend1} break
				if {$pos2 >= $liststart1} {
					lappend list2 $line2
				} else {
					puts $o df\t$name2\t[join $line2 \t]
				}
			}
			set line2 [split [gets $f2] \t]
		}
		unset -nocomplain a
		if {[llength $list2]} {
			# check lists
			set matrix {}
			set pos1 -1
			# puts [join $list1 \n]\n\n[join $list2 \n]
			foreach l1 $list1 {
				incr pos1
				foreach {start1 end1 type1 size1 zyg1} [list_sub $l1 $poss1] break
				set pos2 -1
				foreach l2 $list2 {
					incr pos2
					foreach {start2 end2 type2 size2 zyg2} [list_sub $l2 $poss2] break
					if {$type2 ne $type1} continue
					set diff [expr {abs($size2-$size1)}]
					set sizediff [min [max [expr {round(0.1*$size1)}] 30] 400]
					set overlap [expr {[overlap $start1 $end1 $start2 $end2]/min(double($end1-$start1+1),double($end2-$start2+1))}]
					set dist [expr {abs($end2 - $end1)}]
					#puts "$pos1\t$pos2\t$type1\t$type2\t$diff\t$overlap\t$dist"
					if {($diff < $sizediff) && (
						($overlap > 0.5) || (($overlap > 0) && ($dist < 60))
					)} {
						if {($zyg2 eq $zyg1)} {set s sm} else {set s mm}
						lappend matrix [list $pos1 $pos2 $diff $overlap $dist $s]
					}
				}
			}
			set matrix [lsort -integer -index 2 $matrix]
			list_foreach {pos1 pos2 diff overlap dist s} $matrix {
				if {[info exists a(1,$pos1)]} continue
				if {[info exists a(2,$pos2)]} continue
				set l1 [lindex $list1 $pos1]
				set l2 [lindex $list2 $pos2]
				puts $o $s\t$name1,$name2\t[join $l1 \t]\t[join $l2 \t]
				set a(1,$pos1) 1
				set a(2,$pos2) 1
			}
			set poss {}
			foreach n [array names a 2,*] {
				lappend poss [lindex [split $n ,] end]
			}
			set list2 [list_sub $list2 -exclude $poss]
		}
		# next block
		set pos -1
		foreach l1 $list1 {
			incr pos
			if {[info exists a(1,$pos)]} continue
			puts $o df\t$name1\t[join $l1 \t]
		}
		set list1 {}
		foreach {liststart1 listend1} [list_sub $line1 $poss1] break
	}


	foreach tline2 $list2 {
		puts $o df\t$name2\t[join $tline2 \t]
	}
	while {![eof $f2]} {
		if {![llength $line2]} {
			set line2 [split [gets $f2] \t]
			continue
		}
		puts $o df\t$name2\t[join $line2 \t]
		set line2 [split [gets $f2] \t]
	}
	flush $o
	close $f1
	close $f2
	close $o
	# cg select -s pos < $tempfile >@ stdout
	set chr [lindex [split [file tail $svfile1] -] 1]
	set outfile svcompar_[lindex [split [file tail $svfile1] -] 0]_[lindex [split [file tail $svfile2] -] 0]-$chr.tsv
	cg select -s pos < $tempfile > $outfile
	putslog "finished $outfile"

}
