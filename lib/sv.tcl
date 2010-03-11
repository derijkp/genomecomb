
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
	set o [open $pairfile.numinfo w]
	puts $o key\tvalue
	puts $o mode\t$mode
	puts $o min\t$min
	puts $o max\t$max
	close $o
}

proc svhisto {pairfile} {
	set out [file root $pairfile].disthisto
	set f [open $pairfile]
	set header [tsv_open $f]
	set distpos [lsearch $header dist]
	set num 1
	while {![eof $f]} {
		incr num
		if {![expr $num%1000000]} {puts stderr $num}
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

#proc svcorrectline {line min max} {
#	if {[inlist {n d i} [lindex $line 6]]} {
#		set dist [lindex $line 7]
#		if {$dist < $min} {
#			lset line 6 i
#		} elseif {$dist > $max} {
#			lset line 6 d
#		}
#	}
#	return $line
#}
#
#proc sv_writecov {dbfile chr} {
#	set outfile [file root $dbfile].cov
#	set chrsize [chrsize $chr]
#	array set infoa [exec sqlite3 -separator \t $dbfile {select * from numinfo}]
#	set min $infoa(min)
#	set max $infoa(max)
#	set mode $infoa(mode)
#	set f [open "| sqlite3 -separator \"\t\" $dbfile \"select strand1,start1,end1,strand2,start2,end2,type,dist from hits\""]
#	set o [open $outfile w]
#	set plist {}
#	set list {}
#	set line [split [gets $f] \t]
#	set line [svcorrectline $line $min $max]
#	set curpos [lindex $line 2]
#	for {set pos 0} {$pos < $chrsize} {incr pos 5} {
#		if {![expr {$pos%50000}]} {puts $pos}
#		set epos [expr {$pos+5}]
#		while {![eof $f] && ($curpos < $epos)} {
#			lappend list $line
#			set line [split [gets $f] \t]
#			set line [svcorrectline $line $min $max]
#			set curpos [lindex $line 2]
#		}
#		if {![llength $plist] && ![llength $list] } {
#			puts $o $pos\t0\t0\t0\t0\t0\t0\t0
#		} else {
#			set templist [list_concat $plist $list]
#			histogram [list_subindex $templist 6] typesa
#			set dists [list_remove [list_subindex $templist 7] -1 {}]
#			if {[llength $dists]} {
#				set dist [lmath_majority $dists]
#			} else {
#				set dist -1
#			}
#			puts $o $pos\t[get typesa(n) 0]\t[get typesa(s) 0]\t[get typesa(c) 0]\t[get typesa(r) 0]\t[get typesa(d) 0]\t[get typesa(i) 0]\t[format %.0f $dist]
#		}
#		if {[eof $f]} break
#		set plist $list
#		set list {}
#	}
#}
#
#proc svgroup {glist max} {
#	unset -nocomplain ga
#	set numg 0
#	set type [lindex $glist 0 7]
#	list_foreach {strand1 start1 end1 strand2 chr2 start2 end2 type dist} $glist {
#		for {set g 0} {$g < $numg} {incr g} {
#			if {$chr2 ne $ga($g,chr)} continue
#			if {[expr {$end1-[lindex $ga($g,end1) 0]}] > $max} continue
#			if {[expr {abs($dist-[lindex $ga($g,dist) 0])}] > 100} continue
#			if {($type eq "c") && [expr {abs($start2-[lindex $ga($g,start2) 0])}] > $max} continue
#			lappend ga($g,end1) $end1
#			lappend ga($g,dist) $dist
#			lappend ga($g,start2) $start2
#			break
#		}
#		if {$g == $numg} {
#			set ga($numg,chr) $chr2
#			set ga($numg,dist) $dist
#			set ga($numg,end1) $end1
#			set ga($numg,start2) $start2
#			incr numg
#		} 
#	}
#	set result {}
#	for {set g 0} {$g < $numg} {incr g} {
#		if {[llength $ga($g,end1)] == 1} continue
#		lappend result [list $ga($g,chr) $ga($g,end1) $ga($g,start2) $ga($g,dist)]
#	}
#	return $result
#}
#
#proc svout {o chr type mode max olist} {
#	if {[llength $olist] <= 1} {return {}}
#	set types [list_subindex $olist 7]
#	set poss [list_find $types $type]
#	set result 0
#	if {[llength $poss] > 2} {
#		set glist [list_sub $olist $poss]
#		set groups [svgroup $glist $max]
#		foreach g $groups {
#			foreach {chr2 end1s start2s dists} $g break
#			histogram $types a
#			if {$type == "c"} {
#				set dist -1
#			} else {
#				set dist [expr {abs($mode -[lmath_average $dists])}]
#			}
#			puts $o $chr\t[lindex $end1s end]\t$type\t$chr2\t[lindex [lsort -integer $start2s] 0]\t[format %.0f $dist]\t[llength $start2s]
#			incr result
#		}
#	}
#	return $result
#}
#
#proc svwindow_clearstart {list start} {
#	set pos 0
#	set len [llength $list]
#	while {$pos < $len} {
#		if {[lindex $list $pos 2] >= $start} break
#		incr pos
#	}
#	lrange $list $pos end
#}
#
#proc 2dhisto {table {step 10} {step2 5}} {
#	set histo {}
#	set table [lsort -integer -index 2 $table]
#	set field 2
#	set first [lindex $table 0 $field]
#	set hstart [expr {$first-$first%$step}]
#	set next [expr {$hstart + $step}]
#	set ellist {}
#	foreach line $table {
#		set pos [lindex $line $field]
#		if {$pos > $next} {
#			if {[llength $ellist]} {
#				set phisto [extract_histo $ellist end $step2]
#				list_foreach {size num} $phisto {
#					if {$num > 2} {
#						lappend histo [list $hstart $size $num]
#					}
#				}
#				set ellist {}
#			}
#			set hstart [expr {$pos-$pos%$step}]
#			set next [expr {$hstart + $step}]
#		}
#		lappend ellist $line
#	}
#	# file_write /tmp/histo2d.tsv [table2tsv $histo]
#	return $histo
#}
#
#proc sv_2dhisto {pairfile} {
#if 0 {
#	package require Tclx
#	signal -restart error SIGINT
#	lappend auto_path ~/dev/completegenomics/lib
#	package require Extral
#	cd /media/passport/complgen/sv
#	cd /complgen/sv
#	set pairfile sv78-20s.tsv
#	set pairfile sv70-20s.tsv
#	set pairfile sv79-20s.tsv
#	set outfile [file root $pairfile]-2d.sv
#	catch {close $o}
#	catch {close $f}
#	set o [open $outfile w]
#}
#	set o  stdout
#
#	set windowsize 1000
#	set lognum [expr {100*$windowsize}]
#	set f [open $pairfile.numinfo]
#	tsv_open $f
#	array set infoa [split [string trim [read $f]] \n\t]
#	close $f
#	set infoa(step) 5
#	set min $infoa(min)
#	set max $infoa(max)
#	set mode $infoa(mode)
#	set infoa(hmode) [expr {$infoa(mode)-$infoa(mode)%$infoa(step)}]
#	#set infoa(hremove) [list [expr {$infoa(hmode)-$infoa(step)}] $infoa(hmode) [expr {$infoa(hmode)+$infoa(step)}]]
#	set infoa(hremove) [list_fill 5 [expr {$infoa(hmode)-2*$infoa(step)}] $infoa(step)]
#	# set f [open "| sqlite3 -separator \"\t\" $dbfile \"select strand1,start1,end1,numl,strand2,chr2,start2,end2,numr,type,dist from hits\""]
#	set f [open $pairfile]
#	set header [tsv_open $f]
#	set poss [list_cor $header {strand1 start1 end1 numl strand2 chr2 start2 end2 numr type dist}]
#	puts $o pos\tsize\tnum
#	set list {}
#	set line [split [gets $f] \t]
#	set line [list_sub $line $poss]
#	set start [lindex $line 2]
#	set start [expr {$start-$start%$windowsize}]
#	set end [expr {$start + $windowsize -1}]
#	while {![eof $f]} {
#		if {![expr $start % $lognum]} {puts stderr $start}
#		while {![eof $f]} {
#			if {[llength $line]} {
#				set type [lindex $line 9]
#				set refpos [lindex $line 2]
#				if {$refpos >= $end} break
#				switch  $type {
#					r {
#						set size [lindex $line 10]
#						if {$size > 0} {
#							lset line 10 [expr {-$size}]
#						}
#						lappend list $line
#					}
#					s - c {
#					}
#					i - d - n {
#						lappend list $line
#					}
#					default {
#						error "unkown type $type"
#					}
#				}
#			}
#			set line [split [gets $f] \t]
#			set line [list_sub $line $poss]
#		}
#		if {[llength $list]} {
#			set l [2dhisto $list 20]
#			if {[llength $l]} {
#				puts $o [table2tsv $l]
#			}
#		}
#		set list {}
#		incr start $windowsize
#		incr end $windowsize
#	}
#
#	close $o
#	close $f
#}
#
#proc histo {sortedlist {step 5}} {
#	set result {}
#	set first [lindex $sortedlist 0]
#	set start [expr {$first-$first%$step}]
#	set next [expr {$start + $step}]
#	set line {}
#	foreach el $sortedlist {
#		if {$el > $next} {
#			lappend result [list $start [llength $line]]
#			set line {}
#			incr start $step
#			incr next $step
#		}
#		lappend line $el
#	}
#	lappend result [list $start [llength $line]]
#	return $result
#}

proc table2tsv {table} {
	set result {}
	foreach line $table {
		append result [join $line \t]\n
	}
	return $result
}

#proc list_maxima {gllist {step 5} {min 3}} {
#	set histo [histo $gllist $step]
#	set hs [list_subindex $histo 1]
#	set averagefilter {0.2 0.2 0.2 0.2 0.2}
#	set hs [lmath_filter $hs $averagefilter]
#	set hs [lmath_filter $hs $averagefilter]
#	set maxima {}
#	set pnum 0
#	set num 0
#	set pos -1
#	list_foreach {len} $histo nnum $hs {
#		if {($num > $min) && ($num >= $pnum) && ($num > $nnum)} {
#			lappend maxima [list [expr {$len-$step}] $num]
#		}
#		incr pos
#		set pnum $num
#		set num $nnum
#	}
#	return $maxima
#}

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

proc svtools_aprgoto {f start} {
	global svi
	set xfieldpos $svi(xfieldpos)
	if {$start < $svi(findex)} {set start $svi(findex)}
	set index $svi(index)
	set fpos [expr {round([lindex $index [expr {($start-$svi(findex))/10000}]])}]
	seek $f $fpos
	while {![eof $f]} {
		set fpos [tell $f]
		set line [gets $f]
		set pos [lindex [split $line \t] $xfieldpos]
		if {$pos > $start} break
	}
	seek $f $fpos
}

proc e {table} {
	[edit].editor set [table2tsv $table]
}

proc extract_histo {data field {step 5}} {
	set histo {}
	set first [lindex $data 0 $field]
	set hstart [expr {$first-$first%$step}]
	set next [expr {$hstart + $step}]
	set ellist {}
	foreach line $data {
		set pos [lindex $line $field]
		if {$pos > $next} {
			lappend histo [list $hstart [llength $ellist] $ellist]
			set ellist {}
			set keep [expr {$hstart+$step}]
			set hstart [expr {$pos-$pos%$step}]
			set next [expr {$hstart + $step}]
			set num 0
			while {$keep < $hstart} {
				lappend histo [list $keep 0 {}]
				incr keep $step
				incr num
				if {$num == 4} break
			}
		}
		lappend ellist $line
	}
	lappend histo [list $hstart [llength $ellist] $ellist]
	set num 4
	while {[incr num -1]} {
		lappend histo [list $next 0 {}]
		incr next $step
	}
	set temp {}
	set start [expr {[lindex $histo 0 0]-3*$step}]
	set num 4
	while {[incr num -1]} {
		lappend temp [list $start 0 {}]
		incr start $step
	}
	set histo [list_concat $temp $histo]
	# join [list_subindex $histo {0 1}] \n
	return $histo
}

proc pf {list} {
	set temp {}
	foreach el $list {
		lappend temp [format %0.1f $el]
	}
	return $temp
}

proc histo_maxima {histo step {min 3}} {
	global infoa
	set distpos $infoa(distpos)
	set result {}
	while 1 {
		if {[llength $histo] < 2} break
		set hs [list_subindex $histo 1]
		if {[llength $hs] > 5} {
			set averagefilter {0.1 0.2 0.4 0.2 0.1}
			set hs [lmath_filter $hs $averagefilter]
#			set hs [lmath_filter $hs $averagefilter]
#			set hs [lmath_filter $hs $averagefilter]
		}
		# join [pf $hs] \n
		set maxpos [maxpos $hs]
		set ph [lindex $hs $maxpos]
		if {$ph < $min} break
		set peak [lindex $hs $maxpos]
		set poss {}
		set elements {}
		set pos $maxpos
		set minpos [expr {$pos-1}]
		set minh $peak
		while {$pos > 0} {
			incr pos -1
			if {$pos < 0} break
			set h [lindex $hs $pos]
			if {$h < $minh} {
				set minpos $pos
				set minh $h
			} elseif {$h > [expr {$minh+0.1}]} {
				break
			}
		}
		set start $minpos
		set pos $maxpos
		set minpos [expr {$pos-1}]
		set minh $peak
		set len [llength $histo]
		while {$pos < $len} {
			incr pos
			set h [lindex $hs $pos]
			if {$h < $minh} {
				set minpos $pos
				set minh $h
			} elseif {$h > [expr {$minh+0.1}]} {
				break
			}
		}
		set end $minpos
		set elements [list_subindex [lrange $histo $start $end] 2]
		if {[llength $elements] <= 1} {
			set elements [lindex $elements 0]
		} else {
			set elements [list_concat {*}$elements]
		}
		if {![llength $elements]} break
		set gapsize [lmath_average [list_subindex $elements $distpos]]
		lappend result [list $gapsize [llength $elements] $elements]
		incr start -1
		incr end
		set histo [list_concat [lrange $histo 0 $start] [lrange $histo $end end]]
	}
	return $result
}

proc maxpos {vector} {
	set maxval [lmath_max $vector]
	set maxpos [lsearch $vector $maxval]
}

proc cluster_fts {data} {
	global infoa
	set checknum 20
	set maxdist 15
	set match 3
	set mismatch -5
	set result {}
	set end1pos $infoa(end1pos)
	while 1 {
		if {[llength $data] < 10} break
		set xs [list_subindex $data $end1pos]
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
			if {$diff < $maxdist} {incr score $match} else {set score [expr {$score + ($diff/$maxdist)*$mismatch}]}
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
		foreach diff [lrange $ds $pos end] {
			if {$diff < $maxdist} {incr score $match} else {set score [expr {$score + ($diff/$maxdist)*$mismatch}]}
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

proc svwindow_check {plist pclist prlist list clist rlist start windowsize} {
	global infoa
	set distpos $infoa(distpos)
	set end1pos $infoa(end1pos)
	set weight1pos $infoa(weight1pos)
	set weight2pos $infoa(weight2pos)
	set chr1pos $infoa(chr1pos)
	set trfpos $infoa(trfpos)
	set clustmin 10
	set step 5
	set chr [lindex [lindex $list 0] $chr1pos]
	set result {}
	# check for insertions and deletions
	set table [list_concat $plist $list]
	if {[llength $table] >= 20} {
		set minadd [expr {$start - $windowsize/2}]
		set maxadd [expr {$start + $windowsize/2}]
		# make histogram
		set step $infoa(step)
		set mode $infoa(mode)
		set dtable [lsort -integer -index $distpos $table]
		set histo [extract_histo $dtable $distpos $step]
		# join [list_subindex $histo {0 1}] \n
		set gapsizehisto [histo_maxima $histo $step 1]
		# list_subindex $gapsizehisto {0 1}
		if {[llength $gapsizehisto] > 1} {
			set closest -1
			set bestdist 1000000000
			set num 0
			list_foreach {len} $gapsizehisto {
				set dist [expr {abs($len-$mode)}]
				if {$dist < $bestdist} {
					set closest $num
					set bestdist $dist
				}
				incr num
			}
			if {$bestdist > 100} {
				set closest -1
			}
			# list_subindex $gapsizehisto {0 1}
			if {$closest != -1} {
				set lmode [lindex $gapsizehisto $closest 0]
				list_pop gapsizehisto $closest
			} else {
				set lmode $mode
			}
			array set qtrans {1 bad 2 poor 3 average 4  average 5 good 6 excellent}
			list_foreach {len num data} $gapsizehisto {
				if {[llength $data] < 10} continue
				set data [lsort -integer -index $end1pos $data]
				foreach list [cluster_fts $data] {
					if {[llength $list] < $clustmin} continue
					set cstart [lindex $list 0 2]
					set cend [lindex $list end 2]
					set patchsize [expr {$cend-$cstart+1}]
					if {($cend < $minadd) || ($cend >= $maxadd)} continue
					# test heterozygosity (and check if len exists)
					# ---------------------------------------------
					set hlist {}
					foreach line $table {
						set pos [lindex $line 2]
						if {($pos >= $cstart) && ($pos <= $cend)} {
							lappend hlist $line
						}
					}
					set histo [extract_histo $hlist $distpos 5]
					# join [list_subindex $histo {0 1}] \n
					set maxima [histo_maxima $histo $step 1]
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
					set gapsize [expr {round( [lmath_average [list_subindex $list $distpos]] )}]
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
					set quality unk
					set score 3
					if {$patchsize > 1000} {
						incr score -2
					} elseif {$type eq "del"} {
						if {($patchsize < 500) && ($patchsize > 200)} {
							incr score
						} elseif {$patchsize < 30} {
							incr score -2
						} elseif {$patchsize < 100} {
							incr score -1
						}
					} else {
						if {$patchsize < 500} {
							incr score
						}
					}
					if {$num >= 50} {
						incr score
					} elseif {$num < 16} {
						incr score -1
					}
					if {$weight >= 30} {
						incr score
					} elseif {$weight < 10} {
						incr score -1
					} elseif {$weight < 5} {
						incr score -2
					} elseif {$weight < 2} {
						if {$score > 2} {set score 2}
					}
					if {($type eq "ins") && ($score > 4) && ($diff < 100)} {set score 4}
					if {($weight < 2) && ($num <= 20)} {set score 1}
					if {$score <= 0} {
						set quality artefact
					} else {
						set quality $qtrans($score)
					}
					# test trf
					# --------
					set numnontrf [llength [list_find -exact [list_subindex $list $trfpos] 0]]
					if {$numnontrf < $clustmin} {
						set quality trfartefact
					} elseif {$numnontrf < [expr {0.2*$num}]} {
						set quality trfbad
					} elseif {$numnontrf < [expr {0.4*$num}]} {
						set quality trfpoor
					}
					# smallins
					# -------
					if {($type eq "ins") && ($diff < 100)} {
						if {$zyg eq "hom"} {
							set type msmins
						} else {
							set type hsmins
						}
					}
					if {($type eq "del") && ($diff < 50)} {
						if {$zyg eq "hom"} {
							set type msmdel
						} else {
							set type hsmdel
						}
					}
					# add to result
					# -------------
					lappend result [list $chr $cstart $cend $type [expr {round($diff)}] $zyg $gapsize $quality $num $numnontrf $weight $patchsize]
				}
			}
		}
	}
	# check for inversions
	if {[expr {[llength $prlist] + [llength $rlist]}] > 10} {
		set rtable [list_concat $prlist $rlist]
		# join [lsort -integer -index 2 $rtable] \n
		set srtables [cluster_fts $rtable]
		foreach rtable $srtables {
			if {[llength $rtable] < $clustmin} continue
			set size [expr {[lindex $rtable 0 $distpos]-$mode}]
			set cstart [lindex $rtable 0 $end1pos]
			set cend [lindex $rtable end $end1pos]
			if {($cend >= $minadd) && ($cend < $maxadd)} {
				set num [llength $rtable]
				# test heterozygosity
				set list {}
				foreach line $table {
					set pos [lindex $line $end1pos]
					if {($pos >= $cstart) && ($pos < $cend)} {
						lappend list $line
					}
				}
				if {[llength $list] > [expr {$num*0.4}]} {set zyg het} else {set zyg hom}
				# add to result
				set patchsize [expr {$cend-$cstart+1}]
				set weight1 [lmath_average [list_subindex $rtable $weight1pos]]
				set weight2 [lmath_average [list_subindex $rtable $weight2pos]]
				set weight [expr {round ([min $weight1 $weight2])}]
				set quality unk
				lappend result [list $chr $cstart $cend inv $size $zyg {} $quality $num $weight $patchsize]
			}
		}
	}
	# check for translocations
	if {[expr {[llength $pclist] + [llength $clist]}] > 10} {
		set chr2pos $infoa(chr2pos)
		set start2pos $infoa(start2pos)
		set ctable [list_concat $pclist $clist]
		unset -nocomplain a
		foreach line $ctable {
			set chr2 [lindex $line $chr2pos]
			lappend a($chr2) $line
		}
		foreach chr2 [array names a] {
			set rtable $a($chr2)
			if {[llength $rtable] < 5} continue
			# join [lsort -integer -index 2 $rtable] \n
			set srtables [cluster_fts $rtable]
			foreach rtable $srtables {
				if {[llength $rtable] < $clustmin} continue
				set size [expr {[lindex $rtable 0 $distpos]-$mode}]
				set cstart [lindex $rtable 0 $end1pos]
				set cend [lindex $rtable end $end1pos]
				if {($cend >= $minadd) && ($cend < $maxadd)} {
					set num [llength $rtable]
					# test heterozygosity
					set list {}
					foreach line $table {
						set pos [lindex $line $end1pos]
						if {($pos >= $cstart) && ($pos < $cend)} {
							lappend list $line
						}
					}
					if {[llength $list] > [expr {$num*0.4}]} {set zyg het} else {set zyg hom}
					# add to result
					set patchsize [expr {$cend-$cstart+1}]
					set weight1 [lmath_average [list_subindex $rtable $weight1pos]]
					set weight2 [lmath_average [list_subindex $rtable $weight2pos]]
					set weight [expr {round ([min $weight1 $weight2])}]
					set quality unk
					set pos2 [lindex $rtable 0 $start2pos]
					lappend result [list $chr $cstart $cend trans $pos2 $zyg $chr2 $quality $num $weight $patchsize]
				}
			}
		}
	}
	return $result
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
set trffile /data/db/regdb-simple_repeats.tsv

}

proc svfind {pairfile trffile} {

	catch {close $trf}
	catch {close $o}
	catch {close $f}
	if {[inlist {.rz} [file extension $pairfile]]} {
		set bpairfile [file root $pairfile]
	} else {
		set bpairfile $pairfile
	}
	set outfile [file root $bpairfile]-sv.tsv
#set outfile test-sv.tsv
	set windowsize 1000
	set lognum [expr {1000000 - 100000%$windowsize}]
svtools_loadindex $pairfile end1
	global infoa
	set dir [file dir [file normalize $outfile]]
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
#svtools_aprgoto $f 2811000
	set o [open $outfile w]
	puts $o [join {chr patchstart pos type size zyg gapsize/chr2 quality numreads numnontrf weight patchsize} \t]
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
seek $trf 10914705 start
set trfline [get_region $trf $trfposs]
	while {![eof $trf]} {
		set trfline [get_region $trf $trfposs]
		foreach {trfchr trfstart trfend} $trfline break
		set trfchr [chr2num $trfchr]
		if {$trfchr >= $chr} break
	}
	# go over file and find sv
	while {![eof $f]} {
		if {![expr $start % $lognum]} {puts stderr $start}
		# load trf regions
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
#error STOPPED
#if {$start == 173000} {error STOPPED}
			# check for sv
#puts $pstart-$end:\n[join $trflist \n]
			set result [svwindow_check $plist $pclist $prlist $list $clist $rlist $start $windowsize]
#if {[llength $result]} {error STOPPED}
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

}

if 0 {
	cd /complgen/sv
	set svfile1 GS102/GS102-20-paired-sv.tsv
	set svfile2 GS103/GS103-20-paired-sv.tsv
	set outfile svcompar_GS102_GS103/svcompar_GS102_GS103-20.tsv
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
	set header2 [split [gets $f2] \t]
	puts $o match\tsample\t[join $header1 \t]
	set poss1 [list_cor $header1 {patchstart pos type size zyg}]
	set poss2 [list_cor $header2 {patchstart pos type size zyg}]
	set list2 {}
	set line2 [split [gets $f2] \t]
	while {![eof $f1]} {
		set line1 [split [gets $f1] \t]
		if {![llength $line1]} continue
		foreach {start1 pos1 type1 size1 zyg1} [list_sub $line1 $poss1] break
		set min [expr {$start1 - 200}]
		set mstart1 [expr {$start1 - 10}]
		set mpos1 [expr {$pos1 + 10}]
		set sizediff [max [expr {round(0.1*$size1)}] 20]
		set minsize1 [expr {$size1 - $sizediff}]
		set maxsize1 [expr {$size1 + $sizediff}]
		set found 0
		foreach tline2 $list2 {
			foreach {start2 pos2 type2 size2 zyg2} [list_sub $tline2 $poss2] break
			if {$pos2 < $min} {
				puts $o df\t$name2\t[join $tline2 \t]
				set list2 [list_remove $list2 $tline2]
				flush $o
			} elseif {($type2 eq $type1) && ($zyg2 eq $zyg1)
					&& ($pos2 >= $mstart1) && ($pos2 < $mpos1)
					&& ($size2 >= $minsize1) && ($size2 < $maxsize1)} {
				puts $o sm\t$name1,$name2\t[join $line1 \t]\t[join $tline2 \t]
				set found 1
				set list2 [list_remove $list2 $tline2]
			}
		}
		while {![eof $f2] && !$found} {
			if {![llength $line2]} {
				set line2 [split [gets $f2] \t]
				continue
			}
			foreach {start2 pos2 type2 size2 zyg2} [list_sub $line2 $poss2] break
			if {$pos2 < $min} {
				puts $o df\t$name2\t[join $line2 \t]
				flush $o
			} elseif {$pos2 < $mstart1} {
			} elseif {$start2 >= $mpos1} {
				break
			} elseif {($type2 eq $type1) && ($zyg2 eq $zyg1) && ($size2 >= $minsize1) && ($size2 < $maxsize1)} {
				puts $o sm\t$name1,$name2\t[join $line1 \t]\t[join $line2 \t]
				set found 1
			} else {
				lappend list2 $line2
			}
			set line2 [split [gets $f2] \t]
		}
		if {!$found} {
			puts $o df\t$name1\t[join $line1 \t]
		}
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
	cg select -s pos < $tempfile >@ stdout

}

