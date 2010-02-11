
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
	exec touch [file dir prefix]/map2sv_FINISHED
	foreach file [glob $prefix-*] {
		puts $file
		set f [open ${file}-paired.tsv w]
		puts $f [join {chr1 bin strand1 start1 end1 weight1 numl type chr2 strand2start2 end2 weight2 numr dist num fnum side} \t]
		close $f
		exec sort -nk 5 $file >> ${file}-paired.tsv
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

proc svdbinfo {dbfile} {
	catch {
		exec sqlite3 $dbfile {
			create table info (
				key text primary key,
				value text)
		}
	}
	catch {
		exec sqlite3 $dbfile {
			create table numinfo (
				key text primary key,
				value integer)
		}
	}
	svhisto $dbfile
	set list [list_remove [split [file_read [file root $dbfile].disthisto] \n] {}]
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
	catch {
		exec sqlite3 $dbfile [subst {
			insert into numinfo (key,value) values ('mode',$mode)
		}]
	}
	catch {
		exec sqlite3 $dbfile [subst {
			insert into numinfo (key,value) values ('min',$min)
		}]
	}
	catch {
		exec sqlite3 $dbfile [subst {
			insert into numinfo (key,value) values ('max',$max)
		}]
	}
}

proc svhisto {dbfile} {
	set out [file root $dbfile].disthisto
	set f [opensqlite3 $dbfile {select dist from hits}]
	set num 1
	while {![eof $f]} {
		incr num
		if {![expr $num%1000000]} {puts stderr $num}
		set dist [gets $f]
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

proc svcorrectline {line min max} {
	if {[inlist {n d i} [lindex $line 6]]} {
		set dist [lindex $line 7]
		if {$dist < $min} {
			lset line 6 i
		} elseif {$dist > $max} {
			lset line 6 d
		}
	}
	return $line
}

proc sv_writecov {dbfile chr} {
	set outfile [file root $dbfile].cov
	set chrsize [chrsize $chr]
	array set infoa [exec sqlite3 -separator \t $dbfile {select * from numinfo}]
	set min $infoa(min)
	set max $infoa(max)
	set mode $infoa(mode)
	set f [open "| sqlite3 -separator \"\t\" $dbfile \"select strand1,start1,end1,strand2,start2,end2,type,dist from hits\""]
	set o [open $outfile w]
	set plist {}
	set list {}
	set line [split [gets $f] \t]
	set line [svcorrectline $line $min $max]
	set curpos [lindex $line 2]
	for {set pos 0} {$pos < $chrsize} {incr pos 5} {
		if {![expr {$pos%50000}]} {puts $pos}
		set epos [expr {$pos+5}]
		while {![eof $f] && ($curpos < $epos)} {
			lappend list $line
			set line [split [gets $f] \t]
			set line [svcorrectline $line $min $max]
			set curpos [lindex $line 2]
		}
		if {![llength $plist] && ![llength $list] } {
			puts $o $pos\t0\t0\t0\t0\t0\t0\t0
		} else {
			set templist [list_concat $plist $list]
			histogram [list_subindex $templist 6] typesa
			set dists [list_remove [list_subindex $templist 7] -1 {}]
			if {[llength $dists]} {
				set dist [lmath_majority $dists]
			} else {
				set dist -1
			}
			puts $o $pos\t[get typesa(n) 0]\t[get typesa(s) 0]\t[get typesa(c) 0]\t[get typesa(r) 0]\t[get typesa(d) 0]\t[get typesa(i) 0]\t[format %.0f $dist]
		}
		if {[eof $f]} break
		set plist $list
		set list {}
	}
}

proc svgroup {glist max} {
	unset -nocomplain ga
	set numg 0
	set type [lindex $glist 0 7]
	list_foreach {strand1 start1 end1 strand2 chr2 start2 end2 type dist} $glist {
		for {set g 0} {$g < $numg} {incr g} {
			if {$chr2 ne $ga($g,chr)} continue
			if {[expr {$end1-[lindex $ga($g,end1) 0]}] > $max} continue
			if {[expr {abs($dist-[lindex $ga($g,dist) 0])}] > 100} continue
			if {($type eq "c") && [expr {abs($start2-[lindex $ga($g,start2) 0])}] > $max} continue
			lappend ga($g,end1) $end1
			lappend ga($g,dist) $dist
			lappend ga($g,start2) $start2
			break
		}
		if {$g == $numg} {
			set ga($numg,chr) $chr2
			set ga($numg,dist) $dist
			set ga($numg,end1) $end1
			set ga($numg,start2) $start2
			incr numg
		} 
	}
	set result {}
	for {set g 0} {$g < $numg} {incr g} {
		if {[llength $ga($g,end1)] == 1} continue
		lappend result [list $ga($g,chr) $ga($g,end1) $ga($g,start2) $ga($g,dist)]
	}
	return $result
}

proc svout {o chr type mode max olist} {
	if {[llength $olist] <= 1} {return {}}
	set types [list_subindex $olist 7]
	set poss [list_find $types $type]
	set result 0
	if {[llength $poss] > 2} {
		set glist [list_sub $olist $poss]
		set groups [svgroup $glist $max]
		foreach g $groups {
			foreach {chr2 end1s start2s dists} $g break
			histogram $types a
			if {$type == "c"} {
				set dist -1
			} else {
				set dist [expr {abs($mode -[lmath_average $dists])}]
			}
			puts $o $chr\t[lindex $end1s end]\t$type\t$chr2\t[lindex [lsort -integer $start2s] 0]\t[format %.0f $dist]\t[llength $start2s]
			incr result
		}
	}
	return $result
}

proc svwindow_clearstart {list start} {
	set pos 0
	set len [llength $list]
	while {$pos < $len} {
		if {[lindex $list $pos 2] >= $start} break
		incr pos
	}
	lrange $list $pos end
}

proc histo {sortedlist {step 5}} {
	set result {}
	set first [lindex $sortedlist 0]
	set start [expr {$first-$first%$step}]
	set next [expr {$start + $step}]
	set line {}
	foreach el $sortedlist {
		if {$el > $next} {
			lappend result [list $start [llength $line]]
			set line {}
			incr start $step
			incr next $step
		}
		lappend line $el
	}
	lappend result [list $start [llength $line]]
	return $result
}

proc table2tsv {table} {
	set result {}
	foreach line $table {
		append result [join $line \t]\n
	}
	return $result
}

proc list_maxima {gllist {step 5} {min 3}} {
	set histo [histo $gllist $step]
	set hs [list_subindex $histo 1]
	set averagefilter {0.2 0.2 0.2 0.2 0.2}
	set hs [lmath_filter $hs $averagefilter]
	set hs [lmath_filter $hs $averagefilter]
	set maxima {}
	set pnum 0
	set num 0
	set pos -1
	list_foreach {len} $histo nnum $hs {
		if {($num > $min) && ($num >= $pnum) && ($num > $nnum)} {
			lappend maxima [list [expr {$len-$step}] $num]
		}
		incr pos
		set pnum $num
		set num $nnum
	}
	return $maxima
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
			set hstart [expr {$pos-$pos%$step}]
			set next [expr {$hstart + $step}]
		}
		lappend ellist $line
	}
	return $histo
}

proc extract_peaks {histo field {step 5} {min 3}} {
	# find maxima
	if {![llength $histo]} {return {}}
	set hs [list_subindex $histo 1]
	set averagefilter {0.2 0.2 0.2 0.2 0.2}
	set hs [lmath_filter $hs $averagefilter]
	set hs [lmath_filter $hs $averagefilter]
# write_file /tmp/histo.tsv [table2tsv $histo]
# 
if 0 {
	set temp {}; foreach b [list_subindex $histo {0}] h $hs {append temp $b\t$h\n}
	file_write /tmp/histo.tsv $temp
}
	# find maxima
	set maximapos {}
	set maxima {}
	set pnum 0
	set num 0
	set pos -1
	list_foreach {len ellen el} $histo nnum $hs {
		if {($num > 3) && ($num >= $pnum) && ($num > $nnum)} {
			lappend maximapos $pos
			lappend maxima [list [expr {$len-$step}] $num]
		}
		incr pos
		set pnum $num
		set num $nnum
	}
	if {[llength $maxima] <= 1} {return {}}
	# gather peaks starting from maxima
	set elimits [lmath_calc [lmath_calc [lrange $maximapos 1 end] + [lrange $maximapos 0 end-1]] / 2]
	lappend elimits [llength $histo]
	set slimit 0
	set result {}
	foreach pos $maximapos elimit $elimits {
		set len [lindex [lindex $histo $pos] 0]
		set data {}
		set i $pos
		while {$i >= $slimit} {
			set line [lindex $histo $i]
			set num [lindex $line 1]
			if {$num < $min} break
			lappend data [lindex $line end]
			incr i -1
		}
		set i [expr {$pos+1}]
		while {$i < $elimit} {
			set line [lindex $histo $i]
			set num [lindex $line 1]
			if {$num < $min} break
			lappend data [lindex $line end]
			incr i
		}
		set data [list_concat $data]
		if {[llength $data]} {
			lappend result [list $len [llength $data] $data]
		}
		set slimit $elimit
	}
	return $result
}

proc extract_blobs {histo field {step 5} {min 3}} {
	# find maxima
	if {![llength $histo]} {return {}}
	set hs [list_subindex $histo 1]
	if {[llength $histo] > 6} {
		set averagefilter {0.2 0.2 0.2 0.2 0.2}
		set hs [lmath_filter $hs $averagefilter]
		set hs [lmath_filter $hs $averagefilter]
	}
# write_file /tmp/histo.tsv [table2tsv $histo]
# 
if 0 {
	set temp {}; foreach b [list_subindex $histo {0}] h $hs {append temp $b\t$h\n}
	file_write /tmp/histo.tsv $temp
}
	# find maxima
	set maximapos {}
	set maxima {}
	set pnum 0
	set num 0
	set pos -1
	list_foreach {len ellen el} $histo nnum $hs {
		if {($num > 3) && ($num >= $pnum) && ($num > $nnum)} {
			lappend maximapos $pos
			lappend maxima [list [expr {$len-$step}] $num]
		}
		incr pos
		set pnum $num
		set num $nnum
	}
	if {[llength $maxima] == 0} {return {}}
	# gather blobs starting from maxima
	set prev 0
	set prevdata {}
	foreach pos $maximapos {
		set len [lindex [lindex $histo $pos] 0]
		set data {}
		set i $pos
		while {$i >= $prev} {
			set line [lindex $histo $i]
			set num [lindex $line 1]
			if {$num < $min} break
			lappend data [lindex $line end]
			incr i -1
		}
		if {[expr {$i-$prev}] > 2} {
			if {[llength $prevdata]} {
				lappend result [list $len [llength $prevdata] $prevdata]
			}
		} else {
			set data [list_union $prevdata $data]
		}
		set i [expr {$pos+1}]
		while {$i < $len} {
			set line [lindex $histo $i]
			set num [lindex $line 1]
			if {$num < $min} break
			lappend data [lindex $line end]
			incr i
		}
		set data [list_concat $data]
		lappend result [list $len [llength $data] $data]
	}
	return $result
}

proc 2dhisto {table {step 10} {step2 5}} {
	set histo {}
	set table [lsort -integer -index 2 $table]
	set field 2
	set first [lindex $table 0 $field]
	set hstart [expr {$first-$first%$step}]
	set next [expr {$hstart + $step}]
	set ellist {}
	foreach line $table {
		set pos [lindex $line $field]
		if {$pos > $next} {
			if {[llength $ellist]} {
				set phisto [extract_histo $ellist end $step2]
				list_foreach {size num} $phisto {
					if {$num > 2} {
						lappend histo [list $hstart $size $num]
					}
				}
				set ellist {}
			}
			set hstart [expr {$pos-$pos%$step}]
			set next [expr {$hstart + $step}]
		}
		lappend ellist $line
	}
	# file_write /tmp/histo2d.tsv [table2tsv $histo]
	return $histo
}


proc svwindow_check {plist pclist prlist list clist rlist mid} {
	global infoa
	if {![llength $plist] && ![llength $list]} {return {}}
	# make histogram
	set step $infoa(step)
	set mode $infoa(mode)
	set table [list_concat $plist $list]
	if {[llength $table] < 20} {return {}}
	set table [lsort -integer -index end $table]
	set histo [extract_histo $table end 5]
	# list_subindex $histo {0 1}
	set gapsizehisto [extract_peaks $histo end 5 3]
	if {[llength $gapsizehisto] <= 1} {return {}}
	# list_subindex $gapsizehisto {0 1}
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
		list_pop gapsizehisto $closest
	}
	set result {}
	list_foreach {len num data} $gapsizehisto {
		set data [lsort -integer -index 2 $data]
		set histo [extract_histo $data 2 5]
		set poshisto [extract_blobs $histo 2 5 3]
		# list_subindex $poshisto {0 1}
		list_foreach {elpos elnum eldata} $poshisto {
			if {$elnum < 20} continue
			set eldata [lsort -integer -index 2 $eldata]
			set elchr [lindex $eldata 0 5]
			set elstart [lindex $eldata 0 2]
			set elend [lindex $eldata end 2]
			if {$elstart > $mid} continue
			lappend result [list $elchr $elstart $elend $len $elnum]
		}
	}
	return $result
}

proc sv_2dhisto {dbfile} {
if 0 {
	package require Tclx
	signal -restart error SIGINT
	lappend auto_path ~/dev/completegenomics/lib
	package require Extral
	cd /complgen/sv
	cd /media/passport/complgen/sv
	set dbfile sv79-20s.sqlite
	set dbfile sv78-20s.sqlite
	set dbfile sv70-20s.sqlite
	set outfile [file root $dbfile]-2d.sv
	catch {close $o}
	catch {close $f}
	set o [open $outfile w]
}
	set o  stdout
	set windowsize 1000
	set lognum [expr {100*$windowsize}]
	array set infoa [exec sqlite3 -separator \t $dbfile {select * from numinfo}]
	set infoa(step) 5
	set min $infoa(min)
	set max $infoa(max)
	set mode $infoa(mode)
	set infoa(hmode) [expr {$infoa(mode)-$infoa(mode)%$infoa(step)}]
	#set infoa(hremove) [list [expr {$infoa(hmode)-$infoa(step)}] $infoa(hmode) [expr {$infoa(hmode)+$infoa(step)}]]
	set infoa(hremove) [list_fill 5 [expr {$infoa(hmode)-2*$infoa(step)}] $infoa(step)]
	set f [open "| sqlite3 -separator \"\t\" $dbfile \"select strand1,start1,end1,numl,strand2,chr2,start2,end2,numr,type,dist from hits\""]
	puts $o pos\tsize\tnum
	set list {}
	set line [split [gets $f] \t]
	set start [lindex $line 2]
	set start [expr {$start-$start%$windowsize}]
	set end [expr {$start + $windowsize -1}]
	while {![eof $f]} {
		if {![expr $start % $lognum]} {puts stderr $start}
		while {![eof $f]} {
			if {[llength $line]} {
				set type [lindex $line 9]
				set refpos [lindex $line 2]
				if {$refpos >= $end} break
				if {$type eq "c"} {
					# lset line end [expr {-[lindex $line end]}]
					lappend list $line
				} elseif {$type eq "r"} {
				} elseif {$type eq "s"} {
				} else {
					lappend list $line
				}
			}
			set line [split [gets $f] \t]
		}
		if {[llength $list]} {
			set l [2dhisto $list 20]
			if {[llength $l]} {
				puts $o [table2tsv $l]
			}
		}
		set list {}
		incr start $windowsize
		incr end $windowsize
	}
	close $o
	close $f
}

if 0 {

877850 878250 mz all
1230000 1230500
19008500 19010000

proc findsv {dbfile chr outfile} {

package require Tclx
signal -restart error SIGINT
lappend auto_path ~/dev/completegenomics/lib
package require Extral
cd /complgen/sv
cd /media/passport/complgen/sv
set dbfile sv79-20s.sqlite
set dbfile sv78-20s.sqlite
set outfile [file root $dbfile].sv
set chr 20
set windowsize 1000
set lognum [expr {100*$windowsize}]

catch {close $o}
catch {close $f}
	set chrsize [chrsize $chr]
	array set infoa [exec sqlite3 -separator \t $dbfile {select * from numinfo}]
	set infoa(step) 5
	set min $infoa(min)
	set max $infoa(max)
	set mode $infoa(mode)
	set infoa(hmode) [expr {$infoa(mode)-$infoa(mode)%$infoa(step)}]
	#set infoa(hremove) [list [expr {$infoa(hmode)-$infoa(step)}] $infoa(hmode) [expr {$infoa(hmode)+$infoa(step)}]]
	set infoa(hremove) [list_fill 5 [expr {$infoa(hmode)-2*$infoa(step)}] $infoa(step)]
	set f [open "| sqlite3 -separator \"\t\" $dbfile \"select strand1,start1,end1,numl,strand2,chr2,start2,end2,numr,type,dist from hits\""]
	set o [open $outfile w]
	set list {}
	set clist {}
	set rlist {}
	unset -nocomplain plist
	set line [split [gets $f] \t]
	set start [lindex $line 2]
	set start [expr {$start-$start%$windowsize}]
	set end [expr {$start + $windowsize -1}]
	while {![eof $f]} {
		if {![expr $start % $lognum]} {puts stderr $start}
		while {![eof $f]} {
			if {[llength $line]} {
				set type [lindex $line 9]
				set refpos [lindex $line 2]
				if {$refpos >= $end} break
				if {$type eq "c"} {
					lappend clist $line
				} elseif {$type eq "r"} {
					lappend rlist $line
				} elseif {$type eq "s"} {
				} else {
					lappend list $line
				}
			}
			set line [split [gets $f] \t]
		}
if {$start >= 125000} {error STOPPED}
#if {$start >= 1230000} {error STOPPED}
#if {$start >= 19009000} {error STOPPED}
# file_write /tmp/temp.tsv [join $plist \n]\n*****\n[join $list \n]
# join $list \n
		if {[info exists plist]} {
			if {[llength $plist] && [llength $list]} {
				set result [svwindow_check $plist $pclist $prlist $list $clist $rlist $start]
				# if {[llength $result]} {error STOPPED}
				foreach l $result {
					puts $o [join $l \t]
				}
				flush $o
			}
		}
		set plist $list
		set pclist $clist
		set prlist $rlist
		set list {}
		set clist {}
		set rlist {}
		incr start $windowsize
		incr end $windowsize
	}

	close $o
	close $f

}

lappend auto_path ~/dev/completegenomics/lib
package require Extral

cd /complgen/GS00103/MAP/GS000005015-MAP/
set libfile [glob lib_DNB_*.tsv]
set file mapping.tsv.gz
set rfile reads.tsv.gz

# set mategap [getmategap $libfile]


zcat mapping.tsv.gz | ~/dev/completegenomics/bin/checkmap.tcl | less

zcat mapping.tsv.gz | ~/dev/completegenomics/bin/map2sv 2 | less

~/dev/completegenomics/cg map2sv mapping.tsv.gz test/sv78

cd /complgen/sv
set file sv79-20sn
set dbfile sv79-20s.sqlite

- 10736000 10736037 + 10737078 10737117 r 1041
18539168 18539207 - 18543826 18543866 r 4619
26158566 26158604 - 26159062 26159100 r 458

sqlite3 -separator $'\t' sv79-20s.sqlite "select * from hits where type = 'i' and dist < 1000000"

set list [exec sqlite3 

}