
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
		exec sort -nk 4 $file > ${file}s
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
		if {![expr $num%1000000]} {puts $num}
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

if 0 {

proc findsv {dbfile chr outfile} {

lappend auto_path ~/dev/completegenomics/lib
package require Extral
cd /complgen/sv
set dbfile sv79-20s.sqlite
set outfile [file root $dbfile].sv
set chr 20
catch {close $o}
catch {close $f}
	set chrsize [chrsize $chr]
	array set infoa [exec sqlite3 -separator \t $dbfile {select * from numinfo}]
	set min $infoa(min)
	set max $infoa(max)
	set mode $infoa(mode)
	set f [open "| sqlite3 -separator \"\t\" $dbfile \"select strand1,start1,end1,strand2,chr2,start2,end2,type,dist from hits\""]
	set o [open $outfile w]
	set plist {}
	set list {}
	array set lista {c {} r {} d {} i {}}
	array set nexta {c 0 r 0 d 0 i 0}
	array set ns {c {} r {} d {} i {}}
	while {![eof $f]} {
		set line [split [gets $f] \t]
		set type [lindex $line 7]
		set dist [lindex $line 8]
		if {[inlist {n d i} $type]} {
			if {$dist < $min} {
				set type i
				lset line 7 i
			} elseif {$dist > $max} {
				set type d
				lset line 7 d
			}
		}
		if {$type eq "n"} {
			foreach temp {c r d i} {
				if {$nexta($temp)} {
					lappend ns($temp) [list [lindex $line 2] $dist]
				}
			}
			continue
		} elseif {$type eq "s"} {
			continue
		}
		set curpos [lindex $line 2]
		if {![expr {$curpos%10000}]} {puts $curpos}
		lappend lista($type) $line
		set end1 [lindex $line 2]
		set olist $lista($type)
		set len [llength $olist]
		if {($len > 1) && (($end1 > $nexta($type)) || ([llength $ns($type)] > [expr {10*$len}]))} {
			# svout $o $chr $type $mode $max $olist
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
					set numv [llength $start2s]
					set s [lindex $end1s 0]
					set e [lindex $end1s end]
					set numn 0
					list_foreach {p d} $ns($type) {
						if {$p >= $s} {
							incr numn
						}
						if {$p > $e} break
					}
					#putsvars type numv numn ns($type)
					if {$numn < [expr {10*$numv}]} {
						puts $o $chr\t[lindex $end1s end]\t$type\t$chr2\t[lindex [lsort -integer $start2s] 0]\t[format %.0f $dist]\t$numv\t$numn
					}
					incr result
				}
			}
			set nexta($type) 0
			set ns($type) {}
			set lista($type) {}
		}
		set next [min [expr {$end1 + $max}] [lindex $line 5]]
		set ns($type,n) 0
		if {[llength $lista($type)] == 1} {
			set nexta($type) $next
		} elseif {$next > $nexta($type)} {
			set nexta($type) $next
		}
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
