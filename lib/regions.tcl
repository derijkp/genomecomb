package require Extral

# get data from channel, but use cache etc to keep the results together and in order
proc filtervarget {v} {
	global cache
	if {[llength [get cache($v,r) ""]]} {
		return [list_shift cache($v,r)]
	}
	if {![info exists cache($v)]} {
		set temp [gets $v]
		set vline [split $temp \t]
	} else {
		set vline $cache($v)
	}
	if {![llength $vline]} {return ""}
	set curid [lindex $vline 0]
	set todo {}
	while 1 {
		set id [lindex $vline 0]
		if {$id != $curid} {
			set cache($v) $vline
			break
		}
		lappend todo $vline
		if {[eof $v]} break
		set vline [split [gets $v] \t]
	}
	set todo [lsort -integer -index 3 $todo]
	set result [list [list_sub [list_shift todo] {2 3 4}]]
	foreach {temp cstart cend} [lindex $result 0] break
	list_foreach {temp temp chr start end} $todo {
		if {$end <= $cend} {
		} elseif {$start <= $cend} {
			lset result end 2 $end
			set cend $end
		} else {
			lappend result [list $chr $start $end]
			set cstart $start
			set cend $end
		}
	}
	set cache($v,r) $result
	return [list_shift cache($v,r)]
}

proc refconsregions {varfile} {
	global cache
	puts stderr "getting ref-(in)consistent regions from $varfile"
	set v [open "| grep ref- $varfile"]
	unset -nocomplain cache($v); unset -nocomplain cache($v,r);
	set vline [filtervarget $v]
	set num 0
	puts "chromosome\tbegin\tend"
	flush stdout
	if {![llength $vline]} {
		close $v
		return
	}
	foreach {curvchr curvstart curvend} $vline break
	set num 0
	while 1 {
		incr num
		if {![expr $num%100000]} {puts stderr $num}
		set vline [filtervarget $v]
		if {![llength $vline]} break
		foreach {vchr vstart vend} $vline break
		if {($vchr == $curvchr) && ($vstart <  $curvstart)} {
			error "This should not happen: not sorted"
		}
		if {($vchr == $curvchr) && ($vstart <=  $curvend)} {
			set curvend $vend
			continue
		}
		puts $curvchr\t$curvstart\t$curvend
		set curvchr $vchr
		set curvstart $vstart
		set curvend $vend
	}
	puts $curvchr\t$curvstart\t$curvend
	close $v
}

proc regsubtract {regfile1 regfile2} {
	global cache

	array set trans {X 90 Y 91 M 92}
	set f1 [opencgifile $regfile1 header]
	if {[lrange [split $header \t] 0 2] ne "chromosome begin end"} {
		error "header error in $regfile"
	}
	set f2 [opencgifile $regfile2 header]
	if {[lrange [split $header \t] 0 2] ne "chromosome begin end"} {
		error "header error in $regfile"
	}
	set num 0
	set line1 [cggets $f1]
	foreach {chr1 start1 end1} $line1 break
	set nchr1 [get trans($chr1) $chr1]
	set line2 [cggets $f2]
	foreach {chr2 start2 end2} $line2 break
	set nchr2 [get trans($chr2) $chr2]
	puts "chromosome\tbegin\tend"
	while 1 {
		incr num
		if {![expr $num%100000]} {puts stderr $num}
		while 1 {
			# putsvars chr1 chr2 start1 end1 start2 end2
			if {$nchr2 > $nchr1} break
			if {$start2 >= $end1} break
			if {$start2 == $end2} {
			} elseif {$start2 > $start1} {
				puts $chr1\t$start1\t$start2
				set start1 $end2
			} elseif {$end2 >= $start1} {
				set start1 $end2
			}
			set line2 [cggets $f2]
			if {![llength $line2]} break
			foreach {chr2 start2 end2} $line2 break
			set nchr2 [get trans($chr2) $chr2]
		}
		if {$end1 > $start1} {
			puts $chr1\t$start1\t$end1
		}
		set line1 [cggets $f1]
		if {![llength $line1]} break
		foreach {chr1 start1 end1} $line1 break
		set nchr1 [get trans($chr1) $chr1]
	}
	close $f1
	close $f2

}

proc reghisto {regfile} {
	global cache
	unset -nocomplain a
	set v [opencgifile $regfile header]
	set num 0
	while 1 {
		incr num
		if {![expr $num%100000]} {puts stderr $num}
		set vline [cggets $v]
		if {![llength $vline]} break
		foreach {vchr vstart vend} $vline break
		set len [expr {$vend-$vstart}]
		if {![info exists a($vchr,$len)]} {
			set a($vchr,$len) 1
		} else {
			incr a($vchr,$len)
		}
	}
	close $v
	unset -nocomplain tota
	set o stdout
	set names [lsort -dictionary [array names a]]
	foreach name $names {
		foreach {chr len} [split $name ,] break
		puts $o $chr\t$len\t$a($name)
		if {![info exists tota($len)]} {
			set tota($len) 1
		} else {
			incr tota($len) $a($name)
		}
	}
	set names [lsort -dictionary [array names tota]]
	set tot 0
	foreach len $names {
		puts $o total\t$len\t$tota($len)
		incr tot $tota($len)
	}
	puts $o all\tall\t$tot
	puts $o all\tallnonull\t[expr {$tot-$tota(0)}]
}

proc regjoin {regfile1 regfile2} {
	global cache

	catch {close $f1}; catch {close $f2}
	set f1 [open $regfile1]
	set poss1 [open_region $f1]
	set num 0
	set line1 [gets $f1]
	foreach {chr start end} [list_sub $line1 $poss1] break
	set nchr [chr2num $chr]
	set line1 [gets $f1]
	foreach {chr1 start1 end1} [list_sub $line1 $poss1] break
	set nchr1 [chr2num $chr1]
	if {$regfile2 ne ""} {
		set f2 [open $regfile2]
		set poss2 [open_region $f2]
		set line2 [gets $f2]
		foreach {chr2 start2 end2} [list_sub $line2 $poss2] break
		set nchr2 [chr2num $chr2]
	} else {
		set nchr2 -1
	}
	puts "chromosome\tbegin\tend"
	while 1 {
		incr num
		if {![expr $num%100000]} {puts stderr $num}
		# putsvars nchr1 start1 end1 nchr2 start2 end2
		if {($nchr1 == -1) && ($nchr2 == -1)} break
		if {($nchr2 == -1) || ($nchr1 < $nchr2) || (($nchr1 == $nchr2) && ($start1 < $start2))} {
			foreach {ntchr tchr tstart tend} [list $nchr1 $chr1 $start1 $end1] break
			while {![eof $f1]} {
				set line1 [split [gets $f1] \t]
				if {[llength $line1]} break
			}
			if {[llength $line1]} {
				foreach {chr1 start1 end1} [list_sub $line1 $poss1] break
				set nchr1 [chr2num $chr1]
				set line1 {}
			} else {
				set nchr1 -1
			}
		} else {
			foreach {ntchr tchr tstart tend} [list $nchr2 $chr2 $start2 $end2] break
			while {![eof $f1]} {
				set line2 [split [gets $f2] \t]
				if {[llength $line1]} break
			}
			if {[llength $line2]} {
				foreach {chr2 start2 end2} [list_sub $line2 $poss2] break
				set nchr2 [chr2num $chr2]
				set line2 {}
			} else {
				set nchr2 -1
			}
		}
		if {($ntchr > $nchr) || ($tstart > $end)} {
			puts $chr\t$start\t$end
			foreach {nchr chr start end} [list $ntchr $tchr $tstart $tend] break
		} elseif {$tend > $end} {
			set end $tend
			#puts "---->$chr:$start-$end"
		}
		#puts "$chr1:$start1-$end1"
	}
	puts $chr\t$start\t$end
	close $f1
	if {$regfile2 ne ""} {
		close $f2
	}
}

if 0 {

while {![eof $v]} {
	set vline [gets $v]
	foreach {vstart vend} [lrange $vline 3 4] break
	incr removed [expr {$vend-$vstart}]
}
67125976


package require Tclx
signal -restart error SIGINT
	lappend auto_path /home/peter/dev/completegenomics/lib
	package require Extral

	set regfile1 /data/db/_data_db_ucsc-selfchain.tsv
	set regfile2 {}

	set regfile /home/peter/dev/completegenomics/test/testreg.tsv
	set varfile /home/peter/dev/completegenomics/test/testvar.tsv
	set regfile /media/passport/complgen/GS00102/reg-GS000000078-ASM.tsv
	set varfile /media/passport/complgen/GS00102/var-GS000000078-ASM.tsv
refconsregions $varfile

close $v;unset cache($v);unset cache($v,r);
	set v [open "| grep ref- $varfile"]
filtervarget $v

	cd /media/passport/complgen/GS00102
	set regfile1 reg-GS000000078-ASM.tsv
	set regfile2 reg-refcons-GS000000078-ASM.tsv
	
}

