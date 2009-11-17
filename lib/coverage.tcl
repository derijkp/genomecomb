proc covered regfile {
	set f [opencgifile $regfile header]
	if {[lrange [split $header \t] 0 2] ne "chromosome begin end"} {
		error "header error in $regfile"
	}
	set num 0
	unset -nocomplain a
	while {![eof $f]} {
		incr num
		if {![expr $num%100000]} {puts stderr $num}
		set line [cggets $f]
		foreach {chr start end} $line break
		if {[info exists a($chr)]} {
			set a($chr) [expr {$a($chr) + $end - $start}]
		} else {
			set a($chr) [expr {$end - $start}]
		}
	}
	close $f
	puts chr\tbases
	set list [list_remove [lsort -dict [array names a]] M]
	lappend list M
	set total 0
	foreach chr $list {
		puts $chr\t$a($chr)
		incr total $a($chr)
	}
	puts ""
	puts total\t$total
}

proc filtervarget {v} {
	global cache
	if {[llength [get cache($v,r) ""]]} {
		return [list_shift cache($v,r)]
	}
	if {![info exists cache($v)]} {
		set vline [split [gets $v] \t]
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
		} else {
			lappend result [list $chr $start $end]
			set cstart $start
			set cend $end
		}
	}
	set cache($v,r) $result
	return [list_shift cache($v,r)]
}

proc filterreg {regfile varfile} {
	global cache
	array set trans {X 100 Y 101 M 102}
	set f [opencgifile $regfile header]
	if {[lrange [split $header \t] 0 2] ne "chromosome begin end"} {
		error "header error in $regfile"
	}
	set v [open "| grep ref- $varfile"]
	unset -nocomplain cache($v); unset -nocomplain cache($v,r);
	set vline [filtervarget $v]
	set num 0
	puts "chromosome\tbegin\tend"
	while {![eof $f]} {
		incr num
		if {![expr $num%100000]} {puts stderr $num}
		set line [cggets $f]
		foreach {chr start end} $line break
		set chr [get trans)($chr) $chr]
		set nstart $start
		while 1 {
			if {![llength $vline]} break
			foreach {vchr vstart vend} $vline break
			set vchr [get trans)($vchr) $vchr]
			if {$vchr > $chr} break
			if {($vchr < $chr) || ($vend == $vstart) || ($vstart < $start)} {
				set vline [filtervarget $v]
				foreach {vchr vstart vend} $vline break
				continue
			}
			if {$vstart > $end} break
			if {$vstart > $nstart} {
				puts $chr\t$nstart\t$vstart
			}
			set nstart $vend
			set vline [filtervarget $v]
		}
		if {$end >= $nstart} {
			puts $chr\t$nstart\t$end
		}
	}
	close $f
	close $v
}

proc gapsrefcons {varfile} {
	global cache
	set v [open "| grep ref- $varfile"]
	unset -nocomplain cache($v); unset -nocomplain cache($v,r);
	unset -nocomplain a
	set num 0
	while 1 {
		incr num
		if {![expr $num%100000]} {puts stderr $num}
		set vline [filtervarget $v]
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

if 0 {

while {![eof $v]} {
	set vline [gets $v]
	foreach {vstart vend} [lrange $vline 3 4] break
	incr removed [expr {$vend-$vstart}]
}
67125976


	lappend auto_path /home/peter/dev/completegenomics/lib
	package require Extral
	set regfile /media/passport/complgen/GS00102/reg-GS000000078-ASM.tsv
	set varfile /media/passport/complgen/GS00102/var-GS000000078-ASM.tsv
	set regfile /home/peter/dev/completegenomics/test/testreg.tsv
	set varfile /home/peter/dev/completegenomics/test/testvar.tsv
filterreg $regfile $varfile

close $v;unset cache($v);unset cache($v,r);
	set v [open "| grep ref- $varfile"]
filtervarget $v

}
