package require Extral

proc open_region {f {headerVar {}}} {
	if {[string length $headerVar]} {
		upvar $headerVar header
	}
	set header [tsv_open $f]
	if {[string index $header 0] eq "#"} {
		set header [string range $header 1 end]
	}
	set tests {
		{chromosome begin end}
		{chrom start end}
		{chrom chromStart chromEnd}
		{genoName genoStart genoEnd}
		{genoName genoStart genoEnd}
		{tName tStart tEnd}
		{chr begin end}
		{chromosome start end}
		{chrom txStart txEnd}
	}
	foreach test $tests {
		set poss2 [list_cor $header $test]
		if {[lsearch $poss2 -1] == -1} {
			return $poss2
		}
	}
	error "not a region file"
}

proc get_region {f poss} {
	while 1 {
		set line [split [gets $f] \t]
		if {[llength $line]} break
	}
	set result [list_sub $line $poss]
#	set chr [lindex $result 0]
#	if {![isint $chr]} {
#		lset result 0 [chr2num $chr]
#	}
#	return $result
}

proc refconsregions {varfile} {
	putslog "getting ref-(in)consistent regions from $varfile"
	cg select -q {$varType == "ref-consistent" || $varType == "ref-inconsistent" || $varType == "no-call-rc" || $varType == "no-call-ri"} $varfile rctemp.tsv
	cg regjoin rctemp.tsv >@stdout
	file delete rctemp.tsv
}

proc nocallregions {varfile outfile} {
	putslog "getting partial no-call regions from $varfile"
	set h [cg select -h $varfile]
	if {[inlist $h allele]} {
		cg select -q {$varType == "no-call" && $allele != "all"} $varfile nctemp.tsv
	} else {
		cg select -q {$varType == "no-call" && $haplotype != "all"} $varfile nctemp.tsv
	}
	cg regjoin nctemp.tsv > $outfile
	file delete nctemp.tsv
}

proc regsubtract {regfile1 regfile2} {
	global cache
	# catch {close $f1}
	# catch {close $f2}
	set f1 [open $regfile1]
	set poss1 [open_region $f1]
	set f2 [open $regfile2]
	set poss2 [open_region $f2]
	close $f1; close $f2
	exec reg_subtract $regfile1 {*}$poss1 $regfile2 {*}$poss2 >@ stdout 2>@ stderr
}

proc reghisto {regfile} {
	global cache
	unset -nocomplain a
	set v [opencgifile $regfile header]
	set num 0
	while 1 {
		incr num
		if {![expr $num%100000]} {putslog $num}
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
	puts $o all\tallnonull\t[expr {$tot-[get tota(0) 0]}]
}

proc regjoin {regfile1 regfile2} {
	global cache
	# catch {close $f1}
	# catch {close $f2}
	set f1 [open $regfile1]
	set poss1 [open_region $f1]
	close $f1
	if {$regfile2 ne ""} {
		set f2 [open $regfile2]
		set poss2 [open_region $f2]
		close $f2
	} else {
		set poss2 {0 0 0}
	}
	# puts [list reg_join $regfile1 {*}$poss1 $regfile2 {*}$poss2]
	exec reg_join $regfile1 {*}$poss1 $regfile2 {*}$poss2 >@ stdout 2>@ stderr
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

	cd /complgen/compar_GS102_GS103
	set regfile1 reg-GS000000078-ASM.tsv
	set regfile2 reg-refcons-GS000000078-ASM.tsv
	
	set regfile1 /complgen/compar_GS102_GS103/GS102_GS103_regcommon-rc.tsv
	set regfile2 /data/db/regdb-repeatmasker.tsv

set name1 GS102
set name2 GS103
	set regfile1 ../${name2}/sreg-${name2}.tsv
	set regfile2 temp.tsv

cg regsubtract temp1 temp2 > temp3
tail -2000 GS102_GS103_regcommon-rc.tsv > temp1
tail -2000 GS102_GS103_regcommon-rc-repeat.tsv > temp2
kdiff3 temp1 temp2

}

