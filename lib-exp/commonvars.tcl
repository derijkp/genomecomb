#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_commonvars args {
	# read grouping data
	if {[llength $args] != 3} {
		error "Wrong # args: should be:\ncg commonvars comparfile groupfile prefix"
	}
	foreach {comparfile groupfile prefix} $args break
	# read compar file
	set f [gzopen [gzfile $comparfile]]
	set header [tsv_open $f]
	set gposs [list_find -glob $header alleleSeq*]
	set gsamples [list_regsub {^alleleSeq[12]-} [list_sub $header $gposs] {}]
	set fams {}
	# read groupfile
	set data [file_read $groupfile]
	set data [split [string trim $data] \n]
	set gheader [split [list_pop data 0] \t]
	set poss [list_cor $gheader {Group Project Family SampleDir}]
	if {[lsearch $poss -1] != -1} {error "$groupfile should contain fields Group, Project, Family and SampleDir"}
	unset -nocomplain groupsa
	unset -nocomplain famsa
	set fnum 1
	foreach line $data {
		set line [split $line \t]
		foreach {group project family sampledir} [list_sub $line $poss] break
		if {$sampledir ni $gsamples} continue
		if {$family == ""} {
			set family $fnum
			incr fnum
		}
		set groupsa($sampledir) $group,$family
		lappend famsa($group) $family
	}
	set groups [ssort -natural [array names famsa]]
	#
	foreach sample $gsamples {
		if {[info exists groupsa($sample)]} {
			lappend fams $groupsa($sample)
		} else {
			puts "Warning: sample $sample is not used"
			lappend fams {}
		}
	}

	set rh {chromosome begin end type ref alt}
	foreach group $groups {
		set oa($group) [open $prefix$group.tsv w]
		puts $oa($group) [join {chromosome begin end type ref alt} \t]\tfreq\thfreq\ttotal
		lappend rh ${group}_freq ${group}_hfreq ${group}_tot
	}
	set o [open $prefix.tsv w]
	puts $o [join $rh \t]
	set bposs [list_cor $header {chromosome begin end type ref alt}]
	set next 10000000
	set prevchromosome {}
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set base [list_sub $line $bposs]
		foreach {chromosome begin} $base break
		if {$chromosome != $prevchromosome} {
			puts stderr $chromosome:1
			set next 10000000
			set prevchromosome $chromosome
		} elseif {$begin > $next} {
			puts stderr $chromosome-$next
			incr next 10000000
		}
		set alleles [split [lindex $base end] ,]
		if {![llength $alleles]} {set alleles [list {}]}
		set genos [list_sub $line $gposs]
		unset -nocomplain a
		unset -nocomplain ah
		unset -nocomplain tota
		foreach {geno1 geno2} $genos {fam temp} $fams {
			if {[inlist {- ?} $geno1] && [inlist {- ?} $geno2]} continue
			set a($geno1,$fam) 1
			set a($geno2,$fam) 1
			if {$geno1 eq $geno2} {
				set ah($geno1,$fam) 1
			}
			set tota($fam) 1
		}
		set result $base
		foreach group $groups {
			set temp {}
			foreach allele $alleles {
				set list [array names a $allele,$group,*]
				lappend temp [llength $list]
			}
			set temp2 {}
			foreach allele $alleles {
				set list [array names ah $allele,$group,*]
				lappend temp2 [llength $list]
			}
			set tot [llength [array names tota $group,*]]
			lappend result $temp $temp2 $tot
			set len [llength $temp]
			if {!$len || $temp == 0} continue
			if {$len > 1} {set tot [join [list_fill $len $tot] ,]}
			puts $oa($group) [join $base \t]\t[join $temp ,]\t[join $temp2 ,]\t$tot
		}
		puts $o [join $result \t]
	}
	close $o
	foreach group $groups {
		close $oa($group)
	}
}
