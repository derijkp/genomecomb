#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require Extral
package require cindex

proc primercheck_search {db searchseq add maxnum} {
	if {[catch {
		foreach {numhits hits} [cindex_searchgenome $db $searchseq $add $maxnum] break
	} errmsg]} {
		set numhits many
		set hits many
	}
	return [list $numhits $hits]
}

proc primercheck_merge {d1 d2} {
        set result {} 
	foreach d [list $d1 $d2] n {1 2} {
	        dict for {key value} $d {
	                if {[catch {dict get $result $key} resultvalues]} {
				set resultvalues {}
			}
			foreach v $value {
				lappend resultvalues [list $n $v]
			}
                        dict set result $key $resultvalues
	        }
	}
        return $result
}

proc primercheck_epcr {hits1 hits2 rhits1 rhits2 maxsize maxamplicons} {
	# return [list [list 10 10000 10200 1 2]]
	if {$hits1 eq "many" || $hits2 eq "many" || $rhits1 eq "many" || $rhits2 eq "many"} {
		return many
	}
	set flist [primercheck_merge $hits1 $hits2]
	set rlist [primercheck_merge $rhits1 $rhits2]
	# This can be done a lot more efficient, no time for it
	set chrs [list_common [dict keys $flist] [dict keys $rlist]]
	set results {}
	set numresults 0
	foreach chr $chrs {
		list_foreach {fsrc fend} [dict get $flist $chr] {
			list_foreach {rsrc rstart} [dict get $rlist $chr] {
				if {$fend > $rstart} continue
				set asize [expr {$rstart-$fend}]
				if {$asize > $maxsize} continue
				lappend results [list $chr $fend $rstart $asize $fsrc $rsrc]
				incr numresults
				if {$numresults > $maxamplicons} {return $results}
			}
		}
	}
	return $results
}

proc primercheck_overlappingamplicons {amplicons overlapVar} {
	upvar $overlapVar overlap
	set amplicons [ssort -natural $amplicons]
	foreach {pchr pbegin pend psize pf pr} [lindex $amplicons 0] break
	set overlap 0
	set result {}
	list_foreach {chr begin end size f r} [lrange $amplicons 1 end] {
		if {$pchr ne $chr} {
			lappend result [list $pchr $pbegin $pend $psize $pf $pr]
			set pchr $chr; set pbegin $begin; set pend $end; set psize $size; set pf $f; set pr $r
		} elseif {$begin >= $pend} {
			lappend result [list $pchr $pbegin $pend $psize $pf $pr]
			set pchr $chr; set pbegin $begin; set pend $end; set psize $size; set pf $f; set pr $r
		} else {
			set overlap 1
			set pbegin $begin
			set pf $f
			if {$end < $pend} {
				set pend $end
				set pr $r
			}
			set psize [expr {$pend-$pbegin}]
		}
	}
	lappend result [list $pchr $pbegin $pend $psize $pf $pr]
	return $result
}

proc cg_primercheck {args} {
	set maxnum 5000
	set maxsize 1000
	set maxamplicons 100
	set pos 0
	foreach {key value} $args {
		switch -- $key {
			-m - --maxnum {
				set maxnum $value
			}
			-s - --maxsize {
				set maxsize $value
			}
			-a - --maxamplicons {
				set maxamplicons $value
			}
			-- break
			default {
				break
			}
		}
		incr pos 2
	}
	set args [lrange $args $pos end]
	if {([llength $args] < 2 || [llength $args] > 3)} {
		errorformat primercheck
		exit 1
	}
	set resultfile {}
	foreach {primerfile dbdir resultfile} $args break
	#
	catch {close $f}; catch {close $fg}
	set fg [genome_open [lindex [glob $dbdir/genome_*.ifas] 0]]
	set f [open $primerfile]
	set header [tsv_open $f]
	set poss [list_cor $header {name primer1 primer2}]
	if {$resultfile ne ""} {
		set o [open $resultfile w]
	} else {
		set o stdout
	}
	puts $o [join [list chromosome begin end name \
		primer1 primer2 numamplicons amplicons \
		primer1_hits primer1_snpsmaxfreq primer1_snps \
		primer2_hits primer2_snpsmaxfreq primer2_snps \
		amplicon_fts remarks \
	] \t]
	set dbsnpfiles [gzfiles $dbdir/var_*snp*.tsv.gz]
	set dbsnpposs {}
	foreach dbsnp $dbsnpfiles {
		set dbsnpheader [cg select -h $dbsnp]
		set temp [tsv_basicfields $dbsnpheader 4]
		lappend temp {*}[list_cor $dbsnpheader {name freq valid weight func submitterCount submitters bitfields}]
		lappend dbsnpposs $temp
	}
	set db [lindex [glob $dbdir/genome_*.ssa] 0]
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set sub [list_sub $line $poss]
		putslog $sub
		foreach {name primer1 primer2} $sub break
		set primer1 [string toupper $primer1]
		set primer2 [string toupper $primer2]
		set primer(1) $primer1
		set primer(2) $primer2
		set plen(1) [string length $primer1]
		set plen(2) [string length $primer2]
		# find target amplicon (full primer sequences)
		foreach {targetchrom targetbegin targetend targetfwd targetrev} {? ? ? ? ?} break
		foreach {maxfreq(1) maxfreq(2) primersnps(1) primersnps(2) primerrep(1) primerrep(2) ampliconfts} {? ? ? ? ? ? ?} break
		set analysis {}
		foreach p {1 2} {
			foreach [list fnumhits($p,f) fhits($p,f)] [primercheck_search $db $primer($p) $plen($p) $maxnum] break
			foreach [list fnumhits($p,r) fhits($p,r)] [primercheck_search $db [seq_complement $primer($p)] 0 $maxnum] break
		}
		set amplicons [primercheck_epcr $fhits(1,f) $fhits(2,f) $fhits(1,r) $fhits(2,r) $maxsize $maxamplicons]
		set overlap 0
		if {$amplicons ne "many" && [llength $amplicons]} {
			set amplicons [primercheck_overlappingamplicons $amplicons overlap]
		}
		if {$amplicons eq "many"} {
			lappend analysis "Main target not found (one primer has > $maxnum hits)"
		} elseif {[llength $amplicons] == 0} {
			lappend analysis "Main target not found (no amplicon with full primer match)"
		} elseif {[llength $amplicons] > $maxamplicons} {
			lappend analysis "Main target not found (> $maxamplicons amplicon with full primer match)"
		} elseif {[llength $amplicons] > 1} {
			lappend analysis "Main target not found (> 1 amplicon ([llength $amplicons]) with full primer match)"
		} else {
			foreach {targetchrom targetbegin targetend targetsize targetfwd targetrev} [lindex $amplicons 0] break
			# analyse for snps in primers
			set primerpos($targetfwd,start) [expr {$targetbegin-$plen($targetfwd)}]
			set primerpos($targetfwd,end) $targetbegin
			set primerpos($targetrev,start) $targetend 
			set primerpos($targetrev,end) [expr {$targetend+$plen($targetrev)}]
			foreach p {1 2} {
				set maxfreq($p) -
				set primersnps($p) {}
				foreach snpposs $dbsnpposs dbsnp $dbsnpfiles {
					set temp [split [exec tabix $dbsnp chr$targetchrom:$primerpos($p,start)-$primerpos($p,end)] \n]
					foreach line $temp {
						set line [split $line \t]
						foreach {chrom begin end type snpname freq valid weight func submitterCount submitters bitfields} [list_sub $line $snpposs] break
						if {[isdouble $freq] && $freq > $maxfreq($p)} {
							set maxfreq($p) $freq
						}
						set descr "${snpname}($type@$chrom:$begin-${end}\;freq=${freq}\;valid=${valid}\;submitterCount=$submitterCount)"
						lappend primersnps($p) $descr
					}
				}
				
			}
			# annotate amplicon
			set ampliconfts {}
			foreach db [glob -nocomplain $dbdir/reg_*_homopolymer.tsv.gz] {
				set temp [split [exec tabix $db $targetchrom:$targetbegin-$targetend] \n]
				set type [lindex [split [file root [gzroot [file tail $db]]] _] 2]
				foreach line $temp {
					foreach {c b e base num} $line break
					lappend ampliconfts $base${num}($c:$b-$e)
				}
			}
		}
		# find unwanted amplicons (use last 15 bases of primer)
		foreach p {1 2} {
			if {$plen($p) >= 15} {
				set prelen($p) 15
				set endseq($p) [string range $primer($p) end-14 end]
			} else {
				set prelen($p) $plen($p)
				set endseq($p) $primer($p)
			}
			foreach [list numhits($p,f) hits($p,f)] [primercheck_search $db $endseq($p) $prelen($p) $maxnum] break
			foreach [list numhits($p,r) hits($p,r)] [primercheck_search $db [seq_complement $endseq($p)] 0 $maxnum] break
			if {$numhits($p,f) eq "many" || $numhits($p,r) eq "many"} {
				set numhits($p) many
				lappend analysis "primer$p has too many hits for further analysis (> $maxnum)"
			} else {
				set numhits($p) [expr {$numhits($p,f) + $numhits($p,r)}]
			}
		}
		set amplicons [primercheck_epcr $hits(1,f) $hits(2,f) $hits(1,r) $hits(2,r) $maxsize $maxamplicons]
		if {$amplicons eq "many"} {
			set numamplicons many
			set resultamplicons many
			lappend analysis "one primer has > $maxnum hits"
		} else {
			if {[llength $amplicons] > $maxamplicons} {
				lappend analysis "too many amplicons (> $maxamplicons) not all are listed"
				set resultamplicons many
			} else {
				set numamplicons [llength $amplicons]
				set resultamplicons {}
				list_foreach {chr begin end fwd rev} $amplicons {
					lappend resultamplicons $chr:$begin-$end
				}
			}
		}
		set resultline {}
		lappend resultline $targetchrom $targetbegin $targetend \
			$name $primer1 $primer2 $numamplicons $resultamplicons \
			$numhits(1) $maxfreq(1) [join $primersnps(1) " "] \
			$numhits(2) $maxfreq(2) [join $primersnps(2) " "] \
			[join $ampliconfts " "] \
			[join $analysis ","]
		puts $o [join $resultline \t]
		flush $o
	}
	close $o
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base $scriptname
	cg_primercheck {*}$argv
}
