#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require cindex
catch {package require dbi}
package require dbi_sqlite3
package require BioTcl
setlog {}

set temperature 65
set extraseq 1000
set prefmargin 30
set maxnum 1000
if {![info exists threads]} {set threads 1}

set temp 100
foreach rtype {
	partrepeat,low repeat,low 3repeat,low completerepeat,low
	partrepeat,multi repeat,multi 3repeat,multi completerepeat,multi
	no3dbsnp dbsnp
} {
	set rscore($rtype) $temp
	incr temp 100
}

proc makeprimers_primer3 {seq rstart rend {size 600} {temperature 60}} {
	global primer3 cachedir
	if {![info exists cachedir] || [file isdir $cachedir]} {
		set tempdir [tempdir]
	} else {
		set tempdir $cachedir
	}
	set f [open $tempdir/primer3input.txt w]
	puts $f "PRIMER_SEQUENCE_ID=temp"
	puts $f "SEQUENCE=$seq"
	puts $f "PRIMER_MIN_SIZE=18"
	puts $f "PRIMER_OPT_SIZE=22"
	puts $f "PRIMER_MAX_SIZE=26"
	puts $f "PRIMER_MIN_TM=[expr {$temperature-1}]"
	puts $f "PRIMER_OPT_TM=$temperature"
	puts $f "PRIMER_MAX_TM=[expr {$temperature+1}]"
	puts $f "PRIMER_SELF_ANY=6"
	puts $f "PRIMER_MAX_END_STABILITY=9"
	puts $f "PRIMER_NUM_NS_ACCEPTED=0"
	puts $f "PRIMER_SELF_END=2"
	puts $f "PRIMER_MAX_POLY_X=3"
	puts $f "PRIMER_MIN_GC=20"
	puts $f "PRIMER_FILE_FLAG=1"
	puts $f "PRIMER_LOWERCASE_MASKING=1"
	puts $f "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0"
#	puts $f "PRIMER_TM_SANTALUCIA=1"
#	puts $f "PRIMER_SALT_CORRECTIONS=1"
#	puts $f "PRIMER_SALT_CONC="
#	puts $f "PRIMER_DIVALENT_CONC="
#	puts $f "PRIMER_DNTP_CONC="
#	puts $f "PRIMER_DNA_CONC=50.0"
	puts $f "PRIMER_NUM_RETURN=5"
	puts $f "PRIMER_PRODUCT_SIZE_RANGE=75-$size"
	puts $f "TARGET=$rstart,[expr {$rend-$rstart}]"
	puts $f "="
	close $f
	set keep [pwd]
	cd $tempdir
	set result [exec primer3 -strict_tags < $tempdir/primer3input.txt]
	set table {}
	foreach {key value} [split $result =\n] {
		switch -regexp $key {
			{^PRIMER_LEFT_?[0-9]*$} {set left [lindex [split $value ,] 0]}
			{^PRIMER_RIGHT_?[0-9]*$} {set right [lindex [split $value ,] 0]}
			{^PRIMER_LEFT_?[0-9]*_SEQUENCE$} {set leftseq $value}
			{^PRIMER_RIGHT_?[0-9]*_SEQUENCE$} {set rightseq $value}
			{^PRIMER_PRODUCT_SIZE_?[0-9]*$} {
				lappend table [list $left [expr {$right-[string length $rightseq]+1}] $leftseq $rightseq]
			}
		}
	}
	set left [lrange [split [string trim [file_read temp.for]] \n] 3 end]
	set right [lrange [split [string trim [file_read temp.rev]] \n] 3 end]
	set primerlist {}
	set result {}
	set min [expr {$rend - $size}]
	set max [expr {$rstart + $size}]
	foreach list [list $left $right] strand {+ -} {
		set temp {}
		foreach line $list {
			foreach {num seq start len} $line break
			if {$strand eq "-"} {
				set start [expr {$start-$len+1}]
			}
			lset line 2 $start
			set end [expr {$start + $len - 1}]
			if {$start < $min} continue
			if {$end > $max} continue
			set line [linsert $line 3 $end $strand]
			set line [linsert $line 6 0]
			lappend line ? ? clean {}
			lappend temp $line
		}
		lappend result $temp
	}
	cd $keep
	return $result
}

# scores 
# none 0
# no3repeat 100
# partrepeat 200
# completerepeat 300
# no3dbsnp 400
# dbsnp 500
proc makeprimers_annotate {line shift} {
	global rscore
	# join [db exec {select * from ft}] \n
	set q [lindex $line end]
	set score 0
	set itemfts {}
	set code clean
	foreach {num seq start end strand len} $line break
	# repeats
	set fts [db exec {
		select start,end,name from ft
		where start <= ? and end >= ? and type = 'repeat'
	} $end $start]
	if {[llength $fts]} {
		# join $fts \n
		set rcode clean
		list_foreach {fstart fend name} $fts {
			set repeatname [lindex $name 0 0]
			if {$repeatname eq "trf"} {
				set repeattype low
			} else {
				set repeattype multi
			}
			set cfstart [expr {$fstart-$start}]
			if {$cfstart < 0} {set cfstart 0}
			set cfend [expr {$fend-$start}]
			if {$cfend >= $len} {set cfend [expr {$len-1}]}
			set pct [expr {double($cfend+$cfstart)/$len}]
			if {($pct > 0.8) && ($repeatname eq "dust") || ([string range $repeatname 0 2] eq "Alu")} {
				set repeatloc completerepeat
				set rcode completerepeat,multi
				set score 100000
			} else {
				if {$pct > 0.9} {
					set repeatloc completerepeat
				} elseif {($cfend >= [expr {$len-1}]) && ($cfstart <= [expr {$len-6}])} {
					set repeatloc 3repeat
				} elseif {[expr {double($len - ($cfend - $cfstart))/$len}] > 0.6} {
					set repeatloc partrepeat
				} else {
					set repeatloc repeat
				}
				if {$rscore($repeatloc,$repeattype) > $score} {
					set score $rscore($repeatloc,$repeattype)
					set rcode $repeatloc,$repeattype
				}
			}
			lappend itemfts [list repeat $cfstart $cfend $repeatloc,$repeattype $repeatname]
		}
		set code $rcode
	}
	# dbsnp
	set fts [db exec {
		select start,end,name from ft
		where start <= ? and end >= ? and type = 'dbsnp'
	} $end $start]
	if {[llength $fts]} {
		# join $fts \n
		set primer [list_fill $len 1]
		list_foreach {fstart fend name} $fts {
			set cfstart [expr {$fstart-$start}]
			if {$cfstart < 0} {set cfstart 0}
			set cfend [expr {$fend-$start}]
			if {$cfend >= $len} {set cfend [expr {$len-1}]}
			set primer [eval {lreplace $primer $cfstart $cfend} [list_fill [expr {$cfend-$cfstart+1}] 0]]
			lappend itemfts [list dbsnp $cfstart $cfend [lindex $name 0]]
		}
		if {[list_remdup [lrange $primer end-6 end]] == 1} {
			if {$score < $rscore(no3dbsnp)} {
				set score $rscore(no3dbsnp)
				set code no3dbsnp
			}
		} else {
			if {$score < $rscore(dbsnp)} {
				set score $rscore(dbsnp)
				set code dbsnp
			}
		}
	}
#	set score [expr {$score+$q}]
	lset line 2 [expr {$start + $shift}]
	lset line 3 [expr {$end + $shift}]
	lset line 6 $score
	lset line end-1 $code
	lset line end $itemfts
	return $line
}

# db: must point to a dir with cindex chromsome files
# pseq: is the sequence to be found
# add: will be added to each position where pseq was found in the results
# nummax: maximum number of hits, if more are encountered, an error will be generated
# result will be a list with
#   * the total number of hits
#   * a dictionary with chromosomes as keys, and a list of the positions found on that chromosome in the value
proc cindex_searchgenome {db pseq {add 0} {nummax {}} {verbose 0}} {
	global cindex_genome maxnum
	if {$nummax eq {}} {set nummax $maxnum}
	if {![info exists cindex_genome]} {
		if {$verbose} {
			putslog "loading genome database"
		}
		set cindex_genome {}
		foreach file [ssort -natural [glob $db/*]] {
			set file [file root $file]
			set chr [lindex [split $file -] end]
			if {$verbose} {
				putslog "loading chr $chr"
			}
			dict set cindex_genome $chr [cindex load $file]
		}
	}
	set results {}
	set numresults 0
	dict for {chr cindex} $cindex_genome {
		set temp [cindex locate $cindex $pseq $nummax]
		incr numresults [llength $temp]
		if {$numresults > $nummax} {error "found > $nummax hits"}
		if {[llength $temp]} {
			if {$add} {set temp [lmath_calc $temp + $add]}
			dict set results $chr $temp
		}
	}
	return [list $numresults $results]
}

proc makeprimers_cindex {name left right {db /complgen/refseq/hg18/genome_hg18.ssa} {verbose 0}} {
	global rscore cachedir a maxnum
	unset -nocomplain a
	set list [list_concat $left $right]
	set new {}
	foreach line $list {
		set pseq [lindex $line 1]
		set pseq [string toupper $pseq]
		if {[info exists a($pseq-)]} continue
		set score [lindex $line 6]
		# if {$score == 100000} continue
		set a($pseq-) {}
		set a($pseq+) {}
		lappend new $pseq
	}
	putslog "search [llength $new]"
	foreach pseq $new {
		puts -nonewline stderr .
		if {[string length $pseq] >= 15} {
			# prelen is the size that was cut off
			set prelen [expr {[string length $pseq]-14}]
			set endseq [string range $pseq end-14 end]
		} else {
			set prelen 0
			set endseq $pseq
		}
		if {![catch {
			foreach {numhits hits} [cindex_searchgenome $db $endseq $prelen {} $verbose] break
		}]} {
			set a(${pseq}+) [list $numhits $hits]
		} else {
			set a(${pseq}+) [list $maxnum {}]
		}
		if {![catch {
			foreach {numhits hits} [cindex_searchgenome $db [seq_complement $endseq] 0 {} $verbose] break
		}]} {
			set a(${pseq}-) [list $numhits $hits]
		} else {
			set a(${pseq}-) [list $maxnum {}]
		}
	}
	putslog ""
	set result {}
	foreach list [list $left $right] strand {+ -} {
		set temp {}
		foreach line $list {
			foreach {num pseq start end strand len score} $line break
			set code [lindex $line end-1]
			if {![info exists a($pseq+)] && ![info exists a($pseq-)]} {
				lappend temp $line
				continue
			}
			set numf [lindex [get a($pseq+) ""] 0]
			set numr [lindex [get a($pseq-) ""] 0]
			set num [expr {$numf+$numr}]
			if {$num >= $maxnum} {
				set score 100000
				set code completerepeat,multi
			} elseif {$num > 100} {
				if {$score < $rscore(completerepeat,multi)} {
					set score [expr {$rscore(completerepeat,multi) + $num}]
					set code completerepeat,multi
				}
			} elseif {$num > 9} {
				if {$score < $rscore(completerepeat,low)} {
					set score [expr {$rscore(completerepeat,low) + $num}]
					set code completerepeat,low
				}
			}
			lset line 6 $score
			lset line end-1 $code
			lset line end-3 $num
			lappend temp $line
		}
		lappend result $temp
	}
	return $result
}

proc makeprimers_dimers {lseq rseq} {
	set test [seq_complement [string range $rseq end-3 end]]
	if {[string first $lseq $test] != -1} {return 1}
	set test [seq_complement [string range $lseq end-3 end]]
	if {[string first $rseq $test] != -1} {return 1}
	return 0
}

proc mergedicts {d1 d2} {
	set result $d1
	dict for {key value} $d2 {
		if {![catch {dict get $d1 $key} value1]} {
			dict set result $key [list_concat $value1 $value]
		} else {
			dict set result $key $value
		}
	}
	return $result
}

proc makeprimers_numhits {aVar lp rp maxsize cchr cstart cend} {
	global maxnum
	upvar $aVar a
	foreach {temp fseq fs fe} $lp break
	foreach {temp rseq rs re} $rp break
	# set impstart [expr {$fe+30}]
	# set impend [expr {$rs-30}]
	set impstart $fe
	set impend $rs
	set esize [expr {$rs-$fe}]
	set fseq [string toupper $fseq]
	set rseq [string toupper $rseq]
	set fseq+ [get a(${fseq}+) ""]
	set rseq+ [get a(${rseq}+) ""]
	set fseq- [get a(${fseq}-) ""]
	set rseq- [get a(${rseq}-) ""]
	if {[lindex ${fseq+} 0] >= $maxnum} {return 100000}
	if {[lindex ${fseq-} 0] >= $maxnum} {return 100000}
	if {[lindex ${rseq+} 0] >= $maxnum} {return 100000}
	if {[lindex ${rseq-} 0] >= $maxnum} {return 100000}
	set fnum [expr {[lindex ${fseq-} 0] + [lindex ${rseq+} 0]}]
	set rnum [expr {[lindex ${fseq-} 0] + [lindex ${rseq-} 0]}]
	if {($fnum == 1) && ($rnum == 1)} {
		return 1
	}
	if {($fnum >= 10) && ($rnum >= 1000)} {
		return 100000
	}
	if {($fnum >= 1000) && ($rnum >= 10)} {
		return 100000
	}
	set rchr [list_union [dict keys ${fseq-}] [dict keys ${rseq-}]]
	set flist [mergedicts [lindex ${fseq+} 1] [lindex ${rseq+} 1]]
	set rlist [mergedicts [lindex ${fseq-} 1] [lindex ${rseq-} 1]]
	# I know this can be done a lot more efficient, no time for it
	set chrs [list_common [dict keys $flist] [dict keys $rlist]]
	set diffs {}
	foreach chr $chrs {
		list_foreach fend [dict get $flist $chr] {
			list_foreach rstart [dict get $rlist $chr] {
				if {$fend > $rstart} continue
				set asize [expr {$rstart-$fend}]
				if {$asize < $esize} {
					if {($chr eq $cchr) && ($fend < $impend) && ($rstart > $impstart)} {
						if {($fend < $cend) || ($rstart > $cstart)} {
							return 999999
						}
					}
				}
				lappend diffs [list $chr $asize $fend $rstart]
			}
		}
	}
	set diffs [lsort -integer [list_subindex $diffs 1]]
	set num 0
	list_foreach diff $diffs {
		if {$diff > [expr {2000+$maxsize}]} break
		incr num
	}
	set minreg [lindex $diffs 1]
	if {($num > 1) && ([llength $flist] == 1) || ([llength $rlist] == 1)} {
		return 1.5
	}
	return $num
}

proc ucsc_epcr {p1 p2} {
	global cachedir
	if {[file exists $cachedir/$p1-$p2.epcr]} {
		return [file_read $cachedir/$p1-$p2.epcr]
	}
	package require http
	set h [http::geturl "http://genome.ucsc.edu/cgi-bin/hgPcr?hgsid=147397568&org=Human&db=hg18&wp_target=genome&wp_f=$p1&wp_r=$p2&Submit=submit&wp_size=4000&wp_perfect=15&wp_good=15&boolshad.wp_flipReverse=0"]
	set data [http::data $h]
	http::cleanup $h
	regexp {<PRE>(.*)</PRE>} $data temp pre
	set result {}
	set list [regexp -inline -all {>chr([^:]+):([0-9]+)[+-]([0-9]+)[^\n]*\n([A-Za-z\n]+)} $pre]
	foreach {temp chr start end seq} $list {
		regsub -all \n $seq {} seq
		lappend result [list $chr $start $end [expr {$end-$start+1}] $seq]
	}
	file_write $cachedir/$p1-$p2.epcr $result
	return $result
}

proc makeprimers_check {filename} {
	# set filename /complgen/compar/primers.tsv
package require Tclx
signal -restart error SIGINT
	set f [open $filename]
	set line [gets $f]
	while {![eof $f]} {
		set line [gets $f]
		if {![llength $line]} break
		if {[regexp {no primers} $line]} continue
		set e [dict get [lindex [split $line \t] 0] pairhits]
		set reg [string range [lindex [split $line \t] 1] 0 end-1]
		set p1 [lindex [split $line \t] 2]
		set line [gets $f]
		set p2 [lindex [split $line \t] 2]
		set hits [ucsc_epcr $p1 $p2]
		set poss {}
		foreach hit $hits {
			lappend poss [join [lrange $hit 0 2] -]
		}
		puts "$reg\t$e\t[llength $hits]\t[join $poss ,]"
	}
	close $f
}

proc makeprimers_filterseq {fseq ftype} {
	set fts [e features 0]
	switch $ftype {
		none {set fvars 1 ; set frepeats hard}
		repeatsoft {set fvars 1 ; set frepeats soft}
		repeatends {set fvars 1 ; set frepeats ends}
		vars {set fvars 0 ; set frepeats ends}
		repeats {set fvars 1 ; set frepeats no}
		all {return $fseq}
	}
	foreach ft $fts {
		foreach {type loc descr} $ft break
		set floc [lindex $loc 0]
		set start [dict get $floc start]
		set end [dict get $floc end]
		if {$end < $start} continue
		if {[dict exists $floc acc]} continue
		set filter 0
		if {$type eq "variation" && $fvars} {
			if {![regexp dbSNP [dict get $descr db_xref]]} continue
			set fseq [string_replace $fseq $start $end [string_fill N [expr {$end-$start+1}]]]
		} elseif {$type eq "repeat_region"} {
			switch $frepeats {
				hard {
					set fseq [string_replace $fseq $start $end [string_fill N [expr {$end-$start+1}]]]
				}
				soft {
					set fseq [string_replace $fseq $start $end [string tolower [string range $fseq $start $end]]]
				}
				ends {
					incr start 12
					incr end -12
					if {$end >= $start} {
						set fseq [string_replace $fseq $start $end [string_fill N [expr {$end-$start+1}]]]
					}
				}
			}
		}
	}
	return $fseq
}

proc makeprimers_region {name maxsize prefsize temperature dbdir db extraseq minfreq} {
	global cachedir a prefmargin fg
	unset -nocomplain a
	foreach {cchr cstart cend} [split $name -] break
	set fg [genome_open $dbdir]
	set regstart [expr {$cstart-$extraseq-1}]
	set regend [expr {$cend+$extraseq}]
	set seq [string toupper [genome_get $fg $cchr $regstart $regend 1]]
	set rstart $extraseq
	set rend [expr {$extraseq+$cend-$cstart}]
	# mask snps and repeats, and find extended region
	# fts in db
	catch {db destroy}
	dbi_sqlite3 db
	db open :memory:
	db exec {create table ft (id integer primary key, type text, chromosome text, start integer, end integer, name text, freq double)}
	# get snps in region
	set dbsnpfiles [gzfiles $dbdir/var_*snp*.tsv.gz]
	set id 1
	foreach dbsnpfile $dbsnpfiles {
		set dbsnpheader [cg select -h $dbsnpfile]
		set poss [tsv_basicfields $dbsnpheader 3]
		lappend poss [lsearch $dbsnpheader name]
		lappend poss [lsearch $dbsnpheader freq]
		set temp [tabix $dbsnpfile chr$cchr $regstart $regend]
		foreach dbsnpline $temp {
			set dbsnpline [split $dbsnpline \t]
			foreach {chr start end name freq} [list_sub $dbsnpline $poss] break
			set start [expr {$start - $regstart}]
			if {$start < 0} {set start 0}
			incr end -1
			set end [expr {$end - $regstart}]
			if {[expr {$end - $start}] > 2} continue
			set freq [max {*}[split $freq ,]]
			if {$freq <= $minfreq && $minfreq >= 0} continue
			db set [list ft $id] type dbsnp name $name chromosome $cchr start $start end $end freq $freq
			incr id
		}
	}
	foreach repeatfile [gzfiles $dbdir/reg_*rmsk.tsv.gz $dbdir/reg_*simpleRepeat.tsv.gz] {
		set temp [tabix $repeatfile chr$cchr $regstart $regend]
		set tempheader [cg select -h $repeatfile]
		set poss [tsv_basicfields $tempheader 3]
		lappend poss [lsearch $tempheader name]
		foreach templine $temp {
			set templine [split $templine \t]
			foreach {chr start end name} [list_sub $templine $poss] break
			set start [expr {$start - $regstart}]
			if {$start < 0} {set start 0}
			incr end -1
			set end [expr {$end - $regstart}]
			db set [list ft $id] type repeat name $name chromosome $cchr start $start end $end
			incr id
		}
	}
	# join [db exec {select * from ft}] \n
	foreach {left right} [makeprimers_primer3 $seq $rstart $rend $maxsize $temperature] break
	if {![llength $left] || ![llength $right]} {
		return {}
	}
	set fleft {}
	foreach line $left {
		set line [makeprimers_annotate $line [expr {$cstart-$extraseq}]]
		lappend fleft $line
	}
	set numl [llength [list_find [list_subindex $fleft 15] clean]]
	set left [lsort -real -index 6 $fleft]
	set fright {}
	foreach line $right {
		set line [makeprimers_annotate $line [expr {$cstart-$extraseq}]]
		lappend fright $line
	}
	set numr [llength [list_find [list_subindex $fright 15] clean]]
	set right [lsort -real -index 6 $fright]
	if {$numl < 5} {set numl 5}
	if {$numr < 5} {set numr 5}
	# join [list_subindex $left 1] \n
	foreach {left right} [makeprimers_cindex $name $left $right $db] break
	set left [lsort -real -index 6 $left]
	set right [lsort -real -index 6 $right]
	set numl [llength [list_find [list_subindex $left 15] clean]]
	set num2 [llength [list_find [list_subindex $right 15] clean]]
	if {$numl < 5} {set numl 5}
	if {$numr < 5} {set numr 5}
	# make pairs
	set bestscore 1000000
	set besthits 1000000
	set bestnum 1000000
	set bestdbsnp 1000000
	set bestsizediff $maxsize
	set bestmargin 0
	set bestpair {}
	set num 1
	foreach lp $left {
		set lseq [lindex $lp 1]
		set tstart [lindex $lp 3]
		set lnum [lindex $lp end-3]
		if {![info exists a($lseq+)]} {puts "$lseq+ not found"; continue}
		foreach rp $right {
			if {![expr {$num%1000}]} {putsprogress $num}
			incr num
			set asize [expr {[lindex $rp 3]-[lindex $lp 2]}]
			if {$asize >= $maxsize} continue
			# putslog "[lindex $lp 2]-[lindex $rp 3] ($asize)"
			set rseq [lindex $rp 1]
			if {![info exists a($rseq+)]} {puts "$rseq+ not found"; continue}
			if {[makeprimers_dimers $lseq $rseq]} continue
			# check if this one is better
			# numhits
			set numhits [makeprimers_numhits a $lp $rp $maxsize $cchr $cstart $cend]
			if {$numhits == 0} {error "must find one hit"}
			set tend [lindex $rp 2]
			set margin [min [expr {$tend-$cend}] [expr {$cstart-$tstart}] $prefmargin]
			set score [expr {[lindex $lp 6]+[lindex $rp 6]}]
			set sizediff [expr {abs ($prefsize - $asize)}]
			set rnum [lindex $rp end-3]
			if {![isint $lnum]} {set lnum 100000}
			if {![isint $rnum]} {set rnum 100000}
			set tnum [expr {$rnum+$lnum}]
			set dbsnp [llength [list_find [list_subindex [lindex $lp end] 0] dbsnp]]
			incr dbsnp [llength [list_find [list_subindex [lindex $rp end] 0] dbsnp]]
			if {$numhits > $besthits} continue
			if {$numhits == $besthits} {
				if {$dbsnp > $bestdbsnp} continue
				if {$dbsnp == $bestdbsnp} {
					if {$margin < $bestmargin} continue
					if {$margin == $bestmargin} {
						if {$margin == $bestmargin} {
							# size
							if {$tnum > $bestnum} continue
							if {$tnum == $bestnum} {
								if {$sizediff > $bestsizediff} continue
								if {$sizediff == $bestsizediff} {
									# score
									if {$score > $bestscore} continue
								}
							}
						}
					}
				}
			}
			# This one is better
			set bestpair [list $lp $rp $score $numhits $sizediff]
			set besthits $numhits
			set bestscore $score
			set bestsizediff $sizediff
			set bestmargin $margin
			set bestnum $tnum
			set bestdbsnp $dbsnp
		}
	}
	# foreach {lp rp} $bestpair break
	# putsvars bestpair besthits bestscore bestsizediff
	return $bestpair
}

proc makeprimers {regionfile dbdir maxsize prefsize db {minfreq -1} {numthreads 1} {o stdout}} {
	global cachedir threads temperature extraseq
	set cachedir [pwd]/cache 
	set threads $numthreads
	file mkdir $cachedir
	catch {e destroy}
	catch {rename e {}}
	EmblFile new e
	puts $o [join {remark label sequence modification scale purification project pair species chromosome cyto target contig pos temperature mg} \t]
	set f [gzopen  $regionfile]
	set header [tsv_open $f]
	set order [tsv_basicfields $header 3]
	set regionlist [split [string trim [read $f]] \n]
	close $f
	foreach region $regionlist {
		foreach {cchr cstart cend} [list_sub $region $order] break
		set cchr [chr_clip $cchr]
		set name "${cchr}-${cstart}-${cend}"
		putslog $name
		set bestpair [makeprimers_region $name $maxsize $prefsize $temperature $dbdir $db $extraseq $minfreq]
		if {[llength $bestpair] == 1} {
			puts $o "$bestpair\t$name"
		} elseif {[llength $bestpair] < 2} {
			puts $o "no primers\t$name"
		} else {
			foreach {lp rp score numhits sizediff} $bestpair break
			set astart [lindex $lp 3]
			set aend [lindex $rp 2]
			set asize [expr {$aend-$astart+1}]
			set margins [expr {$cstart-$astart}]
			lappend margins [expr {$aend-$cend}]
			if {[catch {
				set ucschits [llength [ucsc_epcr [lindex $lp 1] [lindex $rp 1]]]
			}]} {set ucschits ?}
			foreach p [list $lp $rp] d {f r} m $margins {
				set seq [lindex $p 1]
				set start [lindex $p 2]
				set ptemp [lindex $p 9]
				set num [lindex $p end-3]
				set code [lindex $p end-1]
				set info [lindex $p end]
				set comment "pairhits $numhits ucschits $ucschits asize $asize hits $num m $m"
				if {$info ne ""} {lappend comment fts $info}
				puts $o $comment\t${name}$d\t$seq\t\t\t\t\t1\tHs\t$cchr\t\t\t\t$start\t$ptemp
			}
		}
		flush $o
	}
}

proc cg_makeprimers {args} {
	global scriptname action
	set len [llength $args]
	if {$len < 4 || $len > 6} {
		errorformat makeprimers
		exit 1
	}
	set threads 1
	set minfreq -1
	foreach {regionfile maxsize prefsize dbdir minfreq threads} $args break
	set db [lindex [glob $dbdir/genome_*.ssa] 0]
	makeprimers $regionfile $dbdir $maxsize $prefsize $db $minfreq $threads
}

proc makeprimers_makeregions {filteredfile maxsize {o stdout}} {
	global cachedir threads temperature extraseq
	set f [open $filteredfile]
	set line [gets $f]
	set chrpos [lsearch $line chromosome]
	set beginpos [lsearch $line begin]
	if {![isint $beginpos]} {
		set beginpos [lsearch $line start]
	}
	set endpos [lsearch $line end]
	set line [gets $f]
	set cchr [lindex $line $chrpos]
	set cstart [expr {[lindex $line $beginpos] - 20}]
	set cend [expr {[lindex $line $endpos] + 20}]
	set next [expr {$cstart + $maxsize}]
	puts $o chromosome\tstart\tend
	while {![eof $f]} {
		while {![eof $f]} {
			set line [gets $f]
			if {![llength $line]} continue
			set chr [lindex $line $chrpos]
			set fstart [expr {[lindex $line $beginpos] - 20}]
			set fend [expr {[lindex $line $endpos] + 20}]
			if {($chr eq $cchr) && ($fend < $next)} {
				set cend $fend
			} else {
				break
			}
		}
		puts $o $cchr\t$cstart\t$cend
		set cchr $chr
		set cstart $fstart
		set cend $fend
		set next [expr {$cstart + $maxsize}]
	}
}

proc cg_makeregions {args} {
	global scriptname action
	if {[llength $args] != 2} {
		puts stderr "format is: $scriptname $action selvariationfile maxsize"
		exit 1
	}
	foreach {selvariationfile maxsize} $args break
	makeprimers_makeregions $selvariationfile $maxsize
}

