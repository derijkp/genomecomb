package require BioTcl
catch {package require dbi}
package require dbi_sqlite3
setlog {}

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
	if {[file isdir $cachedir]} {
		set tempdir $cachedir
	} else {
		set tempdir [tempdir]
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
			if {($cfend >= [expr {$len-1}]) && ($cfstart <= 1)} {
				set repeatloc completerepeat
				if {($repeatname eq "dust") || ([string range $repeatname 0 2] eq "Alu")} {
					set score 100000
				}
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

proc makeprimers_blast {name left right {db /data/db/blast/build36} {numl {}} {numr {}}} {
	global rscore cachedir threads a tb
	unset -nocomplain a
	if {$numl eq ""} {
		set list [list_concat $left $right]
	} else {
		set list [list_concat [lrange $left 0 $numl] [lrange $right 0 $numr]]
	}
	if {[file exists $cachedir/$name.hitsa]} {
		array set a [split [file_read $cachedir/$name.hitsa] \n]
	}
	set o [open $cachedir/$name.fas w]
	set num 0
	foreach line $list {
		set pseq [lindex $line 1]
		if {[info exists a($pseq-)]} continue
		set score [lindex $line 6]
		# if {$score == 100000} continue
		set pseq [string toupper $pseq]
		set a($pseq-) {}
		set a($pseq+) {}
		puts $o ">$pseq\n$pseq"
		incr num
	}
	close $o
	if {$num} {
		puts stderr "blast $num"
		if {![file exists $db.nin] || ![file exists $db.nsq] || ![file exists $db.nhr]} {
			error "blast database $db not found"
		}
		exec blastall -p blastn -d $db -i $cachedir/$name.fas -o $cachedir/$name.blast -W 7 -m 7 -F F -U T -v 5000 -b 5000 -e 1 -B 10000 -a $threads
		set tb [open $cachedir/$name.tabblast w]
		set f [open $cachedir/$name.blast]
		blast_parse -cmd {invoke dict {
			global a tb
			if {[dict get $dict dir] == 1} {set dir +} else {set dir -}
			set query [dict get $dict query_name]
			set len [string length $query]
			set start [dict get $dict start]
			set end [dict get $dict end]
			set identities [dict get $dict identities]
			set chr [dict get $dict hit_def]
			if {$chr eq "No definition line found"} {
				set chr [dict get $dict hit_acc]
			}
			set hstart [dict get $dict hit_start]
			set hend [dict get $dict hit_end]
			set use 0
			if {$end != $len} {
				set use 1
			} elseif {[expr {$len-$identities}] > 3} {
				set use 2
			}
			if {!$use && ([llength [get a($query$dir) ""]] < 5000)} {
				lappend a($query$dir) [list $chr $hstart $hend]
			}
			set line [list $query $use $len $start $end $dir $chr $hstart $hend $identities]
			puts $tb [join $line \t]
		}} $f
		close $tb
		file_write $cachedir/$name.hitsa [join [array get a] \n]
	}
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
			set numf [llength [get a($pseq+) ""]]
			set numr [llength [get a($pseq-) ""]]
			set num [expr {$numf+$numr}]
			if {$num >= 5000} {
				set score 100000
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

proc makeprimers_numhits {aVar lp rp size cstart cend} {
	upvar $aVar a
	foreach {temp fseq fs fe} $lp break
	foreach {temp rseq rs re} $rp break
	set esize [expr {$rs-$fe}]
	set fseq [string toupper $fseq]
	set rseq [string toupper $rseq]
	set flist [list_concat [get a($fseq+) ""] [get a($rseq+) ""]]
	set rlist [list_concat [get a($fseq-) ""] [get a($rseq-) ""]]
	if {([llength $flist] == 1) && ([llength $rlist] == 1)} {
		return 1
	}
	if {([llength $flist] >= 10) && ([llength $rlist] >= 1000)} {
		return 100000
	}
	if {([llength $flist] >= 1000) || ([llength $rlist] >= 10)} {
		return 100000
	}
	set diffs {}
	list_foreach {fchr fstart fend} $flist {
		list_foreach {rchr rstart rend} $rlist {
			if {$fchr ne $rchr} continue
			if {$fstart > $rend} continue
			set asize [expr {$rstart-$fend}]
			if {$asize < $esize} {
				if {($fend < $rs) && ($rstart > $fe)} {
					if {($fend < $cend) || ($rstart > $cstart)} {
						return 999999
					}
				}
			}
			lappend diffs [list $asize $fend $rstart]
		}
	}
	set diffs [lsort -integer -index 0 $diffs]
	set num 0
	list_foreach diff $diffs {
		if {$diff > [expr {2000+$size}]} break
		incr num
	}
	set minreg [lindex $diffs 0]
	if {($num > 1) && ([llength $flist] == 1) || ([llength $rlist] == 1)} {
		return 1.5
	}
	return $num
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

proc makeprimers_region {name size ssize temperature archive db extraseq} {
	global cachedir a
	unset -nocomplain a
	foreach {cchr cstart cend} [split $name -] break
	set emblfile $cachedir/$name.embl
	if {![file exists $emblfile]} {
		set embl [ensembl_getregion $cchr [expr {$cstart-$extraseq}] [expr {$cend+$extraseq}] -archive $archive]
		file_write $cachedir/$name.embl $embl
	}
	e open $emblfile
	set seq [e sequence 0]
	set rstart $extraseq
	set rend [expr {$extraseq+$cend-$cstart}]
	# mask snps and repeats, and find extended region
#	set estart $rstart
#	set eend $rend
#	set fts [e features 0]
#	foreach ft $fts {
#		foreach {type loc descr} $ft break
#		set floc [lindex $loc 0]
#		if {[dict exists $floc acc]} continue
#		if {$type eq "exon"} {
#			set s [dict get $floc start]
#			set e [dict get $floc end]
#			if {($s >= $eend)||($e <= $estart)} continue
#			if {$e > $eend} {
#				set eend [min $e [expr {$estart + $ssize}]]
#			}
#			if {$s < $estart} {
#				set estart [max $s [expr {$eend - $ssize}]]
#			}
#		}
#	}
	# fts in db
	catch {db destroy}
	dbi_sqlite3 db
	db open :memory:
	db exec {create table ft (id integer primary key, type text, chromosome text, start integer, end integer, name text)}
	set fts [e features 0]
	set id 1
	foreach ft $fts {
		foreach {type loc descr} $ft break
		set floc [lindex $loc 0]
		set start [dict get $floc start]
		if {[dict exists $floc complement]} {set complement 1} else {set complement 0}
		set end [dict get $floc end]
		if {$end < $start} continue
		if {[dict exists $floc acc]} continue
		set filter 0
		if {$type eq "variation"} {
			if {![regexp dbSNP [dict get $descr db_xref]]} continue
			if {[expr {$end - $start}] > 2} continue
			db set [list ft $id] type dbsnp name [dict get [lindex $ft end] db_xref] chromosome $cchr start $start end $end
			incr id
		} elseif {$type eq "repeat_region"} {
			db set [list ft $id] type repeat name [dict get [lindex $ft end] note] chromosome $cchr start $start end $end
			incr id
		}
	}
	# join [db exec {select * from ft}] \n
	foreach {left right} [makeprimers_primer3 $seq $rstart $rend $size $temperature] break
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
	while 1 {
		if {[lindex $left end 6] != 100000} break
		list_pop left
	}
	#set leftclean [list_sub $fleft [list_find [list_subindex $fleft 15] clean]]
	#set leftother [list_sub $fleft -exclude [list_find [list_subindex $fleft 15] clean]]
	#set left [list_concat [lsort -real -index 12 $leftclean] [lsort -real -index 6 $leftother]]
	set fright {}
	foreach line $right {
		set line [makeprimers_annotate $line [expr {$cstart-$extraseq}]]
		lappend fright $line
	}
	set numr [llength [list_find [list_subindex $fright 15] clean]]
	set right [lsort -real -index 6 $fright]
	while 1 {
		if {[lindex $right end 6] != 100000} break
		list_pop right
	}
	#set rightclean [list_sub $fright [list_find [list_subindex $fright 15] clean]]
	#set rightother [list_sub $fright -exclude [list_find [list_subindex $fright 15] clean]]
	#set right [list_concat [lsort -real -index 12 $rightclean] [lsort -real -index 6 $rightother]]
	if {$numl < 5} {set numl 5}
	if {$numr < 5} {set numr 5}
	foreach {left right} [makeprimers_blast $name $left $right $db {} {}] break
	# make pairs
	set bestscore 1000000
	set besthits 1000000
	set bestsizediff $size
	set bestpair {}
	set go 1
	set doleft [lrange $left 0 $numl]
	set doright [lrange $right 0 $numr]
	while 1 {
		foreach lp $doleft {
			set lseq [lindex $lp 1]
			if {![info exists a($lseq+)]} break
			foreach rp $doright {
				set asize [expr {[lindex $rp 3]-[lindex $lp 2]}]
				if {$asize >= $size} continue
				# puts stderr "[lindex $lp 2]-[lindex $rp 3] ($asize)"
				set rseq [lindex $rp 1]
				if {![info exists a($rseq+)]} break
				if {[makeprimers_dimers $lseq $rseq]} continue
				set numhits [makeprimers_numhits a $lp $rp $size $cstart $cend]
				set score [expr {[lindex $lp 6]+[lindex $rp 6]}]
				set sizediff [expr {$size - $asize}]
				if {($numhits < $besthits)
					|| (($numhits == $besthits) && ($score < $bestscore))
					|| (($numhits == $besthits) && ($score == $bestscore) && ($sizediff < $bestsizediff))
				} {
					set bestpair [list $lp $rp $score $numhits $sizediff]
					set besthits $numhits
					set bestscore $score
					set bestsizediff $sizediff
				}
			}
			if {[llength $bestpair] && ($numhits == 1)} break
		}
		if {[llength $bestpair] && ($numhits == 1)} break
		if {!$go} break
		puts stderr "blast all"
		set doleft $left
		set doright $right
		# foreach {left right} [makeprimers_blast $name $left $right $db {} {}] break
		incr go -1
	}
	# foreach {lp rp} $bestpair break
	# putsvars bestpair besthits bestscore bestsizediff
	return $bestpair
}

proc makeprimers {filteredfile {archive may2009} size {db /data/db/blast/build36} numthreads} {
	global cachedir threads
	set temperature 65
	set extraseq 1000
	set threads $numthreads
	set o stdout
	puts $o [join {remark label sequence modification scale purification project pair species chromosome cyto target contig pos temperature mg} \t]
	set cachedir [pwd]/cache
	file mkdir $cachedir
	set wiggleroom 100
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
	set cstart [expr {[lindex $line $beginpos] - 50}]
	set cend [expr {[lindex $line $endpos] + 50}]
	set ssize [expr {$size - $wiggleroom}]
	set next [expr {$cstart + $ssize}]
	catch {e destroy}
	EmblFile new e
	while {![eof $f]} {
		while {![eof $f]} {
			set line [gets $f]
			if {![llength $line]} continue
			set chr [lindex $line $chrpos]
			set fstart [expr {[lindex $line $beginpos] - 50}]
			set fend [expr {[lindex $line $endpos] + 50}]
			if {($chr eq $cchr) && ($fend < $next)} {
				set cend $fend
			} else {
				break
			}
		}
		set name $cchr-$cstart-$cend
		puts stderr $name
		set bestpair [makeprimers_region $name $size $ssize $temperature $archive $db $extraseq]
		if {[llength $bestpair] < 2} {
			puts $o "no primers\t$name"
		} else {
			foreach {lp rp score numhits sizediff} $bestpair break
			foreach p [list $lp $rp] d {f r} {
				set seq [lindex $p 1]
				set start [lindex $p 2]
				set ptemp [lindex $p 9]
				set num [lindex $p end-3]
				set code [lindex $p end-1]
				set info [lindex $p end]
				set comment "pairhits $numhits hits $num"
				if {$info ne ""} {lappend comment fts $info}
				puts $o $comment\t${name}$d\t$seq\t\t\t\t\t1\tHs\t$cchr\t\t\t\t$start\t$ptemp
			}
		}
		flush $o
		set cchr $chr
		set cstart $fstart
		set cend $fend
		set next [expr {$cstart + $ssize}]
	}
}

if 0 {

lappend auto_path ~/dev/completegenomics/lib
package require Tclx
signal -restart error SIGINT
package require Extral
cd /complgen/compar

	set filteredfile /complgen/compar/78vs79_sel.tsv
	set archive may2009
	set size 600
	set o stdout
	set cachedir [pwd]/cache
	set db /data/db/blast/build36
	set blastnum 50
	set threads 1

cd /mnt/extra/CompleteGenomics/refseq
~/bin/blast/formatdb -i BUILD.36.1.REFSEQ.FA -p F -n build36 -l formatdb.log

775 1349 GTGCTGGATTTGCGGGATGT TGGTCTGCATACAAAGTGCACAA none

5-23562903-23562929f
TCTACGTCTGCAGGGAGTGTGG CCCTGTGTGTGTCCTCTGGTGT

set name 1-12842376-12842477
set bestpair [makeprimers_region $name $size $ssize $temperature $archive $db $extraseq]

set p1 GTGCTGGATTTGCGGGATGT
set p2 TGGTCTGCATACAAAGTGCACAATC
set p1 TCTGAGCCTCTCCCTTACTTTACCC
set p2 CTGTTCCAGGTGCCCACTGA
set pos [lsearch [list_subindex $left 1] $p1]
set lp [lindex $left $pos]
set pos [lsearch [list_subindex $right 1] $p2]
set rp [lindex $right $pos]

}
