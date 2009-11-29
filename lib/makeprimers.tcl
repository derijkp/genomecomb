package require BioTcl
setlog {}

proc primer3 {fseq estart eend {size 600} {temperature 60}} {
	global primer3
	set tempdir [tempdir]
	set f [open $tempdir/primer3input.txt w]
	puts $f "PRIMER_SEQUENCE_ID=temp"
	puts $f "SEQUENCE=$fseq"
	puts $f "PRIMER_MIN_SIZE=18"
	puts $f "PRIMER_OPT_SIZE=22"
	puts $f "PRIMER_MAX_SIZE=26"
	puts $f "PRIMER_MIN_TM=[expr {$temperature-2}]"
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
	puts $f "PRIMER_NUM_RETURN=50"
	puts $f "PRIMER_PRODUCT_SIZE_RANGE=[expr {$size-100}]-$size 75-$size"
	puts $f "TARGET=$estart,[expr {$eend-$estart}]"
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
	cd $keep
	return $table
}

proc blast {cachedir name primerlist size {db /data/db/blast/build36}} {
	if {![llength $primerlist]} {return {}}
	if {![file exists $cachedir/$name.blast]} {
		set o [open $cachedir/$name.fas w]
		unset -nocomplain a
		list_foreach {fstart rstart fseq rseq} $primerlist {
			set fseq [string toupper $fseq]
			set rseq [string toupper $rseq]
			if {![info exists a($fseq)]} {
				puts $o ">$fseq\n$fseq"
				set a($fseq) 1
			}
			if {![info exists a($rseq)]} {
				puts $o ">$rseq\n$rseq"
				set a($rseq) 1
			}
		}
		set len [llength $primerlist]
		close $o
		if {![file exists $db.nin] || ![file exists $db.nsq] || ![file exists $db.nhr]} {
			error "blast database $db not found"
		}
		exec blastall -p blastn -d $db -i $cachedir/$name.fas -o $cachedir/$name.blast -W 7 -m 7 -F F -U T -b 1000000000 -e 1 -B 10000
	}
	set f [open $cachedir/$name.blast]
	set b [blast_parse $f]
	unset -nocomplain a
	set tb [open $cachedir/$name.tabblast w]
	foreach dict $b {
		if {[dict get $dict dir] == 1} {set dir +} else {set dir -}
		set query [dict get $dict query_name]
		set len [string length $query]
		set start [dict get $dict start]
		set end [dict get $dict end]
		set identities [dict get $dict identities]
		set chr [dict get $dict hit_def]
		set hstart [dict get $dict hit_start]
		set hend [dict get $dict hit_end]
		set use 0
		if {$end != $len} {
			set use 1
		} elseif {[expr {$len-$identities}] > 3} {
			set use 2
		}
		set line [list $query $use $len $start $end $dir $chr $hstart $hend $identities]
		if {!$use} {lappend a([lindex $line 0]$dir) [list $chr $hstart $hend]}
		puts $tb [join $line \t]
	}
	close $tb
	set numhits {}
	foreach line $primerlist {
		foreach {s e fseq rseq comment} $line break
		set fseq [string toupper $fseq]
		set rseq [string toupper $rseq]
		set flist [list_concat [get a($fseq+) ""] [get a($rseq+) ""]]
		set rlist [list_concat [get a($fseq-) ""] [get a($rseq-) ""]]
		if {([llength $flist] == 1) || ([llength $rlist] == 1)} {
#			return $line
		}
		set diffs {}
		list_foreach {fchr fstart fend} $flist {
			list_foreach {rchr rstart rend} $rlist {
				if {$fchr ne $rchr} continue
				if {$fstart > $rend} continue
				lappend diffs [expr {$rend-$fstart}]
			}
		}
		set diffs [lsort -integer $diffs]
		set sdiff [lindex $diffs 1]
		if {![isint $sdiff] || ($sdiff > [expr {2000+$size}])} {
#			return $line
		}
		set num 0
		foreach diff $diffs {
			if {$diff > [expr {2000+$size}]} break
			incr num
		}
		lappend numhits $num
	}
	set min [lmath_min $numhits]
	set pos [lsearch $numhits $min]
	set line [lindex $primerlist $pos]
	set temp [lindex $line end]
	lappend temp ${min}hit
	lset line end $temp
	return $line
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

proc makeprimers_region {name size ssize temperature archive db cachedir} {
	foreach {cchr cstart cend} [split $name -] break
	set extraseq 1000
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
	set estart $rstart
	set eend $rend
	set fts [e features 0]
	foreach ft $fts {
		foreach {type loc descr} $ft break
		set floc [lindex $loc 0]
		if {[dict exists $floc acc]} continue
		if {$type eq "exon"} {
			set s [dict get $floc start]
			set e [dict get $floc end]
			if {($s >= $eend)||($e <= $estart)} continue
			if {$e > $eend} {
				set eend [min $e [expr {$estart + $ssize}]]
			}
			if {$s < $estart} {
				set estart [max $s [expr {$eend - $ssize}]]
			}
		}
	}
	set primerlist {}
	unset -nocomplain a
	foreach ftype {none repeatsoft repeatends repeats vars all} {
		set fseq [makeprimers_filterseq $seq $ftype]
		set primers [primer3 $fseq $estart $eend $size $temperature]
		foreach line $primers {
			if {[info exists a([lrange $line 0 1])]} continue
			set a([lrange $line 0 1]) $ftype
			lappend line $ftype
			lappend primerlist $line
		}
		set primers [primer3 $fseq $rstart $rend $size $temperature]
		foreach line $primers {
			if {[info exists a([lrange $line 0 1])]} continue
			set a([lrange $line 0 1]) $ftype
			lappend line $ftype
			lappend primerlist $line
		}
	}
	# join $primerlist \n
	set blast [blast $cachedir $name $primerlist $size $db]
	if {![llength $blast]} {return {}}
	lset blast 0 [expr {$cstart-$extraseq + [lindex $blast 0]}]
	lset blast 1 [expr {$cstart-$extraseq + [lindex $blast 1]}]
	return $blast
}

proc makeprimers {filteredfile {archive may2009} size {db /data/db/blast/build36}} {
	set o stdout
	set temperature 65
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
	set cend [expr {[lindex $line $endpos] + 5}]
	set ssize [expr {$size - $wiggleroom}]
	set next [expr {$cstart + $ssize}]
	catch {e destroy}
	EmblFile new e
	while {![eof $f]} {
		while {![eof $f]} {
			set line [gets $f]
			if {![llength $line]} continue
			set chr [lindex $line $chrpos]
			set fstart [expr {[lindex $line $beginpos] - 20}]
			set fend [expr {[lindex $line $endpos] + 5}]
			if {($chr eq $cchr) && ($fend < $next)} {
				set cend $fend
			} else {
				break
			}
		}
		set name $cchr-$cstart-$cend
		puts stderr $name
		set blast [makeprimers_region $name $size $ssize $temperature $archive $db $cachedir]
		if {[llength $blast] < 2} {
			puts $o "no primers\t$name"
		} else {
			foreach {start1 start2 seq1 seq2 comment} $blast break
			puts $o $comment\t${name}f\t$seq1\t\t\t\t\t1\tHs\t$cchr\t\t\t\t$start1\t$temperature
			puts $o $comment\t${name}r\t$seq2\t\t\t\t\t1\tHs\t$cchr\t\t\t\t$start2\t$temperature
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
#package require Tclx
#signal -restart error SIGINT
package require Extral
cd /complgen/compar

	set filteredfile /complgen/compar/78vs79_sel.tsv
	set archive may2009
	set size 600
	set o stdout
	set cachedir [pwd]/cache
	set db /data/db/blast/build36

cd /mnt/extra/CompleteGenomics/refseq
~/bin/blast/formatdb -i BUILD.36.1.REFSEQ.FA -p F -n build36 -o T -l formatdb.log

775 1349 GTGCTGGATTTGCGGGATGT TGGTCTGCATACAAAGTGCACAA none


}
