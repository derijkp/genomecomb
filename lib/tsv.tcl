proc tsv_select {query {qfields {}} {sortfields {}}} {
	set f stdin
	set header [tsv_open $f]
	set keep [tell $f]
	set awk ""
	set sort ""
	set cut ""
	if {[llength $sortfields]} {
		set poss [list_cor $qfields $sortfields]
		if {[lsearch $poss -1] != -1} {
			set poss [list_cor $header $sortfields]
			if {[lsearch $poss -1] != -1} {error "fields [join [list_sub $sortfields [list_find $poss -1]] ,] not found"}
			if {$qfields ne ""} {
				set qposs [list_cor $header $qfields]
				set qposs [lmath_calc $qposs + 1]
				set cut "cut -d \\t -f [join $qposs ,]"
			}
		}
		set poss [lmath_calc $poss + 1]
		set keys {}
		foreach pos $poss {
			lappend keys $pos,$pos
		}
		set sort "gnusort8 -t \\t -V -s -k[join $keys " -k"]"
	}
	if {$query ne ""} {
		list_unmerge [regexp -all -inline {[$]([a-zA-z0-9]+)} $query] 1 fields
		foreach field [list_remdup $fields] {
			set pos [lsearch $header $field]
			if {$pos == -1} {error "field \"$field\" not present"}
			incr pos
			regsub -all \\\$${field}(\[^A-Za-z\]) $query \$$pos\\1 query
		}
		set awk {BEGIN {FS="\t" ; OFS="\t"} }
		append awk $query
		if {($qfields ne "") && ($cut eq "")} {
			set qposs [list_cor $header $qfields]
		} else {
			set qposs [list_cor $header $header]
		}
		set qposs [lmath_calc $qposs + 1]
		append awk " \{print $[join $qposs ,$]\}"
	} elseif {($qfields ne "") && ($cut eq "")} {
		set qposs [list_cor $header $qfields]
		set qposs [lmath_calc $qposs + 1]
		set cut "cut -d \\t -f [join $qposs ,]"
	}
	set pipe {}
	if {$awk ne ""} {
		lappend pipe [list awk $awk]
	}
	if {$sort ne ""} {
		lappend pipe $sort
	}
	if {$cut ne ""} {
		lappend pipe $cut
	}
	puts stderr pipe:[join $pipe " | "]
	if {$qfields ne ""} {puts [join $qfields \t]} else {puts [join $header \t]}
	seek $f $keep
	if {![llength $pipe]} {
		fcopy stdin stdout
	} else {
		eval exec [join $pipe " | "] <@ stdin >@ stdout
	}
}

if 0 {
	lappend auto_path ~/dev/completegenomics/lib
	package require Tclx
	signal -restart error SIGINT
	package require Extral
	set f [open GS102/ASM/var-GS000000078-ASM.tsv]
	set query "\$begin < 2000"
	set qfields "chromosome begin end"
	set sortfields "haplotype"
	set sortfields "chromosome begin"

cg select -q '' < GS102/ASM/var-GS000000078-ASM.tsv | less
cg select -f 'haplotype chromosome begin' < GS102/ASM/var-GS000000078-ASM.tsv | less
cg select -q '$begin < 2000' < GS102/ASM/var-GS000000078-ASM.tsv | less

cg select -q '$begin < 2000' < GS102/ASM/var-GS000000078-ASM.tsv > /tmp/test1
cg select -q '$begin < 2000' -f 'chromosome begin end' -s 'chromosome begin' < GS102/ASM/var-GS000000078-ASM.tsv > /tmp/test2
cg select -q '$begin < 2000' -f 'chromosome begin end' -s 'haplotype' < GS102/ASM/var-GS000000078-ASM.tsv > /tmp/test3

}

proc tsv_sort {filename fields} {
	set f [open $filename]
	set line [gets $f]
	if {[string index $line 0] eq "#"} {set line [string range $line 1 end]}
	set header [split $line \t]
	set poss [list_cor $header $fields]
	if {[lsearch $poss -1] != -1} {error "fields [join [list_sub $fields [list_find $poss -1]] ,] not found"}
	set poss [lmath_calc $poss + 1]
	puts [join $header \t]
	set command "tail +2 [list $filename] | gnusort8 -t \\t -V -s -k[join $poss " -k"] >@ stdout"
	eval exec $command
}

proc file_rootgz {filename} {
	if {[file extension $filename] eq ".gz"} {
		return [file root [file root $filename]]
	} else {
		file root [file root $filename]
	}
}

proc tsv_open {f} {
	set keep 0
	while {![eof $f]} {
		set line [gets $f]
		if {![string length $line]} continue
		if {[string index $line 0] ne "#"} break
		set keep [tell $f]
		set header $line
	}
	if {[string index $line 0] eq ">"} {
		return [split [string range $line 1 end] \t]
	}
	if {![info exists header]} {
		return [split $line \t]
	} elseif {[string index $header 0] eq "#"} {
		seek $f $keep
		return [split [string range $header 1 end] \t]
	} else {
		seek $f $keep
		return [split $header \t]
	}
}

if 0 {
	lappend auto_path ~/dev/completegenomics/lib
	package require Tclx
	signal -restart error SIGINT
	package require Extral
	cd /complgen/compar

	set filename /data/db/ucsc_ori/_data_db_ucsc-exapted_repeats.tsv
	set f [open $filename]
	set fields {chrom chromStart chromEnd}

	# set f [open /complgen/compar/78vs79_compar-filter-sc.tsv]
	set f [open /complgen/compar/78vs79_compar_pvt.tsv]
	set query {compar df sample "|79 78,79" refcons "" ns "" lowscore "" trf "" str "" rp "" sd "" sc "" dbsnp "" loc EXON}

}

proc tsv_next {f xpos next {shift 100000}} {
	# do a ~ binary search to get at next faster
	set start [tell $f]
	while 1 {
		seek $f $shift current
		gets $f
		set line [split [gets $f] \t]
		set x [lindex $line $xpos]
		if {![isdouble $x] || ($x >= $next)} {
			seek $f $start
			set shift [expr {$shift / 2}]
			if {$shift < 1000} break
		}
		set start [tell $f]
	}
	while {![eof $f]} {
		set fpos [tell $f]
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set x [lindex $line $xpos]
		if {$x >= $next} break
	}
	if {![eof $f]} {
		return $fpos
	} else {
		return $x
	}
}

proc tsv_index {file xfield} {
	if {[inlist {.rz} [file extension $file]]} {
		set indexname [file root $file].${xfield}_index
		set tempfile [tempfile]
		exec razip -d -c $file > $tempfile
		set f [open $tempfile]
	} else {
		set indexname $file.${xfield}_index
		set f [open $file]
	}
	set header [tsv_open $f]
	set xpos [lsearch $header $xfield]
	if {$xpos == -1} {error "field $xfield not present in file $file"}
	set fstart [tell $f]
	set fpos $fstart
	set line [split [gets $f] \t]
	set xmin [lindex $line $xpos]
	set findex [expr {$xmin-$xmin%10000}]
	set prev $findex
	set next [expr {$prev + 10000}]
	set index [list $fpos]
	catch {Classy::Progress start [file size $file] "Making index"}
	while {![eof $f]} {
		set fpos [tsv_next $f $xpos $next]
		if {[eof $f]} {
			set xmax $fpos
			break
		}
		lappend index $fpos
		incr prev 10000
		incr next 10000
		catch {Classy::Progress set [tell $f]}
		if {![expr $next%100000]} {puts stderr $next}
	}
	catch {Classy::Progress stop}
	close $f
	set size [file size $file]
	set f [open $file]
	seek $f [expr {$size-5000}] start
	gets $f
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set xmax [lindex $line $xpos]
	}
	set o [open $indexname w]
	puts $o 10000
	puts $o $findex
	puts $o $xmin
	puts $o $xmax
	puts $o [join $index \n]
	close $o
}
