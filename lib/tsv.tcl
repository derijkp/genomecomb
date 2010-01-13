proc tsv_select {query {qfields {}} {sortfields {}}} {
	set f stdin
	set line [gets $f]
	set pos [tell $f]
	seek $f $pos
	if {[string index $line 0] eq "#"} {set line [string range $line 1 end]}
	set header [split $line \t]
	if {[llength $qfields]} {
		set qposs [list_cor $header $qfields]
		puts [join $qfields \t]
	} else {
		set qposs [list_cor $header $header]
		puts $line
	}
	list_unmerge [regexp -all -inline {[$]([a-zA-z0-9]+)} $query] 1 fields
	foreach field [list_remdup $fields] {
		set pos [lsearch $header $field]
		if {$pos == -1} {error "field \"$field\" not present"}
		incr pos
		regsub -all \\\$${field}(\[^A-Za-z\]) $query \$$pos\\1 query
	}
	set awk {BEGIN {FS="\t" ; OFS="\t"} }
	append awk $query
	set qposs [lmath_calc $qposs + 1]
	append awk " \{print $[join $qposs ,$]\}"
	if {![llength $sortfields]} {
		exec awk $awk <@ stdin >@ stdout
	} else {
		set poss [list_cor $header $sortfields]
		if {[lsearch $poss -1] != -1} {error "fields [join [list_sub $sortfields [list_find $poss -1]] ,] not found"}
		set poss [lmath_calc $poss + 1]
		set sort "gnusort8 -t \\t -V -k[join $poss " -k"]"
		eval exec [list awk $awk] | $sort <@ stdin >@ stdout
	}
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
	set command "tail +2 [list $filename] | gnusort8 -t \\t -V -k[join $poss " -k"] >@ stdout"
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
		if {[string index $line 0] eq "#"} continue
		if {![string length $line]} continue
		break
		set keep [tell $f]
		set header $line
	}
	if {[string index $line 0] eq ">"} {
		return [split [string range $line 1 end] \t]
	}
	if {![info exists header]} {set header $line}
	if {[string index $header 0] eq "#"} {
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

	set filename /data/db/_data_db_ucsc-exapted_repeats.tsv
	set fields {chrom chromStart chromEnd}

	# set f [open /complgen/compar/78vs79_compar-filter-sc.tsv]
	set f [open /complgen/compar/78vs79_compar_pvt.tsv]
	set query {compar df sample "|79 78,79" refcons "" ns "" lowscore "" trf "" str "" rp "" sd "" sc "" dbsnp "" loc EXON}

}
