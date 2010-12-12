if 0 {
	set path /complgen/refseq/hg18
	set build hg18
	set dbname snp130
}

package require BioTcl

proc downloaddb_dbsnp_convline {line} {
	foreach {chr begin end class ref observed name} $line break
	set observed [split $observed /]
	set rp [lsearch $observed $ref]
	if {$rp == -1} {
		set complement [seq_complement $observed]
		set rp [lsearch $complement $ref]
		if {$rp != -1} {
			set observed $complement
		}
	}
	set observed [list_remove $observed $ref]
	if {$ref eq "-"} {set ref ""}
	set observed [list_change $observed {- {}}]
	set rlen [string length $ref]
	set alen {}
	foreach el $observed {
		lappend alen [string length $el]
	}
	if {$rlen == 0} {
		set type ins
	} elseif {[inlist $alen 0]} {
		set type del
	} elseif {$rlen == 1 && [inlist $alen 1]} {
		set type snp
	} else {
		set type sub
	}
	set len [llength $observed]
	set name [list_fill $len $name]
	list $chr $begin $end $type $ref [join $observed ,] [join $name ,]
}

proc downloaddb_dbsnp {path build dbname} {
	set ufilename $path/ucsc_${build}_$dbname.tsv
	set filename $path/var_${build}_$dbname.tsv
	puts "Making $filename"
	if {[file exists $filename]} {
		puts "file $filename exists: skipping"
		return
	}
	if {![file exists $ufilename]} {
		downloaddb $path $build $dbname
	}
	puts "Converting $ufilename"
	catch {close $f} ; catch {close $o}
	set f [open $ufilename]
	set o [open $filename.temp w]
	set header [split [gets $f] \t]
	set poss [list_cor $header {chrom start end class refNCBI observed name}]
	puts $o [join {chrom start end type ref alt name} \t]
	set pline [list_sub [split [gets $f] \t] $poss]
	set pline [downloaddb_dbsnp_convline $pline]
	set num 0 ; set next 100000
	while {![eof $f]} {
		incr num; if {$num >= $next} {puts $num; incr next 100000}
		set line [split [gets $f] \t]
		set line [list_sub $line $poss]
		set line [downloaddb_dbsnp_convline $line]
		if {[lrange $pline 0 3] eq [lrange $line 0 3]} {
			foreach {alt name} [list_sub $pline {5 6}] break
			append alt ,[lindex $line 5]
			append name ,[lindex $line 6]
			lset pline 5 $alt
			lset pline 6 $name
		} else {
			puts $o [join $pline \t]
			set pline $line
		}
	}
	close $f ; close $o
	puts "Sorting $filename"
	cg select -f {chrom start end type} $filename.temp $filename.temp2
	file rename $filename.temp2 $filename
	file delete $filename.temp
}
