#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc varfileversion {} {
	return 0.10.0
}

proc fileversion {} {
	return 0.10.0
}

proc fileinfo {file args} {
	set result {}
	if {[file exists $file.info.tsv]} {
		set result [dict merge $result [infofile_read $file.info.tsv]]
	}
	if {[file isdir $file]} {
		set srcdir [file dir $file]
	}
	if {[file exists $srcdir/sampleinfo.tsv]} {
		set result [dict merge $result [infofile_read $srcdir/sampleinfo.tsv]]
	}
	if {[file exists $srcdir/info.txt]} {
		set c [split [file_read $srcdir/info.txt] \n]
		foreach line $c {
			foreach {key value} [split [string range $line 1 end] \t] break
			dict set result $key $value
		}
	}
	if {![file isdir $file]} {
		set f [gzopen $file]
		set header [tsv_open $f comment]
		set result [dict merge $result [comment2dict $comment]]
		dict set $result header $header
		if {![dict exists $result ref]} {
			if {[dict exists $result dbdir]} {
				dict set result ref [dict get $result dbdir]
			} elseif {[dict exists $result GENOME_REFERENCE]} {
				set ref [lindex [dict get $result GENOME_REFERENCE] end]
				if {$ref < 38} {set ref [expr {$ref - 18}]}
				dict set result ref hg$ref
			}
		}
	}
	return $result
}

# potential values in 
# key	examplevalue	meaning
# genomecomb	0.10	genomecomb version
# filetype	var	type of tsv file (e.g. var for variant list)
# split	1/0	variants with different alternative alleles are on separate lines (split = 1)
# ref	hg19	reference database/genome
# field	name 
# header "chromosome begin end ..."	fields in tsv file

# todo: convert vcf data
proc comment2dict {comment} {
	set result [dict create]
	set order {}
	foreach line [split [string trim $comment] \n] {
		if {[string index $line 0] ne "#"} continue
		set line [string range $line 1 end]
		set pos [string first \t $line]
		if {$pos != -1} {
			set key [string range $line 0 [expr {$pos-1}]]
			set value [string trim [string range $line [expr {$pos+1}] end]]
		} else {
			set pos [string first = $line]
			if {$pos != -1} {
				set key [string range $line 0 [expr {$pos-1}]]
				set value [string range $line [expr {$pos+1}] end]
			} else {
				set key {}
				set value $line
			}
		}
		set key [string trim $key]
		lappend order $key
		dict lappend result $key $value
	}
	dict set result order [list_remdup $order]
	return $result
}

proc dict2comment {infodict} {
	if {[dict exists $infodict order]} {
		set order [dict get $infodict order]
		set order [list_common {filetype fileversion} $order]
	} else {
		set order {filetype fileversion}
	}
	set keys [dict keys $infodict]
	set keys [list_remove $keys order]
	set keys [list_common [list_union $order $keys] $keys]
	set result {}
	foreach key $keys {
		foreach value [dict get $infodict $key] {
			append result \#$key\t$value\n
		}
	}
	
	return $result
}
