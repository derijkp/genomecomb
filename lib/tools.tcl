#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc assert {check message} {
	if {![uplevel expr $check]} {
		error $message
	}
}

proc histogram {list aVar} {
	upvar $aVar a
	unset -nocomplain a
	foreach el $list {
		if {![info exists a($el)]} {
			set a($el) 1
		} else {
			incr a($el)
		}
	}
}

proc max {args} {
	if {[llength $args] == 1} {set args [lindex $args 0]}
	lmath_max [list_remove $args {} - ?]
}

proc min {args} {
	if {[llength $args] == 1} {set args [lindex $args 0]}
	lmath_min [list_remove $args {} - ?]
}

proc opensqlite3 {dbfile query} {
	set f [open "| sqlite3 -separator \"\t\" $dbfile \"$query\""]
}

proc timestamp {} {
	set milliseconds [clock milliseconds]
	set ms [format %03d [expr {$milliseconds % 1000}]]
	set seconds [expr {$milliseconds / 1000}]
	return [clock format $seconds -format "%Y-%m-%d %H:%M:%S"].$ms
}

proc chrindexseek {file f chr} {
	set root [gzroot $file]
	file mkdir $root.index
	set indexfile [indexdir_file $file chrindex ok]
	if {!$ok} {
		set tf [gzopen $file]
		set header [gets $tf]
		set chrpos [tsv_basicfields $header 1 0]
		set prevchr {}
		set list {}
		set o [open $indexfile w]
		while {![eof $tf]} {
			set pos [tell $tf]
			set line [gets $tf]
			set chr [chr_clip [lindex $line $chrpos]]
			if {$chr ne $prevchr} {
				puts $o $chr\t$pos
				set prevchr $chr
			}
		}
		gzclose $tf
		close $o
	}
	set trfchrpos [split [string trim [file_read $indexfile]] \n\t]
	set chr [chr_clip $chr]
	if {[catch {set fpos [dict get $trfchrpos $chr]}]} {
		seek $f 0 end
	} else {
		seek $f $fpos start
	}
}

proc chrgzopen {file chr chrfield headerVar} {
	upvar $headerVar header
	set f [gzopen $file]
	set header [tsv_open $f]
	set chrpos [lsearch $header $chrfield]
	if {$chrpos == -1} {
		return $f
	}
	set root [gzroot $file]
	file mkdir $root.index
	set indexfile [indexdir_file $file chrindex ok]
	if {!$ok} {
		set compressed [gziscompressed $file]
		if {$compressed} {
			gzclose $f
			set tempfile [tempfile]
			cg zcat $file > $tempfile
			set f [gzopen $tempfile]
		}
		set prevchr {}
		set list {}
		set o [open $indexfile w]
		while {![eof $f]} {
			set pos [tell $f]
			set line [gets $f]
			set chr [chr_clip [lindex $line $chrpos]]
			if {$chr ne $prevchr} {
				puts $o $chr\t$pos
				set prevchr $chr
			}
		}
		gzclose $f
		close $o
	}
	set trfchrpos [split [string trim [file_read $indexfile]] \n\t]
	set chr [chr_clip $chr]
	if {[catch {set fpos [dict get $trfchrpos $chr]}]} {
		set f [gzopen $file]
		set header [tsv_open $f]
		return $f
	} else {
		set f [gzopen $file $fpos]
		return $f
	}
}

proc ifcatch {command varName args} {
	upvar $varName result
	set error [uplevel [list catch $command $varName]]
	if {!$error} {
		return $result
	}
	set switchlist [list_pop args]
	if {[lindex $switchlist end-1] ne "default"} {
		lappend switchlist default "error [list $result]"
	}
	uplevel switch $args [list $result $switchlist]
}

if 0 {

	ifcatch {error test} result {
		test1 {puts ERRORtest1}
		test2 {puts ERRORtest2}
		default {puts ERRORdefault}
	}

	ifcatch {error test2} result {
		test1 {puts ERRORtest1}
		test2 {puts ERRORtest2}
	}

	ifcatch {error test} result {
		test1 {puts ERRORtest1}
		test2 {puts ERRORtest2}
	}

	ifcatch {set a 1} result {
		test1 {puts ERRORtest1}
		test2 {puts ERRORtest2}
	}

}

proc dict_get_default {d key {default {}}} {
	if {[catch {dict get $d $key} result]} {
		return $default
	} else {
		return $result
	}
}

# listsamples returns the list of samples, defined as the part after the last - in the fieldname
# e.g. sample2 for zyg-gatk-rdsbwa-sample2
proc listsamples {header {pattern {}}} {
	set names {}
	foreach col $header {
		set pos [string last - $col]
		if {$pos != -1} {
			incr pos
			set name [string range $col $pos end]
			if {$pattern eq "" || [string match $pattern $name]} {
				lappend names $name
			}
		}
	}
	list_remdup $names
}

# listanalyses returns the list of analyses, defined as the part after the first - in the fieldname
# e.g. gatk-rdsbwa-sample2 for zyg-gatk-rdsbwa-sample2
proc listanalyses {header {pattern {}} {analysisfieldsVar {}}} {
	if {$analysisfieldsVar ne ""} {
		upvar $analysisfieldsVar analysisfields
		set analysisfields {}
		set useanalysisfields 1
	} else {
		set useanalysisfields 0
	}
	set names {}
	foreach col $header {
		set pos [string first - $col]
		if {$pos != -1} {
			set name [string range $col [expr {$pos+1}] end]
			if {$pattern eq "" || [string match $pattern $name]} {
				lappend names $name
				if {$useanalysisfields} {
					lappend analysisfields [string range $col 0 [expr {$pos-1}]]
				}
			}
		}
	}
	list_remdup $names
}

# The samples command returns the list of "samples" as originally defined and used in the code 
# in many places (= same as analyses now)
# The samples command is kept this way for compatiblity reasons
proc samples {header {pattern {}}} {
	set names {}
	foreach col $header {
		set pos [string first - $col]
		if {$pos != -1} {
			incr pos
			set name [string range $col $pos end]
			if {$pattern eq "" || [string match $pattern $name]} {
				lappend names $name
			}
		}
	}
	list_remdup $names
}

proc trans {trans value} {
	if {[dict exists $trans $value]} {
		return [dict get $trans $value]
	} {
		return $value
	}
}

proc lforeach {args} {
	set result {}
	set pattern [list_pop args]
	set code "lappend result \"$pattern\""
	foreach {*}$args $code
	return $result
}

proc sourcename base {
	if {![regexp {^[^-]*-(.+)$} $base temp name]} {
		set name $base
	}
	return $name
}

proc trimformat args {
	string trimright [string trimright [::format {*}$args] 0] .
}

proc oargserr {cmd def} {
	set result [list $cmd]
	foreach el $def {
		set field [lindex $el 0]
		if {[llength $el] > 1} {
			lappend result ?$field?
		} else {
			lappend result $field
		}
	}
	return [join $result " "]
}

proc oargs {cmd def arg} {
	set parseopts 1
	set pos 0
	set len [llength $def]
	if {[lindex $def end] eq "args"} {
		set useargs 1
		list_pop def
		incr len -1
	} else {
		set useargs 0
	}
	foreach el $def {
		if {[llength $el] > 1} {
			set defa([lindex $el 0]) [lindex $el 1]
		}
	}
	set todo {}
	set optargs {}
	foreach el $arg {
		if {[info exists field]} {
			set a($field) $el
			if {![info exists defa($field)]} {
				if {$useargs} {
					lappend optargs -$field $el
				} else {
					error "unknown option -$field, cmd should be: [oargserr $cmd $def]"
				}
			}
			unset field
		} elseif {$parseopts && [string index $el 0] eq "-"} {
			if {$el eq "--"} {
				# no more options
				set parseopts 0
			} else {
				set field [string range $el 1 end]
			}
		} else {
			lappend todo $el
		}
	}
	if {[info exists field]} {
		error "option -$field without value, should be: [oargserr $cmd $def]"
	}
	set apos 0
	set len [llength $todo]
	foreach el $def {
		set field [lindex $el 0]
		if {[info exists a($field)]} {
			uplevel [list set $field $a($field)]
		} elseif {$apos < $len} {
			uplevel [list set $field [lindex $todo $apos]]
			incr apos
		} elseif {[info exists defa($field)]} {
			uplevel [list set $field $defa($field)]
		} else {
			error "missing arg(s): $field, should be: [oargserr $cmd $def]"
		}
	}
	set todo [lrange $todo $apos end]
	lappend todo {*}$optargs
	if {[llength $todo] && !$useargs} {
		error "too many args: [join $todo ,], should be: [oargserr $cmd $def]"
	}
	uplevel [list set args $todo]
}

proc multimatch {patterns value} {
	set final 1
	foreach temp $patterns {
		if {![string match $temp $value]} {
			set final 0; break
		}
	}
	return $final
}

proc deindent {text} {
	regsub "\n\t*\$" $text {} text
	if {[string index $text 0] eq "\n"} {set start 1} else {set start 0}
	set text [string range $text $start end]
	if {[regexp {^[\t]*} $text temp]} {
		set text [string_change [string range $text [string length $temp] end] [list \n$temp \n]]
	}
}

proc unquote {string} {
	if {[string index $string 0] eq {"} && [string index $string end] eq {"}} {
		return [string range $string 1 end-1]
	}
	return $string
}

proc catchexpr {expr {default {}}} {
	if {![catch {uplevel expr $expr} result]} {
		return $result
	} else {
		return $default
	}
}

proc catchformat {format value {default {}}} {
	if {![catch {format $format $value} result]} {
		return $result
	} else {
		return $default
	}
}

proc catchdef {cmd errorvalue} {
	if {[catch {uplevel $cmd} result]} {
		return $errorvalue
	} else {
		return $result
	}
}

proc list_replace {list replacelist} {
	foreach {old new} $replacelist {
		set pos [lsearch $list $old]
		if {$pos != -1} {
			set list [lreplace $list $pos $pos $new]
		}
	}
	return $list
}

proc putsrvars {args} {
	foreach var $args {
		if {[catch {uplevel [list set $var]} value]} {
			puts "rm($var)"
		} else {
			if {$value ne "" && [string is double $value]} {
				puts $var=$value
			} else {
				puts $var=\"$value\"
			}
		}
	}
}

proc command_exists {command} {
	set len [llength [info commands $command]]
	if {$len} {return $len}
	auto_load $command
	llength [info commands $command]
}

proc command_list {pattern} {
	list_union [array names ::auto_index $pattern] [info commands $pattern]
}

