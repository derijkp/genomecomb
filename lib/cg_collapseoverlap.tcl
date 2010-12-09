#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

package require Extral


if 0 {

	package require Tclx
	signal -restart error SIGINT
	lappend auto_path /home/peter/dev/completegenomics/lib
	package require Extral
	cd /complgen/refseq/hg18test
	set file /complgen/refseq/hg18test/ucsc_hg18_oreganno.tsv
	set resultfile /complgen/refseq/hg18test/reg_hg18_oreganno.tsv
	set file /complgen/refseq/hg18test/ucsc_hg18_rmsk.tsv
	set resultfile /complgen/refseq/hg18test/reg_hg18_rmsk.tsv
	set file /complgen/refseq/hg18test/ucsc_hg18_snp130.tsv
	set resultfile /complgen/refseq/hg18test/reg_hg18_snp130.tsv

cg collapseoverlap /complgen/refseq/hg18test/ucsc_hg18_oreganno.tsv

}

proc collapseoverlap_join {cur scorepos} {
	if {[llength $cur] == 1} {return [lindex $cur 0]}
	if {$scorepos != -1} {
		set cur [lsort -dict -decreasing -index $scorepos $cur]
	}
	set result {}
	set len [llength [lindex $cur 0]]
	for {set i 0} {$i < $len} {incr i} {
		set temp [list_subindex $cur $i]
		set nodup [list_remdup $temp]
		if {[llength $nodup] < 2} {
			lappend result $nodup
		} else {
			lappend result [join $temp ,]
		}
	}
	return $result
}

proc collapseoverlap {file resultfile} {
	puts "making $resultfile"

	catch {close $f} ; catch {close $o}
	if {[catch {open $file} f]} {
		error "Could not open file $file"
	}
	set cor [open_region $f header]
	foreach {chrpos startpos endpos} $cor break
	set scorepos [lsearch $header score]
	if {[catch {open $resultfile.temp w} o]} {
		error "Could not write outputfile $resultfile"
	}
	puts $o [join $header \t]
	set line [split [gets $f] "\t"]
	foreach {chr start end} [list_sub $line $cor] break
	set cur [list $line]
	set num 0; set next 100000
	while {![eof $f]} {
		incr num; if {$num >= $next} {puts $num; incr next 100000}
		set line [split [gets $f] "\t"]
		foreach {cchr cstart cend} [list_sub $line $cor] break
		set newchr [expr {$cchr ne $chr}]
		if {$newchr || ($cstart > $start)} {
			# write preceeding
			set stop 0
			set keepend -1
			foreach l $cur {
				set end [lindex $l $endpos]
				if {!$newchr && ($end > $cstart)} {
					set end $cstart
					set stop 1
				}
				if {$end != $keepend} {
					set joined [collapseoverlap_join $cur $scorepos]
					lset joined $startpos $start
					lset joined $endpos $end
					puts $o [join $joined \t]
					set start $end
				}
				if {$stop} break
				set keepend $end
				list_shift cur
			}
		}
		# add new line
		lappend cur $line
		if {[llength $cur] > 1} {
			set cur [lsort -real -index $endpos $cur]
		}
		set chr $cchr
		set start $cstart
	}
	file rename $resultfile.temp $resultfile
	close $f
	close $o
	puts "Finished $resultfile"
}

proc cg_collapseoverlap {args} {
	if {([llength $args] < 1)} {
		puts stderr "format is: $::base file ..."
		puts stderr " - Collapses overlapping regions in a region file."
		puts stderr " - makes a new file with reg_ prepended to the original filename"
		puts stderr " - Removal of overlap can be done by taking only the highest"
		puts stderr " - scoring region (this is always done when score is available)"
		puts stderr " - or taking all regions in 1 line (if score is not available)"
		exit 1
	}
	foreach {path} $args break
	puts "----------------------------------------------------"
	foreach file $args {
		set path [file dir $file]
		set tail [file tail $file]
		if {[string range $tail 0 4] eq "ucsc_"} {set tail [string range $tail 5 end]}
		set resultfile ${path}/reg_$tail
		if {[file exists $resultfile]} {
			puts "Skipping $resultfile: already exists"
			continue
		}
		collapseoverlap $file $resultfile
	}
}

if {[info exists argv]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base $scriptname
	cg_downloaddb {*}$argv
}


