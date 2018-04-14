#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

package require Extral

if 0 {
	set dbdir /complgen/refseq/hg19/
	set regionfile data/reg_genome_seq.tsv
}

proc cg_genome_seq {args} {
	set extraseq 124
	set freql 0
	set freqN 0.2
	set snpdbpatterns snp
	set delsize 5
	set repeats s
	set gc -1
	set gcsplit {}
	set displaygc 1
	set split 0
	set id {}
	set makemap 0
	set concatlen -1
	set econcatlen 0
	set aconcat ""
	set aconcatlen 0
	set outfile {}
	cg_options genome_seq args {
		-f - -freq {
			set freql $value
		}
		-fp - -freqp {
			if {$value != -1} {set value [expr {$value/100.0}]}
			set freql $value
		}
		-n - -freqn {
			set freqN $value
		}
		-np - -freqnp {
			if {$value != -1} {set value [expr {$value/100.0}]}
			set freqN $value
		}
		-p - -snpdbpattern {
			set snpdbpatterns $value
		}
		-d - -delsize {
			set delsize $value
		}
		-r - -repeatmasker {
			set repeats $value
		}
		-i - -id {
			set id $value
		}
		-g - -gc {
			set gc $value
		}
		-gs - -gcsplit {
			if {![isdouble $value]} {error "$value is not a number"}
			set gcsplit $value
			if {$gc == -1} {set gc 100}
		}
		-gd - -gcdisplay {
			set displaygc [true $value]
			if {$gc == -1} {set gc 0}
		}
		-s - -split {
			set split [true $value]
		}
		-c - -concat {
			set concat $value
			set concatlen [string length $concat]
		}
		-e - -ce - -concatend {
			set econcat $value
			set econcatlen [string length $econcat]
		}
		-ca - -concatadj {
			set aconcat $value
			set aconcatlen [string length $econcat]
		}
		-cn - -concatname {
			set concatname $value
		}
		-m - -mapfile {
			set mapfile $value
			set makemap 1
		}
		-l - -limitchars {
			set limitchars $value
		}
		-namefield {
			set namefield $value
		}
	}  {regionfile dbdir outfile} 2 3
	set dbdir [dbdir $dbdir]
	if {$outfile ne ""} {
		set root [file root $outfile]
		set ext [file extension $outfile]
		set tail [file root [file tail $outfile]]
	}
	if {[isdouble $gcsplit]} {
		if {$outfile eq ""} {error "In order to use the gcsplit option, you have to provide outfile"}
		if {$concatlen != -1} {error "Cannot combine concat and gcsplit options"}
	}
	if {$split} {
		if {$outfile eq ""} {error "In order to use the split option, you have to provide outfile"}
		if {$concatlen != -1} {error "Cannot combine concat and split options"}
		if {[isdouble $gcsplit]} {
			file mkdir $root-lowgc
			file mkdir $root-highgc
		}
	}
	#
	catch {close $f}; catch {close $fg}; catch {close $fo}; catch {close $foh}
	if {$outfile eq ""} {
		set fo stdout
		set outf $fo
	} elseif {$split} {
		set filenum 1
	} else {
		if {![isdouble $gcsplit]} {
			set fo [open $outfile w]
		} else {
			set fo [open $root-lowgc$ext w]
			set foh [open $root-highgc$ext w]
		}
		set outf $fo
	}
	set fg [genome_open [lindex [glob $dbdir/genome_*.ifas] 0]]
	if {![file exists $regionfile]} {
		set regionlist [list_remove [split $regionfile ":-, \n"] {}]
		set regionfile [tempfile]
		set f [open $regionfile w]
		puts $f [join {chromosome begin end} \t]
		foreach {chr begin end} $regionlist {
			if {![isint $begin] || ![isint $end]} {
				error "file $regionfile does not exist, and the argument is also not a properly formatted regionlist"
			}
			puts $f $chr\t$begin\t$end
		}
		close $f
	}
	set f [gzopen $regionfile]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 3]
	if {$id ne ""} {
		set idpos [lsearch $header $id]
		if {$idpos == -1} {
			error "id Column $id not found in header"
		}
	} else {
		set idpos -1
	}
	if {[info exists namefield]} {
		set namepos [lsearch $header $namefield]
		if {$namepos == -1} {
			error "namefield $namefield not found"
		}
	} else {
		set namepos [lsearch $header name]
	}
	if {$concatlen >= 0} {
		if {[info exists concatname]} {set name $concatname} else {set name [file root [file tail $regionfile]]}
		puts $fo "\>$name concatenated"
		set firstline 1
	} else {
		set name concat
	}
	if {$makemap} {
		set fm [open $mapfile w]
		puts $fm [join {chromosome begin end destchromosome destbegin destend name} \t]
	}
	set fstart 0
	if {$concatlen >= 0 && $econcatlen} {
		puts -nonewline $fo $econcat
		incr fstart $econcatlen
	}
	set pchr {}
	set pend {}
	set nextrecord ""
	set nextrecordh ""
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		set sub [list_sub $line $poss]
		putslog $sub
		foreach {chr estart eend} $sub break
		set chr [chr_clip $chr]
		set seq [genome_get $fg $chr [expr {$estart}] [expr {$eend}]]
		set mseq [genome_mask $dbdir $seq $chr [expr {$estart}] [expr {$eend}] $freql $freqN $delsize $repeats $snpdbpatterns]
		if {$concatlen == -1} {
			set name [join [list_sub $sub {0 1 2}] -]
			if {$idpos != -1} {
				set name "[lindex $line $idpos] $name"
			}
			if {$gc == 0} {
				set gcval [seq_gc $seq]
				if {$displaygc} {append name " GC:[format %.1f $gcval]"}
			} elseif {$gc != -1} {
				set gcval [lmath_max [seq_gc $seq $gc]]
				if {$displaygc} {append name " GC:[format %.1f [seq_gc $seq]] maxGC($gc):[format %.1f $gcval]"}
			}
			if {[info exists limitchars]} {
				regsub -all {[^A-Za-z0-9_.-]} $name $limitchars name
			}
			if {$split} {
				if {[isdouble $gcsplit]} {
					if {$gcval >= $gcsplit} {
						set outf [open $root-highgc/[lindex $name 0].fas w]
					} else {
						set outf [open $root-lowgc/[lindex $name 0].fas w]
					}
				} else {
					set outf [open $root[lindex $name 0].fas w]
				}
				incr filenum
				puts $outf \>$name
			} elseif {[isdouble $gcsplit] && $gcval >= $gcsplit} {
				set outf $foh
				puts $outf $nextrecordh\>$name
				set nextrecordh \n
			} else {
				set outf $fo
				puts $outf $nextrecord\>$name
				set nextrecord \n
			}
		} elseif {!$firstline} {
			foreach {chr begin end} $sub break
			if {[chr_compare $chr $pchr] == 0 && $begin == $pend} {
				puts -nonewline $outf $aconcat
				incr fstart $aconcatlen
			} else {
				puts -nonewline $outf $concat
				incr fstart $concatlen
			}
			set pchr $chr
			set pend $end
		}
		if {$split} {
			puts $outf $mseq
			close $outf
		} else {
			puts -nonewline $outf $mseq
		}
		if {$makemap} {
			puts $fm $name\t$fstart\t[expr {$fstart+[string length $mseq]}]\t[join $sub \t]\t[lindex $line $namepos]
		}
		incr fstart [string length $mseq]
		set firstline 0
	}
	if {$concatlen >= 0 && $econcatlen} {
		puts -nonewline $fo $econcat
		incr fstart $econcatlen
	}
	if {!$split} {
		puts $fo ""
		if {[isdouble $gcsplit]} {puts $foh ""}
	}
	if {$makemap} {
		close $fm
	}
	gzclose $f; genome_close $fg
}
