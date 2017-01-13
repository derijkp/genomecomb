#!/bin/sh
# the next line restarts using wish \
exec tclsh "$0" ${1+"$@"}

#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc annotatereg {file dbfile name annotfile near dbinfo} {
# putslog [list annotatereg $file $dbfile $name $annotfile $near $dbinfo]
	set newh [dict get $dbinfo newh]
	set dataposs [dict get $dbinfo dataposs]
	catch {gzclose $f}
	set f [gzopen $file]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 3]
	gzclose $f
	set fields [list_sub $header $poss]
	if {[inlist $poss -1]} {
		error "Cannot annotate $file: wrong fields"
	}
	set f [gzopen $dbfile]
	set dbposs [open_region $f dbheader]
	gzclose $f
	set o [open $annotfile.temp w]
	puts $o [join $newh \t]
	close $o
	if {[gziscompressed $file]} {
		set file "|[gzcat $file] '$file'"
	}
	# putslog [list {*}[gzcat $dbfile] $dbfile | reg_annot $file {*}$poss - {*}$dbposs $near {*}$dataposs]
	if {[gziscompressed $dbfile]} {
		gzcatch {
			exec {*}[gzcat $dbfile] $dbfile | reg_annot $file {*}$poss - {*}$dbposs $near {*}$dataposs >> $annotfile.temp 2>@ stderr
		}
	} else {
		exec reg_annot $file {*}$poss $dbfile {*}$dbposs $near {*}$dataposs >> $annotfile.temp 2>@ stderr
	}
	file rename -force $annotfile.temp $annotfile

}

proc annotatevar {file dbfile name annotfile dbinfo} {
# putsvars file dbfile name annotfile dbinfo
	set newh [dict get $dbinfo newh]
	set dataposs [dict get $dbinfo dataposs]
	catch {gzclose $f}
	set f [gzopen $file]
	set header [tsv_open $f]
	set poss [tsv_basicfields $header 3]
	gzclose $f
	set fields [list_sub $header $poss]
	if {[inlist $poss -1]} {
		error "Cannot annotate $file: wrong fields"
	}
	set type1pos [lsearch $header type]
	if {$type1pos == -1} {
		error "$file has no type field"
	}
	set alt1pos [lsearch $header alt]
	if {$alt1pos == -1} {
		error "$file has no alt field"
	}
	set f [gzopen $dbfile]
	set header [tsv_open $f]
	gzclose $f
	set tempdbposs [tsv_basicfields $header 6 0]
	set dbposs [lrange $tempdbposs 0 2]
	set type2pos [lindex $tempdbposs 3]
	if {$type2pos == -1} {
		error "$dbfile has no type field"
	}
	set alt2pos [lindex $tempdbposs 5]
	if {$alt2pos == -1} {
		error "$dbfile has no alt field"
	}
	set o [open $annotfile.temp w]
	puts $o [join $newh \t]
	close $o
	if {[gziscompressed $file]} {
		set file "|[gzcat $file] '$file'"
	}
	set notfound -
	if {[gziscompressed $dbfile]} {
		gzcatch { 
			exec {*}[gzcat $dbfile] $dbfile | var_annot $file {*}$poss $type1pos $alt1pos - {*}$dbposs $type2pos $alt2pos $notfound {*}$dataposs >> $annotfile.temp 2>@ stderr
		}
	} else {
		# puts [list ../bin/var_annot $file {*}$poss $type1pos $alt1pos $dbfile {*}$dbposs $type2pos $alt2pos $notfound {*}$dataposs]
		exec var_annot $file {*}$poss $type1pos $alt1pos $dbfile {*}$dbposs $type2pos $alt2pos $notfound {*}$dataposs >> $annotfile.temp 2>@ stderr
	}
	file rename -force $annotfile.temp $annotfile
}

proc cg_annotatedb_info {dbfile {near -1}} {
	if {[file exists [gzroot $dbfile].opt]} {
		set a [dict create {*}[file_read [gzroot $dbfile].opt]]
	} else {
		set a [dict create]
	}
	if {[dict exists $a name]} {
		set name [dict get $a name]
	} else {
		set split [split [lindex [split [file root [file tail [gzroot $dbfile]]] -] 0] _]
		set name [lindex $split end]
	}
	dict set a name $name
	set dbtype [lindex [split [file tail $dbfile] _] 0]
	dict set a dbtype $dbtype
	set f [gzopen $dbfile]
	set header [tsv_open $f comment]
	gzclose $f
	dict set a header $header
	if {$dbtype eq "gene"} {
		set outfields {}
		set poss {}
	} elseif {$dbtype eq "var"} {
		set outfields [dict_get_default $a fields {name freq score}]
		set poss [tsv_basicfields $header 3]
	} elseif {$dbtype eq "bcol"} {
		set outfields {}
		set poss {}
	} else {
		if {[dict exists $a fields]} {
			set outfields [dict get $a fields]
		} else {
			switch -glob $name {
				1000g* {set outfields freq}
				rmsk {set outfields name}
				evofold {set outfields score}
				rnaGene {set outfields name}
				simpleRepeat {set outfields score}
				tRNAs {set outfields name}
				default {set outfields {name freq score annotation}}
			}
		}
		set poss [tsv_basicfields $header 3]
	}
	set outfields [list_common [list_remdup $outfields] $header]
	set dataposs [list_cor $header $outfields]
	if {[inlist $dataposs -1]} {
		set poss [list_find $dataposs -1]
		error "required output field(s) [join [list_sub $outfields $poss] ,] not found in file $dbfile"
	}
	if {[dict exists $a headerfields]} {
		set headerfields [dict get $a headerfields]
		set olen [llength $outfields]
		set hlen [llength $headerfields]
		if {$olen != $hlen && !($olen == 0 && $hlen == 1)} {
			error "error: number of headerfields in opt file differs from number of outfields"
		}
		set newh $headerfields
		if {$near != -1} {
			lappend newh ${name}_dist
		}
	} else {
		set newh $name
		if {$dbtype eq "gene"} {
			set newh [list ${name}_impact ${name}_gene ${name}_descr]
		} elseif {$dbtype eq "bcol"} {
			set newh [list $name]
		} else {
			switch [llength $outfields] {
			0 {
				if {$near != -1} {
					set newh ${name}_dist
				} else {
					set newh $name
				}
			}
			1 {
				set newh $name
				if {$near != -1} {
					lappend newh ${name}_dist
				}
			}
			default {
				set newh {}
				foreach key $outfields {
					lappend newh ${name}_$key
				}
				if {$near != -1} {
					lappend newh ${name}_dist
				}
			}
			}
		}
	}
	dict set a outfields $outfields
	dict set a newh $newh
	dict set a dataposs $dataposs
	dict set a basicfields $poss
	if {![dict exists $a fields]} {
		dict set a fields [dict get $a name]
	}
	return $a
}

proc cg_annotate_job {args} {
	set near -1
	set dbdir {}
	set replace e
	set multidb 0
	set upstreamsize 2000
	cg_options annotate args {
		-near {
			set near $value
		}
		-dbdir {
			set dbdir $value
		}
		-name {
			set namefield $value
		}
		-replace {
			if {$value ni "y n e"} {error "invalid value $value for -replace, must be one of y n e"}
			set replace $value
		}
		-multidb {
			set multidb $value
		}
		-u - --upstreamsize {
			set upstreamsize $value
		}
	} {orifile resultfile} 3
	set dbfiles {}
	foreach testfile $args {
		if {[file isdir $testfile]} {
			lappend dbfiles {*}[glob -nocomplain $testfile/var_*.tsv $testfile/gene_*.tsv $testfile/mir_*.tsv $testfile/reg_*.tsv $testfile/bcol_*.tsv]
		} elseif {![file exists $testfile]} {
			error "File $testfile does not exist"
		} else {
			lappend dbfiles $testfile
		}
	}
	set names {}
	set newh {}
	foreach dbfile $dbfiles {
		set dbinfo [cg_annotatedb_info $dbfile $near]
		lappend names [dict get $dbinfo name]
		lappend newh {*}[dict get $dbinfo newh]
	}
	set ext [file extension [gzroot $orifile]]
	if {[file exists $orifile]} {
		# only check for existing columns if file already exists
		if {$ext eq ".vcf"} {
			set f [gzopen $orifile]
			set header [vcf2sft_header $f]
			catch {gzclose $f}
		} else {
			set header [header $orifile]
		}
		set poss [tsv_basicfields $header 6 0]
		set common [list_common $header $newh]
		if {[llength $common]} {
			if {$replace eq "e"} {
				error "Error: field(s) [join $common ,] already in file"
			}
			if {$replace eq "n"} {
				foreach name $common {
					set skip($name) 1
				}
			}
		}
	}
	set tempbasefile [indexdir_file $resultfile vars.tsv ok]
	job_logdir $tempbasefile.log_jobs
	# If $orifile is a vcf file, convert
	set ext [file extension [gzroot $orifile]]
	if {$ext eq ".vcf"} {
		set convertedfile [indexdir_file $resultfile cvars.tsv ok]
		job annot-vcf2tsv -deps {$orifile} -targets {$convertedfile} -code {
			cg vcf2tsv -split 1 $dep $target
		}
		set orifile $convertedfile
		set usefile [indexdir_file $resultfile vars.tsv ok]
		set ok 0
	} else {
		set usefile [indexdir_file $orifile vars.tsv ok]
	}
	if {!$ok} {
		job annot-createusefile -deps {$orifile} -targets {$usefile} -vars {ok orifile usefile dbfiles} -code {
			# usefile: smaller file with only variants used for actual annotation; 
			# if orifile is small, a link to it is made.
			# If it contains to many extra columns a cut down version is made
			set f [gzopen $orifile]
			set header [tsv_open $f]
			catch {gzclose $f}
			if {[gziscompressed $orifile] || [file dir $target] ne "$orifile.index" || ([llength $header] > 10 && [llength $dbfiles] >= 4)} {
				tsv_varsfile $orifile $usefile
				puts "Using varfile $usefile"
			} else {
				if {[file dir $target] eq "$orifile.index"} {
					mklink $orifile $target
				}
			}
		}
	}
	set afiles {}
	foreach dbfile $dbfiles {
		set dbinfo [cg_annotatedb_info $dbfile $near]
		set name [dict get $dbinfo name]
		if {[info exists namefield]} {set name $namefield}
		if {[info exists skip($name)]} {
			puts "Skipping $dbfile: $name already in file"
			continue
		}
		set dbtype [lindex [split [file tail $dbfile] _] 0]
		set target $tempbasefile.${name}_annot
		lappend afiles $target
		if {$dbtype eq "gene"} {
			if {$near != -1} {error "-near option does not work with gene dbfiles"}
			if {$dbdir eq ""} {
				set dbdir [file dir [file_absolute $dbfile]]
			}
			set genomefile [lindex [glob -nocomplain $dbdir/genome_*.ifas] 0]
			if {![file exists $genomefile]} {
				error "no genomefile (genome_*.ifas) found in $dbdir, try using the -dbdir option"
			}
			job annot-[file tail $dbfile] -deps {$usefile $genomefile $dbfile} -targets {$target} -vars {genomefile dbfile name dbinfo upstreamsize} -code {
				set genecol [dict_get_default $dbinfo genecol {}]
				set transcriptcol [dict_get_default $dbinfo transcriptcol {}]
				annotategene $dep $genomefile $dbfile $name $target $genecol $transcriptcol $upstreamsize
			}
		} elseif {$dbtype eq "mir"} {
			if {$near != -1} {error "-near option does not work with gene dbfiles"}
			if {$dbdir eq ""} {
				set dbdir [file dir [file_absolute $dbfile]]
			}
			job annot-[file tail $dbfile] -deps {$usefile $dbfile} -targets {$target} -vars {dbfile name dbinfo upstreamsize} -code {
				set genecol [dict_get_default $dbinfo genecol name]
				set transcriptcol [dict_get_default $dbinfo transcriptcol transcript]
				set extracols [dict_get_default $dbinfo extracols status]
				# we do not need the genomefile for plain annotation
				set genomefile {}
				annotatemir $dep $genomefile $dbfile $name $target $genecol $transcriptcol 100 1 0 $extracols $upstreamsize
			}
		} elseif {$dbtype eq "var"} {
			if {$near != -1} {error "-near option does not work with var dbfiles"}
			job annot-[file tail $dbfile] -deps {$usefile $dbfile} -targets {$target} -vars {dbfile name dbinfo upstreamsize} -code {
				set f [gzopen $dep]
				set header [tsv_open $f]
				catch {gzclose $f}
				if {[file extension $dbfile] ne ".bcol"} {
					set altpos [lsearch $header alt]
					if {$altpos == -1} {
						puts "Skipping: $orifile has no alt field"
						file_write $target ""
					}
				}
				set outfields [dict get $dbinfo outfields]
				if {[file extension $dbfile] eq ".bcol"} {
					annotatebcolvar $dep $dbfile $name $target
				} else {
					annotatevar $dep $dbfile $name $target $dbinfo
				}
			}
		} elseif {$dbtype eq "bcol"} {
			if {$near != -1} {error "-near option does not work with bcol dbfiles"}
			job annot-[file tail $dbfile] -deps {$usefile $dbfile} -targets {$target} -vars {dbfile name} -code {
				annotatebcol $dep $dbfile $name $target
			}
		} else {
			job annot-[file tail $dbfile] -deps {$usefile $dbfile} -targets {$target} -vars {dbfile name dbinfo near} -code {
				set outfields [dict get $dbinfo outfields]
				annotatereg $dep $dbfile $name $target.temp $near $dbinfo
				file rename -force $target.temp $target
			}
		}
	}
putsvars orifile header afiles
	job annot-paste -deps [list $orifile {*}$afiles] -targets {$resultfile} -vars {orifile afiles multidb replace newh resultfile} -code {
		if {$multidb} {
			set temp2 [filetemp $resultfile]
			cg select -f id $orifile $temp2
			set temp [filetemp $resultfile]
			exec tsv_paste $temp2 {*}$afiles > $temp
			file delete $temp2
			file rename -force $temp $resultfile
		} elseif {$replace eq "y"} {
			set temp2 [filetemp $resultfile]
			set f [gzopen $orifile]
			set header [tsv_open $f]
			gzclose $f
			cg select -f [list_lremove $header $newh] $orifile $temp2
			set temp [filetemp $resultfile]
			exec tsv_paste $temp2 {*}$afiles > $temp
			file delete $temp2
			file rename -force $temp $resultfile
		} else {
			set temp [filetemp $resultfile]
			exec tsv_paste $orifile {*}$afiles > $temp
			file rename -force $temp $resultfile
		}
		if {[llength $afiles]} {file delete {*}$afiles}
	}
}

proc cg_annotate {args} {
	set args [job_init {*}$args]
	cg_annotate_job {*}$args
	job_wait
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base $scriptname
	cg_annotate {*}$argv
}

