#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc process_sv {cgdir dir dbdir {force 0}} {
	set keepdir [pwd]
	set cgdir [file normalize $cgdir]
	set dir [file normalize $dir]
	set dbdir [file normalize $dbdir]
	set name [file tail $dir]
	file mkdir $dir/sv
	cd $dir/sv

	if {$force || ![file exists $dir/sv/${name}_map2sv_sort_FINISHED]} {
		cg map2sv $cgdir $dir/sv/$name
	}
	set resultfiles {}
	foreach chr {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X M Y} {
		lappend resultfiles $dir/sv/$name-$chr-paired-sv.tsv
		set file [gzfile $dir/sv/$name-$chr-paired.tsv]
		set root [file root $file]
		if {$force || ![file exists $file]} {
			set file $file.rz
		}
		if {$force || ![file exists $root.tsv.end1_index]} {
			puts "Indexing $file"
			cg tsv_index end1 $file
		}
		if {$force || ![file exists $root.tsv.numinfo]} {
			puts "Info on $file"
			cg svinfo $file
		}
		if {$force || ![file exists $root-sv.tsv]} {
			puts "svfind $file"
			cg svfind $file $dbdir/reg_hg18_simpleRepeat.tsv
		}
		if {$force || [file extension $file] ne ".rz"} {
			puts "razipping $file"
			cg_razip $file
		}
	}
	cg cat {*}$resultfiles > $dir/sv-$name.tsv
	putslog "Done: finished finding sv in $dir"
	cd $keepdir
}

proc cg_process_sv {args} {
	global scriptname action
	if {([llength $args] < 2) || ([llength $args] > 3)} {
		puts stderr "format is: $scriptname $action sampledir destdir dbdir ?force?"
		puts stderr " - processes sv for one sample directory."
		puts stderr " - By default, only files that are not present already will be created."
		puts stderr " -When force is given as a parameter, everything will be recalculated and overwritten."
		exit 1
	}
	foreach {dir destdir dbdir force} $args break
	switch $force {
		force {set force 1}
		"" {set force 0}
		default {error "unrecognized option $force"}
	}
	process_sv $dir $destdir $dbdir $force
}
