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
	foreach chr {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X M Y} {
		set file $dir/sv/$name-$chr-paired.tsv
		set root [file root $file]
		if {$force || ![file exists $file]} {
			set file $file.gz
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
		if {$force || [file extension $file] ne ".gz"} {
			puts "bgzipping $file"
			exec bgzip $file
		}
	}
	puts stderr "Finished finding sv in $dir"
	cd $keepdir
}
