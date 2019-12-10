proc multidb_dir_open {compar_dir database} {
	if {![file exists $compar_dir/vars.tsv]} {
		file_write $compar_dir/vars.tsv [join {chromosome begin end type ref alt id} \t]\n
		file_write $compar_dir/vars.tsv.maxid 0
		file_write $compar_dir/vars.tsv.count 0
	} elseif {![file exists $compar_dir/vars.tsv.maxid] || ![file exists $compar_dir/vars.tsv.count]} {
		set temp [cg select -g all -gc {max(id),count(id)} $compar_dir/vars.tsv]
		set maxid [lindex $temp end-1]
		if {![isint $maxid]} {set maxid 0}
		set count [lindex $temp end]
		if {![isint $count]} {set count 0}
		file_write $compar_dir/vars.tsv.maxid $maxid
		file_write $compar_dir/vars.tsv.count $count
	}
	if {![file exists $compar_dir/analysis.tsv]} {
		set fields {id experiment ngsseq gentli_ngsseq reference mapper varcall individual sample}
		file_write $compar_dir/analysis.tsv [join $fields \t]\n
		file_write $compar_dir/analysis.tsv.maxid 0
		file_write $compar_dir/analysis.tsv.count 0
	} elseif {![file exists $compar_dir/analysis.tsv.maxid] || ![file exists $compar_dir/analysis.tsv.count]} {
		set temp [cg select -g all -gc {max(id),count(id)} $compar_dir/analysis.tsv]
		set maxid [lindex $temp end-1]
		if {![isint $maxid]} {set maxid 0}
		set count [lindex $temp end]
		if {![isint $count]} {set count 0}
		file_write $compar_dir/analysis.tsv.maxid $maxid
		file_write $compar_dir/analysis.tsv.count $count
	}
	set genofields {var analysis sequenced zyg alleleSeq1 alleleSeq2 quality coverage}
	if {[file exists $compar_dir/geno.tsv]} {
		set fields [cg select -h $compar_dir/geno.tsv]
		set genofields [list_union $genofields $fields]
	}
	return $genofields
}

proc multidb_dir_import {compar_dir} {
	file mkdir $compar_dir/old
	puts "Importing for $compar_dir"
	if {[file exists $compar_dir/analysis.tsv.insert]} {
		puts "Importing analysis.tsv.insert"
		cg cat -m -c 0 $compar_dir/analysis.tsv $compar_dir/analysis.tsv.insert > $compar_dir/analysis.tsv.temp
		set count1 [file_read $compar_dir/analysis.tsv.count]
		set count2 [file_read $compar_dir/analysis.tsv.insert.count]
		file_write $compar_dir/analysis.tsv.temp.count [expr {$count1 + $count2}]
		file rename -force -- $compar_dir/analysis.tsv.temp $compar_dir/analysis.tsv
		file rename -force -- $compar_dir/analysis.tsv.temp.count $compar_dir/analysis.tsv.count
		file rename -force -- $compar_dir/analysis.tsv.insert.maxid $compar_dir/analysis.tsv.maxid
		file rename -force -- $compar_dir/analysis.tsv.insert $compar_dir/old
		file rename -force -- $compar_dir/analysis.tsv.insert.count $compar_dir/old
	}
	if {[file exists $compar_dir/geno.tsv.insert]} {
		puts "Importing geno.tsv.insert"
		if {[file exists $compar_dir/geno.tsv]} {
			cg cat -m -c 0 $compar_dir/geno.tsv $compar_dir/geno.tsv.insert > $compar_dir/geno.tsv.temp
			set count1 [file_read $compar_dir/geno.tsv.count]
			set count2 [file_read $compar_dir/geno.tsv.insert.count]
			file_write $compar_dir/geno.tsv.temp.count [expr {$count1 + $count2}]
			file rename -force -- $compar_dir/geno.tsv.temp $compar_dir/geno.tsv
			file rename -force -- $compar_dir/geno.tsv.temp.count $compar_dir/geno.tsv.count
		} else {
			file copy -force $compar_dir/geno.tsv.insert $compar_dir/geno.tsv
			file copy -force $compar_dir/geno.tsv.insert.count $compar_dir/geno.tsv.count
		}
		file rename -force -- $compar_dir/geno.tsv.insert $compar_dir/old
		file rename -force -- $compar_dir/geno.tsv.insert.count $compar_dir/old
	}
	if {[file exists $compar_dir/vars.tsv.new]} {
		puts "Importing geno.tsv.insert"
		file rename -force -- $compar_dir/vars.tsv $compar_dir/old
		file rename -force -- $compar_dir/vars.tsv.count $compar_dir/old
		file rename -force -- $compar_dir/vars.tsv.maxid $compar_dir/old
		file rename -force -- $compar_dir/vars.tsv.new $compar_dir/vars.tsv
		file rename -force -- $compar_dir/vars.tsv.new.count $compar_dir/vars.tsv.count
		file rename -force -- $compar_dir/vars.tsv.new.maxid $compar_dir/vars.tsv.maxid
		file rename -force -- $compar_dir/vars.tsv.insert $compar_dir/old
		file rename -force -- $compar_dir/vars.tsv.insert.count $compar_dir/old
	}
}

