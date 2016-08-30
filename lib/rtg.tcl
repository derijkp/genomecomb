#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc rtg2annotvar {file {outfile {}}} {
	if {$outfile eq ""} {
		set outfile [file dir $file]/annotvar-[file tail [file dir $file]].tsv
	}
	catch {close $f}; catch {close $o}
	set f [gzopen $file]
	set keepheader {}
	while {![eof $f]} {
		set line [gets $f]
		if {[string index $line 0] ne "#"} break
		lappend keepheader $line
	}
	if {![llength $keepheader]} {
		set header {name position type reference prediction posterior coverage correction support_statistics}
		set line [split $line \t]
	} else {
		set header [list_pop keepheader]
		set header [split [string range $header 1 end] \t]
	}
	if {$header ne "name position type reference prediction posterior coverage correction support_statistics"} {
		error "header=$header\n$file is not a correct rtg file"
	}
	set o [open $outfile.temp w]
	if {[llength $keepheader]} {
		puts $o [join $keepheader \n]
	}
	puts $o [join {chromosome begin end type reference alleleSeq1 alleleSeq2 posterior coverage correction numA numC numG numT percA percC percG percT nonidentityposterior} \t]
	set poss [list_cor $header {name chromosome position type reference prediction posterior coverage correction nonidentity-posterior support_statistics}]
	set line [split $line \t]
	set num 0
	while 1 {
		incr num
		if {[llength $line] > 3} {
			if {![expr $num%100000]} {putsprogress $num}
			foreach {name position het} $line break
			if {[inlist {e o} $het]} {
				set temp [rtg_line $line $poss]
				puts $o [join $temp \t]
			}
		}
		if {[eof $f]} break
		set line [split [gets $f] \t]
	}
	close $o
	gzclose $f
	file rename -force $outfile.temp $outfile
}

proc cg_rtg2sft {args} {
	if {[llength $args] != 2} {
		errorformat rtg2sft
		exit 1
	}
	foreach {file outfile} $args break
	rtg2annotvar $file $outfile
}

proc rtg_line {cline poss} {
	set line [list_sub $cline $poss]
	set support [lrange $cline [lindex $poss end] end]
	foreach {name chromosome position type reference prediction posterior coverage correction nonidentityposterior} $line break
	if {[inlist {e o} $type]} {
		set type snp
	} elseif {$type ne "="} {
		return {}
	}
	unset -nocomplain stats
	foreach {b n p} $support {
		set stats($b,n) $n
		set stats($b,p) $p
	}
	set chromosome [chr_clip $name]
	set alleles [lrange [split $prediction :] 0 1]
	if {[llength $alleles] == 1} {lappend alleles [lindex $alleles 0]}
	set begin [expr {$position-1}]
	set end $position
	set temp {}
	foreach base {A C G T} {
		lappend temp [get stats($base,n) 0]
	}
	foreach base {A C G T} {
		lappend temp [get stats($base,p) 0.0]
	}
	return [list $chromosome $begin $end $type $reference {*}$alleles $posterior $coverage $correction {*}$temp $nonidentityposterior]
}

proc process_rtgsample {dir destdir {force 0}} {
	set keepdir [pwd]
	set dir [file_absolute $dir]
	set destdir [file_absolute $destdir]
	file mkdir $destdir
	cd $destdir
	putslog "Processing sample $dir"
	set name [file tail $destdir]
	set varfile [gzfile $dir/*snp*.txt]
	# annotated vars file
	if {$force || ![file exists annotvar-$name.tsv]} {
		putslog "Create annotated varfile annotvar-$name.tsv from $varfile"
		if {$force || ![file exists uannotvar-$name.tsv]} {
			rtg2annotvar $varfile uannotvar-$name.tsv.temp
			file rename -force uannotvar-$name.tsv.temp uannotvar-$name.tsv
		}
		putslog "Sorting"
		cg select -s "chromosome begin end" uannotvar-$name.tsv > annotvar-$name.tsv.temp
		file delete uannotvar-$name.tsv
		file rename -force annotvar-$name.tsv.temp annotvar-$name.tsv
	}
	# sample specific filters
	if {$force || ![file exists reg_cluster-$name.tsv]} {
		putslog "Find cluster regions for annotvar-$name.tsv"
		cg clusterregions < annotvar-$name.tsv > temp.tsv
		file rename -force temp.tsv reg_cluster-$name.tsv
	}
	if {$force || ![file exists fannotvar-$name.tsv]} {
		# add filterdata to annotvar
		set todo {}
		lappend todo [list cluster cl reg_cluster-$name.tsv]
		annot_annotvar annotvar-$name.tsv fannotvar-$name.tsv $todo
	}
	putslog "Make allposs files"
	set files [ssort -natural [glob $dir/allpos/chr*/*snps.txt*]]
	file mkdir allpos
	foreach file $files {
		set comments {}
		set f [gzopen $file]
		set buffering [fconfigure $f -buffering]
		fconfigure $f -buffering line
		while {![eof $f]} {
			set line [gets $f]
			if {[string index $line 0] ne "#"} break
			lappend comments $line
			set header $line
		}
		fconfigure $f -buffering $buffering
		set chr [lindex [split $line \t] 0]
		set allposfile allpos/chr${chr}_snps.txt
		if {$force || ![file exists $allposfile.gz]} {
			putslog "Making file $allposfile"
			set o [open $allposfile.temp w]
			foreach l [lrange $comments 0 end-1] {
				puts $o #$l
			}
			puts $o $header
			puts $o $line
			flush $o
			fcopy $f $o
			close $o
			exec bgzip $allposfile.temp
			file rename -force $allposfile.temp.gz $allposfile.gz
		}
	}
	if {$force || ![file exists sreg-$name.tsv]} {
		putslog "Make region file sreg-$name.tsv"
		set files [ssort -natural [glob allpos/chr*snps.txt allpos/chr*snps.txt.gz]]
		file delete sreg-$name.tsv.temp
		set f [open sreg-$name.tsv.temp w]
		puts $f "chromosome\tbegin\tend"
		gzclose $f
		foreach file $files {
			putslog "Processing $file"
			set chr [lindex [split [file tail $file] _] 0]
			set f [gzopen $file]
			set header [tsv_open $f]
			catch {gzclose $f}
			set poscol [lsearch $header position]
			set coveragecol [lsearch $header coverage]
			set cat [gzcat $file]
			exec {*}$cat $file | getregions $chr $poscol $coveragecol 9 1 -1 >> sreg-$name.tsv.temp
		}
		file rename -force sreg-$name.tsv.temp sreg-$name.tsv
	}
	cd $keepdir
}

proc cg_process_rtgsample {args} {
	if {([llength $args] < 2) || ([llength $args] > 3)} {
		errorformat process_rtgsample
		exit 1
	}
	process_rtgsample {*}$args
}

proc annot_rtg_init {dir} {
	global annot
	catch {annot_rtg_close $dir}
	set annot(cov,$dir) {-1 {} 0 {}}
}

proc annot_rtg_get {dir chr begin} {
	global annot
	set chr [chr_clip $chr]
	foreach {curchr chrfile present poss} [get annot(cov,$dir) {{} {} 0}] break
	if {$chr ne $curchr} {
		if {[llength $chrfile]} {
			tsv_index_close $chrfile position
		}
		set chrfile [gzfile $dir/allpos/chr${chr}_snps.txt]
		if {[file exists $chrfile]} {
			tsv_index_open $chrfile position 1
			set present 1
			set header [tsv_index_header $chrfile]
			set poss [list_cor $header {name chromosome position type reference prediction posterior coverage correction nonidentity-posterior support_statistics}]
		} else {
			set present 0
			set poss {}
		}
		set annot(cov,$dir) [list $chr $chrfile $present $poss]
	}
	if {!$present} {return {? ? ? ? ? ? ? ? ? ? ? ? ? ?}}
	ifcatch {tsv_index_get $chrfile position [expr {$begin+1}]} cline -regexp {
		"not found in" {
			set cline {}
		}
	}
	if {![llength $cline]} {return {}}
	return [rtg_line $cline $poss]
}

proc annot_rtg_close {dir} {
	global annot
	foreach {curchr chrfile present} [get annot(cov,$dir) {{} {} 0}] break
	if {[llength $chrfile]} {
		tsv_index_close $chrfile position
	}
	unset annot(cov,$dir)
}

proc rtgregions {cgdir comparfile rtgdir} {
	set cgdir [file_absolute $cgdir]
	set rtgdir [file_absolute $rtgdir]
	set cgsample [file tail $cgdir]
	set rtgsample [file tail $rtgdir]
	if {![file exists $cgdir/reg_rtgnotsame-$cgsample.tsv]} {
		putslog "$cgdir/reg_rtgnotsame-$cgsample.tsv"
		cg select -f {chromosome begin end} -q "!same($cgsample,$rtgsample)" $comparfile $cgdir/reg_rtgnotsame-$cgsample.tsv.temp
		file rename -force $cgdir/reg_rtgnotsame-$cgsample.tsv.temp $cgdir/reg_rtgnotsame-$cgsample.tsv
	}
	if {![file exists $cgdir/reg_posrtg-$cgsample.tsv]} {
		putslog "$cgdir/reg_posrtg-$cgsample.tsv"
		cg regsubtract $rtgdir/sreg-$rtgsample.tsv $cgdir/reg_rtgnotsame-$cgsample.tsv > $cgdir/reg_posrtg-$cgsample.tsv.temp
		file rename -force $cgdir/reg_posrtg-$cgsample.tsv.temp $cgdir/reg_posrtg-$cgsample.tsv
	}
	# make filter for cg by removing good poss in rtg
	if {![file exists $cgdir/reg_rtg-$cgsample.tsv]} {
		putslog "$cgdir/reg_rtg-$cgsample.tsv"
		cg regsubtract $cgdir/sreg-$cgsample.tsv $cgdir/reg_posrtg-$cgsample.tsv > $cgdir/reg_rtg-$cgsample.tsv.temp
		file rename -force $cgdir/reg_rtg-$cgsample.tsv.temp $cgdir/reg_rtg-$cgsample.tsv
	}
	# how much remains after filter
	if {![file exists $cgdir/filteredrtg-$cgsample.tsv]} {
		putslog "$cgdir/filteredrtg-$cgsample.tsv"
		cg regsubtract $cgdir/sreg-$cgsample.tsv $cgdir/reg_rtg-$cgsample.tsv > $cgdir/filteredrtg-$cgsample.tsv.temp
		file rename -force $cgdir/filteredrtg-$cgsample.tsv.temp $cgdir/filteredrtg-$cgsample.tsv
	}
	if {![file exists $cgdir/filteredrtg-$cgsample.covered]} {
		putslog "$cgdir/filteredrtg-$cgsample.covered"
		cg covered $cgdir/filteredrtg-$cgsample.tsv > $cgdir/filteredrtg-$cgsample.covered.temp
		file rename -force $cgdir/filteredrtg-$cgsample.covered.temp $cgdir/filteredrtg-$cgsample.covered
	}
}

proc cg_rtgregions {args} {
	global scriptname action
	if {[llength $args] != 3} {
		puts stderr "format is: $scriptname $action cgdir comparfile rtgdir"
		exit 1
	}
	foreach {cgdir comparfile rtgdir} $args break
	rtgregions $cgdir $comparfile $rtgdir
}
