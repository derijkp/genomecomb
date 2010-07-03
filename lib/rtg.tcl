proc rtg2annotvar {file {outfile {}}} {
	if {$outfile eq ""} {
		set outfile [file dir $file]/annotvar-[file tail [file dir $file]].tsv
	}
	catch {close $f}; catch {close $o}
	set f [rzopen $file]
	set o [open $outfile.temp w]
	puts $o [join {chromosome begin end type alleleSeq1 alleleSeq2 posterior coverage correction numA numC numG numT percA percC percG percT nonidentityposterior reference} \t]
	while {![eof $f]} {
		set line [gets $f]
		if {[string index $line 0] ne "#"} break
		set header $line
	}
	if {![info exists header]} {
		set header {name position type reference prediction posterior coverage correction support_statistics}
		set line [split $line \t]
	} else {
		set header [split [string range $header 1 end] \t]
	}
	if {$header ne "name position type reference prediction posterior coverage correction support_statistics"} {
		error "header=$header\n$file is not a correct rtg file"
	}
	set poss [list_cor $header {name chromosome position type reference prediction posterior coverage correction nonidentity-posterior support_statistics}]
	set line [split $line \t]
	set num 0
	while 1 {
		incr num
		if {![expr $num%100000]} {putslog $num}
		foreach {name position het} $line break
		if {[inlist {e o} $het]} {
			set temp [rtg_line $line $poss]
			puts $o [join $temp \t]
		}
		if {[eof $f]} break
		set line [split [gets $f] \t]
	}
	close $o
	close $f
	file rename -force $outfile.temp $outfile
}

proc rtg_line {cline poss} {
	set line [list_sub $cline $poss]
	set support [lrange $cline [lindex $poss end] end]
	foreach {name chromosome position type reference prediction posterior coverage correction nonidentityposterior} $line break
	if {[inlist {e o} $type]} {set type snp}
	unset -nocomplain stats
	foreach {b n p} $support {
		set stats($b,n) $n
		set stats($b,p) $p
	}
	regsub ^chr $name {} chromosome
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
	return [list $chromosome $begin $end $type {*}$alleles $posterior $coverage $correction {*}$temp $nonidentityposterior $reference]
}

proc process_rtgsample {dir dbdir {force 0}} {
	set keepdir [pwd]
	set dir [file normalize $dir]
	cd $dir
	puts stderr "Processing sample $dir"
	set name [file tail $dir]
	set varfile [lindex [glob $dir/ori/*unfiltered.snp*] 0]
	# annotated vars file
	if {$force || ![file exists annotvar-$name.tsv]} {
		puts stderr "Create annotated varfile annotvar-$name.tsv from $varfile"
		if {$force || ![file exists unsorted.tsv]} {
			rtg2annotvar $varfile unsorted.tsv
		}
		puts stderr "Sorting"
		cg select -s "chromosome begin end" unsorted.tsv > temp.tsv
		file rename -force temp.tsv annotvar-$name.tsv
	}
	# sample specific filters
	if {$force || ![file exists reg_cluster-$name.tsv]} {
		puts stderr "Find cluster regions for annotvar-$name.tsv"
		cg clusterregions < annotvar-$name.tsv > temp.tsv
		file rename -force temp.tsv reg_cluster-$name.tsv
	}
	if {$force || ![file exists fannotvar-$name.tsv]} {
		# add filterdata to annotvar
		set todo {}
		lappend todo [list cluster cl reg_cluster-$name.tsv]
		lappend todo [list trf trf $dbdir/regdb-simple_repeats.tsv]
		lappend todo [list str str $dbdir/regdb-microsatelite.tsv]
		lappend todo [list segdup sd $dbdir/regdb-segdups.tsv]
		lappend todo [list selfchain sc $dbdir/regdb-selfchain.tsv]
		lappend todo [list repeat rp $dbdir/regdb-repeatmasker.tsv]
		lappend todo [list rna rna $dbdir/regdb-rnagenes.tsv]
		foreach file [lsort -dictionary [glob $dbdir/checked*.tsv]] {
			set value [lindex [split [file root [file tail $file]] _] 1]
			if {$value eq ""} {set value checked}
			lappend todo [list checked $value $file]
		}
		annot_annotvar annotvar-$name.tsv fannotvar-$name.tsv $todo
	}
	if {$force || ![file exists sreg-$name.tsv]} {
		puts stderr "Make region file sreg-$name.tsv"
		set files [lsort -dict [glob allpos/chr*_snps.txt*]]
		file delete temp.tsv
		set f [open temp.tsv w]
		puts $f "chromosome\tbegin\tend"
		close $f
		foreach file $files {
			set chr [lindex [split [file tail $file] _] 0]
			puts stderr "Processing $file"
			set f [rzopen $file]
			while {![eof $f]} {
				set line [gets $f]
				if {[string index $line 0] ne "#"} break
				set header $line
			}
			catch {close $f}
			set header [string range $header 1 end]
			set poscol [lsearch $header position]
			set coveragecol [lsearch $header coverage]
			if {[inlist {.rz .gz} [file extension $file]]} {set cat zcat} else {set cat cat}
			exec $cat $file | getregions $chr $poscol $coveragecol 9 1 -1 >> temp.tsv
		}
		file rename -force temp.tsv sreg-$name.tsv
	}
	cd $keepdir
}

proc annot_rtg_init {dir} {
	global annot
	catch {annot_rtg_close $dir}
	set annot(cov,$dir) {-1 {} 0 {}}
}

proc annot_rtg_get {dir chr begin} {
	global annot
	foreach {curchr chrfile present poss} [get annot(cov,$dir) {{} {} 0}] break
	if {$chr ne $curchr} {
		if {[llength $chrfile]} {
			tsv_index_close $chrfile position
		}
		set chrfile [lindex [glob -nocomplain $dir/allpos/chr${chr}_snps.txt $dir/allpos/chr${chr}_snps.txt.rz $dir/allpos/chr${chr}_snps.txt.gz] 0]
		if {[llength $chrfile]} {
			tsv_index_open $chrfile position 1
			set present 1
		} else {
			set present 0
		}
		set header [tsv_index_header $chrfile]
		set poss [list_cor $header {name chromosome position type reference prediction posterior coverage correction nonidentity-posterior support_statistics}]
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

if 0 {

	lappend auto_path ~/dev/completegenomics/lib
	package require Extral
package require Tclx
signal -restart error SIGINT
	set dbdir /complgen/refseq
	set file /complgen/rtg102/GS00102-DNA-D06-unfiltered.snps.gz
	set dir [file dir $file]

}

if 0 {
	set cgdir /complgen/GS102
	set comparfile /complgen/multicompar/compar.tsv
	set rtgdir /complgen/rtg102
	submit.tcl cg rtgregions /complgen/GS102 /complgen/multicompar/compar.tsv /complgen/rtg102
	submit.tcl cg rtgregions /complgen/GS103 /complgen/multicompar/compar.tsv /complgen/rtg103
	submit.tcl cg rtgregions /complgen/GS100 /complgen/multicompar/compar.tsv /complgen/rtg100
	submit.tcl cg rtgregions /complgen/GS101A01 /complgen/multicompar/compar.tsv /complgen/rtg101A01
	submit.tcl cg rtgregions /complgen/GS101A02 /complgen/multicompar/compar.tsv /complgen/rtg101A02
}

proc rtgregions {cgdir comparfile rtgdir} {
	set cgsample [file tail $cgdir]
	set rtgsample [file tail $rtgdir]
	if {![file exists $cgdir/reg_rtgnotsame-$cgsample.tsv]} {
		puts stderr "$cgdir/reg_rtgnotsame-$cgsample.tsv"
		cg select -f {chromosome begin end} -q "!\$same($cgsample,$rtgsample)" $comparfile $cgdir/temp.tsv
		file rename -force $cgdir/temp.tsv $cgdir/reg_rtgnotsame-$cgsample.tsv
	}
	if {![file exists $cgdir/reg_posrtg-$cgsample.tsv]} {
		puts stderr "$cgdir/reg_posrtg-$cgsample.tsv"
		cg regsubtract $rtgdir/sreg-$rtgsample.tsv $cgdir/reg_rtgnotsame-$cgsample.tsv > $cgdir/temp.tsv
		file rename -force $cgdir/temp.tsv $cgdir/reg_posrtg-$cgsample.tsv
	}
	# make filter for cg by removing good poss in rtg
	if {![file exists $cgdir/reg_rtg-$cgsample.tsv]} {
		puts stderr "$cgdir/reg_rtg-$cgsample.tsv"
		cg regsubtract $cgdir/sreg-$cgsample.tsv $cgdir/reg_posrtg-$cgsample.tsv > $cgdir/temp.tsv
		file rename -force $cgdir/temp.tsv $cgdir/reg_rtg-$cgsample.tsv
	}
	# how much remains after filter
	if {![file exists $cgdir/filteredrtg-$cgsample.tsv]} {
		puts stderr "$cgdir/filteredrtg-$cgsample.tsv"
		cg regsubtract $cgdir/sreg-$cgsample.tsv $cgdir/reg_rtg-$cgsample.tsv > $cgdir/temp.tsv
		file rename -force $cgdir/temp.tsv $cgdir/filteredrtg-$cgsample.tsv
	}
	if {![file exists $cgdir/filteredrtg-$cgsample.covered]} {
		puts stderr "$cgdir/filteredrtg-$cgsample.covered"
		cg covered $cgdir/filteredrtg-$cgsample.tsv > $cgdir/temp.tsv
		file rename $cgdir/temp.tsv $cgdir/filteredrtg-$cgsample.covered
	}
}
