proc rtg2annotvar {file {outfile {}}} {
	if {$outfile eq ""} {
		set outfile [file dir $file]/annotvar-[file tail [file dir $file]].tsv
	}
	catch {close $f}; catch {close $o}
	set f [rzopen $file]
	set o [open $outfile.temp w]
	puts $o [join {chromosome begin end type alleleSeq1 alleleSeq2 posterior coverage correction numA numC numG numT percA percC percG percT reference} \t]
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
	set line [split $line \t]
	set num 0
	while 1 {
		incr num
		if {![expr $num%100000]} {putslog $num}
		foreach {name position het reference prediction posterior coverage correction} $line break
		if {[inlist {e o} $het]} {
			unset -nocomplain stats
			foreach {b n p} [lrange $line 8 end] {
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
			puts $o "$chromosome\t$begin\t$end\tsnp\t[join $alleles \t]\t$posterior\t$coverage\t$correction\t[join $temp \t]\t$reference"
		}
		if {[eof $f]} break
		set line [split [gets $f] \t]
	}
	close $o
	close $f
	file rename -force $outfile.temp $outfile
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
	# fake region file
	if {$force || ![file exists sreg-$name.tsv]} {
		set chrs [list_fill 22 1 1]
		lappend chrs M X Y
		set o [open sreg-$name.tsv w]
		puts $o	[join {chromosome begin end} \t]
		foreach chr $chrs {
			puts $o $chr\t1\t1000000000
		}
		close $o
	}
	cd $keepdir
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

