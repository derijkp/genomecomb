proc cg_download_genes {args} {
	cg_options annotate args {
		-geneidcol {set geneidcol $value}
	} {resultfile build geneset}
	if {[file isfile $resultfile]} {
		error "The file '$resultfile' already exists..."
	}
	set resulttail [file tail $resultfile]
	set temp $resultfile.temp
	file mkdir $temp
	set ucscfile $temp/ucsc_${build}_${geneset}.tsv
	cg_download_ucsc $ucscfile $build $geneset
	if {$geneset in "ensGene knownGene" && ![info exists geneidcol]} {
		unset -nocomplain a
		if {$geneset eq "ensGene"} {
			if {![catch {
				cg_download_ucsc $temp/ensemblToGeneName.tsv ${build} ensemblToGeneName
			}]} {
				array set a [split [string trim [file_read $temp/ensemblToGeneName.tsv]] "\n\t"]
			}
		} else {
			if {![catch {
				cg_download_ucsc $temp/kgXref.tsv ${build} kgXref
			}]} {
				array set a [split [string trim [cg select -f {kgID geneSymbol} $temp/kgXref.tsv]] "\n\t"]
			}
		}
		putslog "Adding geneid"
		catch {close $f} ; catch {close $o}
		set f [open $ucscfile]
		set o [open $temp/$resulttail.temp2 w]
		set header [split [gets $f] \t]
		set namepos [lsearch $header name]
		lappend header geneid
		puts $o [join $header \t]
		while {[gets $f line] != -1} {
			set line [split $line \t]
			set name [lindex $line $namepos]
			set geneid [get a($name) $name]
			lappend line $geneid
			puts $o [join $line \t]
		}
		close $o
		close $f
	        cg select -s - -f {chrom start end strand geneid *} $temp/$resulttail.temp2 {*}[compresspipe $resulttail 12] > $temp/$resulttail
	} else {
		if {![info exists geneidcol]} {
			if {$geneset in "genscan acembly"} {
				set geneidcol name
			} else {
				set geneidcol name2
			}
		}
		set fields [cg select -h $ucscfile]
		if {[lsearch $fields $geneidcol] == -1} {
			set geneidcol [lindex [list_common {name2 name geneid id} $fields] 0]
		}
	        cg select -s - -f [list chrom start end strand "geneid=\$$geneidcol" *] $ucscfile {*}[compresspipe $resulttail 12] > $temp/$resulttail
	}
	# move to results
	putslog "move results to $resultfile and [gzroot $resultfile].info"
	file rename -force -- $ucscfile.info [gzroot $resultfile].info
	compress $temp/$resulttail $resultfile
	file delete -force $temp
}
