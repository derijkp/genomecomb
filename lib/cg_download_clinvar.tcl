proc cg_download_clinvar {args} {
	set keep 0
	cg_options download_clinvar args {
		-k {set keep $value}
	} {resultfile build url} 2 3 {
		download data from clinvar
	}
	if {![info exists url]} {
		if {$build eq "hg19"} {
			set url ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
		} elseif {$build eq "hg38"} {
			set url ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
		} else {
			error "Do not know url for build $build, please provide the url"
		}
	}
	regsub {/[^/]+$} $url {} baseurl
	set tempdir $resultfile.temp
	set tempresult $tempdir/result.tsv
	set version [timestamp]
	file mkdir $tempdir
	set tail [file tail $url]
	wgetfile $url $tempdir/$tail
	wgetfile ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/README.txt $tempdir/README.txt
	set tsvfile $tempdir/[file root [gzroot $tail]].tsv
	#
	putslog "Convert $tail to $tsvfile"
	cg vcf2tsv -split ori -typelist ". CAF R" $tempdir/$tail $tsvfile
	# check version
	set f [open $tsvfile]
	set header [tsv_open $f comment]
	close $f
	if {$build eq "hg19" && ![regexp reference=GRCh37 $comment]} {
		error "$url is from a different reference genome version"
	}
	if {$build eq "hg38" && ![regexp reference=GRCh38 $comment]} {
		error "$url is from a different reference genome version"
	}
	# link properly to alleles
	set f [open $tsvfile]
	set header [tsv_open $f comment]
	set apos [lsearch $header alt]
	set clnapos [lsearch $header CLNALLE]
	set cposs [list_find -regexp $header ^CLN]
	set o [open $tempresult.temp w]
	puts -nonewline $o $comment
	puts $o [join $header \t]
	while 1 {
		if {[gets $f line] == -1} break
		set line [split $line \t]
		set alts [split [lindex $line $apos] ,]
		set clna [split [lindex $line $clnapos] ,]
		if {$clna == 1} {
			puts $o [join $line \t]
		} else {
			set pos 0
			foreach a $clna {
				if {$a == -1 || $a == 0} continue
				incr a -1
				set result $line
				lset result $apos [lindex $alts $a]
				foreach cpos $cposs {
					set list [split [lindex $line $cpos] ,] 
					lset result $cpos [lindex $list $pos]
				}
				puts $o [join $result \t]
				incr pos
			}
		}
	}
	close $o
	close $f
	# collapse
	cg collapsealleles $tempresult.temp > $tempresult.temp2
	file rename $tempresult.temp2 $tempresult
	#
	putslog "Write opt and info"
	file_write $tempresult.opt "fields\t{CLNACC CLNDBN}\nheaderfields\t{clinvar_acc clinvar_disease}\n"
	# info
	file_write $tempresult.info [deindent [subst {
		clinvar variants (version $version)
		================
		
		Download info
		-------------
		dbname	clinvar (clinvar variants)
		version	$version
		website	http://www.ncbi.nlm.nih.gov/clinvar
		source	$url
		time	[timestamp]
		summary	ClinVar aggregates information about genomic variation and its relationship to human health.
		
		README
		------
	}]]
	exec echo $tempdir/README.txt >> $tempresult.info
	file rename -force $tempresult.opt $resultfile.opt
	file rename -force $tempresult.info $resultfile.info
	file rename -force $tempresult $resultfile
	if {!$keep} {file delete -force $tempdir}
}


