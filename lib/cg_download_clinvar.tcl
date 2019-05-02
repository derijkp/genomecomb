proc cg_download_clinvar {args} {
	set keep 0
	set papuurl {}
	cg_options download_clinvar args {
		-k {set keep $value}
	} {resultfile build url papuurl} 2 4 {
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
	# init and download
	if {![regexp {clinvar_(.*)\.vcf} [file tail $url] temp version]} {
		set version [timestamp]
	}
	set tempdir $resultfile.temp
	set tempresult $tempdir/result.tsv
	file mkdir $tempdir
	set tail [file tail $url]
	set paputail [file tail $papuurl]
	wgetfile $url $tempdir/$tail
	wgetfile ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/README.txt $tempdir/README.txt
	set tsvfile $tempdir/[file root [gzroot $tail]].tsv
	#
	# opt
	putslog "Write opt and info"
	file_write [gzroot $resultfile].opt "fields\t{CLNDISDB CLNDN CLNSIG}\nheaderfields\t{clinvar_acc clinvar_disease clinvar_sig}\n"
	# info
	file_write [gzroot $resultfile].info [deindent [subst {
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
	}]]\n
	exec cat $tempdir/README.txt >> [gzroot $resultfile].info
	# var file
	putslog "Convert $tail to $tsvfile"
	cg vcf2tsv -split ori -typelist ". CAF R" $tempdir/$tail $tsvfile.temp
	if {$papuurl ne ""} {
		wgetfile $papuurl $tempdir/$paputail
		cg vcf2tsv -split ori -typelist ". CAF R" $tempdir/$paputail $tsvfile-papu.temp
		cg cat -c f $tsvfile.temp $tsvfile.temp $tsvfile-papu.temp | cg select -s - > $tsvfile
	} else {
		file rename $tsvfile.temp $tsvfile
	}
	# check version
	set f [open $tsvfile]
	set header [tsv_open $f comment]
	close $f
	if {$build eq "hg19" && ![regexp reference.GRCh37 $comment]} {
		error "$url is from a different reference genome version than $build"
	}
	if {$build eq "hg38" && ![regexp reference.GRCh38 $comment]} {
		error "$url is from a different reference genome version than $build"
	}
	# link properly to alleles
	set f [open $tsvfile]
	set header [tsv_open $f comment]
	set apos [lsearch $header alt]
	set clnapos [lsearch $header CLNALLE]
	if {$clnapos != -1} {
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
	} else {
		# collapse
		cg collapsealleles $tsvfile > $tempresult.temp2
	}
	compress $tempresult.temp2 $resultfile
	if {!$keep} {file delete -force $tempdir}
}


