#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_download_kaviar {args} {
	set version {}
	set build hg19
	set keep 0
	cg_options download_dbnsfp args {
		-version {set version $value}
		-k {set keep value}
	} {resultfile build url} 1 3 {
		download data from dbnsfp
	}
	if {[file exists $resultfile]} {
		putslog "skipping file $resultfile: exists"
		return
	}
	if {![info exists url]} {
		if {$build eq "hg19"} {
			set url http://s3-us-west-2.amazonaws.com/kaviar-160204-public/Kaviar-160204-Public-hg19-trim.vcf.tar
		} elseif {$build eq "hg38"} {
			set url http://s3-us-west-2.amazonaws.com/kaviar-160204-public/Kaviar-160204-Public-hg38-trim.vcf.tar
		} else {
			error "unsupported build, specify url"
		}
	}
	set tail [file tail $url]
	set rtail [file tail $resultfile]
	if {![regexp -- "-${build}-" $tail]} {
		error "url $url does not seem to be of build $build"
	}
	if {$version eq ""} {
		regexp {[0-9][0-9.]*[0-9][a-z]?} $tail version
	}
	set base [file root [gzroot $tail]]
	set tempdir $resultfile.temp
	file mkdir $tempdir
	if {![file exists $tempdir/$tail]} {
		putslog "Downloading $tail"
		wgetfile $url $tempdir/$tail
	}
	#
	putslog "Extracting $tail"
	exec tar xvf $tempdir/$tail --directory $tempdir
	set vcffile [glob $tempdir/*/*/Kaviar*vcf.gz]
	cg vcf2tsv $vcffile $tempdir/$rtail.temp
	cg select -f {chromosome begin end type ref alt {freqp=vformat("%.5f",100.0 @* $frequency)} allelecount totalallelecount} $tempdir/$rtail.temp $tempdir/$rtail.temp2
	file delete $tempdir/$rtail.temp
	#
	putslog "Extracting info"
	set readme [glob $tempdir/*/*/README_Kaviar*]
	set o [open $tempdir/info w]
	puts $o dbname\tkaviar
	puts $o "version\t$version"
	puts $o "citation\tGlusman G, Caballero J, Mauldin DE, Hood L and Roach J (2011) KAVIAR: an accessible system for testing SNV novelty. Bioinformatics, doi: 10.1093/bioinformatics/btr540"
	puts $o "website\thttp://db.systemsbiology.net/kaviar"
	puts $o "source\t$url"
	puts $o "time\t[timestamp]"
	puts $o ""
	close $o
	exec cat $readme >> $tempdir/info
	file rename $tempdir/info $resultfile.info
	#
	putslog "move result to target"
	file rename -force $tempdir/$rtail.temp2 $resultfile
	if {!$keep} {file delete -force $tempdir}
}
