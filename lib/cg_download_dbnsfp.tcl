#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_download_dbnsfp {args} {
	set url ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv3.3a.zip
	set basebuild hg38
	set version {}
	set usefields {
		MetaSVM_score MetaSVM_pred
		MetaLR_score MetaLR_pred
		SIFT_score SIFT_pred
		REVEL_score REVEL_rankscore
		VEST3_score VEST3_rankscore
		LRT_score LRT_pred
		MutationTaster_score MutationTaster_pred
		M-CAP_score M-CAP_pred
		FATHMM_score FATHMM_pred
		Polyphen2_HDIV_score Polyphen2_HDIV_pred Polyphen2_HVAR_score Polyphen2_HVAR_pred
		GERP++_NR GERP++_RS phyloP100way_vertebrate SiPhy_29way_pi SiPhy_29way_logOdds LRT_Omega 
	}
	set keep 0
	cg_options download_dbnsfp args {
		-version {set version $value}
		-k {set keep value}
		-fields {set usefields $value}
	} {resultfile build url basebuild} 2 4 {
		download data from dbnsfp
	}
	if {[file exists $resultfile]} {
		putslog "skipping file $resultfile: exists"
		return
	}
	set tail [file tail $url]
	if {$version eq ""} {
		regexp {[0-9][0-9.]*[0-9][a-z]?} $tail version
	}
	set base [file root [gzroot $tail]]
	set tempdir $resultfile.temp
	file mkdir $tempdir
	#
	putslog "Unpacking"
	if {![file exists $tempdir/$tail]} {
		putslog "Downloading $tail"
		wgetfile $url $tempdir/$tail
	}
	set time [clock format [file mtime $tempdir/$tail] -format "%Y-%m-%d %H:%M:%S"]
	exec unzip -o $tempdir/$tail -d $tempdir >@ stdout 2>@ stderr
	#
	putslog "Extracting info"
	set readme [glob $tempdir/dbNSFP*.readme.txt]
	set o [open $tempdir/info w]
	puts $o "= dbnsfp =\n"
	puts $o "== Download info =="
	puts $o dbname\tdbnsfp
	puts $o "version\t$version"
	puts $o "citation\tLiu X, Wu C, Li C and Boerwinkle E. 2016. dbNSFP v3.0: A One-Stop Database of Functional Predictions and Annotations for Human Non-synonymous and Splice Site SNVs. Human Mutation. 37(3):235-241."
	puts $o "website\thttps://sites.google.com/site/jpopgen/dbNSFP"
	puts $o "source\t$url"
	puts $o "time\t$time"
	puts $o "license\tRECEX SHARED SOURCE LICENSE; tl;dr version: use but no redistribution"
	puts $o ""
	puts $o "== Description from Readme =="
	close $o
	exec cat $readme >> $tempdir/info
	file rename -force $tempdir/info [gzroot $resultfile].info
	file copy -force $tempdir/LICENSE.txt [gzroot $resultfile].license
	# data
	set files [ssort -natural [glob $tempdir/dbNSFP*_variant.chr*]]
	set f [open [lindex $files 0]]
	set header [split [string range [gets $f] 1 end] \t]
	close $f
	set nh [list_regsub {\([^)]*\)} $header {}]
	set notpresent [list_lremove $usefields $header]
	if {[llength $notpresent]} {
		puts "fields missing in file [file tail [lindex $files 0]]:$notpresent"
	}
	set usefields [list_common $usefields $header]
	set fields {}
	if {$basebuild eq $build} {
		set query {$pos ne "."}
		lappend fields {chromosome=$chr}
		lappend fields {begin=${pos} - 1}
		lappend fields {end=${pos}}
	} else {
		set query "\$\{${build}_pos\} ne \".\""
		lappend fields "chromosome=\$${build}_chr"
		lappend fields "begin=\$\{${build}_pos\} - 1"
		lappend fields "end=\$\{${build}_pos\}"
	}
	lappend fields {type="snp"} ref alt
	foreach field $usefields {
		regsub -all {[^A-Za-z0-9_]} $field {} newfield
		if {$newfield ne $field} {
			lappend fields "$newfield=\$\{$field\}"
		} else {
			lappend fields $field
		}
	}
	#
	set todo {}
	foreach file $files {
		putslog $file
		lappend todo $file.tsv
		if {[file exists $file.tsv]} continue
		cg select -hp $nh -f $fields -q $query $file $file.tsv.temp
		cg select -s {chromosome end} $file.tsv.temp $file.tsv.temp2
		file delete $file.tsv.temp
		cg groupby {chromosome begin end type} $file.tsv.temp2 $file.tsv.temp3
		file delete $file.tsv.temp2
		cg select -f {
			chromosome begin end type {ref=lindex($ref,0)} alt *
		} $file.tsv.temp3 $file.tsv.temp4
		file delete $file.tsv.temp3
		file rename -force $file.tsv.temp4 $file.tsv
	}
	putslog "joining files"
	cg cat -c 0 {*}$todo > $tempdir/var_hg19_dbnsfp.tsv.temp
	set tempresult $tempdir/var_hg19_dbnsfp.tsv.temp2[file extension $resultfile]
	cg select -s - $tempdir/var_hg19_dbnsfp.tsv.temp $tempresult
	putslog "move result to target"
	# move dbNSFPzip files to target
	file_write [gzroot $resultfile].opt "fields\t{SIFT_score Polyphen2_HDIV_score Polyphen2_HDIV_pred Polyphen2_HVAR_score Polyphen2_HVAR_pred MetaSVM_score MetaSVM_pred MetaLR_score MetaLR_pred LRT_score LRT_pred MutationTaster_score MutationTaster_pred FATHMM_score REVEL_score VEST3_score GERP_NR GERP_RS MetaLR_score MetaLR_pred SiPhy_29way_pi SiPhy_29way_logOdds LRT_Omega ESP_AA_AF ESP_EA_AF}"
	# lz4 already handled in last select
	file rename -force $tempresult $resultfile
	if {!$keep} {file delete -force $tempdir}
}

