proc cg_download_ucscinfo {args} {
	cg_options download_ucscinfo args {
		-c {set continue $value}
	} {resultfile build dbname} 3 3 {
		download info on UCSC table dbname
	}
	set temp $resultfile.temp
	file mkdir $temp
	# do documentation
	# ----------------
	set timestamp [timestamp]
	set version	$timestamp
	set lastupdated	?
	if {![file exists $temp/$dbname.html]} {
		puts "Downloading $dbname.html ....."
		wgetfile http://genome.ucsc.edu/cgi-bin/hgTrackUi?db=$build&g=$dbname $temp/$dbname.html
	} else {
		puts "Skipping download $temp/$dbname.html (already there)"
	}
	puts "Making $resultfile"
	set c [file_read $temp/$dbname.html]
	if {[regexp {HGERROR-START} $c]} {
		puts "Skipping $resultfile: Could not download info"
		return
	}
	set title $dbname
	if {![regexp {Data version:(</B>)? *([^<]+)} $c t t version]} {
		regexp {GENCODE Version *([^ ]+)} $c t version
	}
	regexp {Data last updated:(&nbsp;)?(</B>)? *([^<]+)} $c t t t lastupdated
	regexp {<TITLE>([^<]+)</TITLE>} $c t title
	regexp {<B style='font-family:serif; font-size:200%;'>([^<]+)</B>} $c t title
#	regsub {^.*<H2>Description</H2>} $c {} c
	if {![regsub {^.*<A NAME='TRACK_HTML'></A>} $c {} c]} {
		set c {}
	} else {
		regsub {</td><td nowrap>.*$} $c {} c
		regsub -all {</[Hh]([0-5])>\n+} $c "</H\\1>\n" c
		regsub -all {<[Hh]([0-5])>} $c "<H\\1>" c
		set c [string_change $c [list <P> {} </P> {} <p> {} </p> {} <H2> "== " </H2> " ==" <H3> "=== " </H3> " ===" <H4> "==== " </H4> " ====" <H5> "===== " </H5> " ====="]]
	}
file_write $temp/$dbname.info [subst [string trim {
= $dbname: [string trim $title] =

== Download info ==
dbname	$dbname
version	$version
lastupdated	$lastupdated
source	http://genome.ucsc.edu/cgi-bin/hgTrackUi?db=$build&g=$dbname
time	$timestamp

[string trim $c]

== Category ==
Annotation
}]]
	file rename -force $temp/$dbname.info $resultfile
	file delete -force $temp
}

proc cg_download_ucsc {args} {
	set continue 1
	cg_options download_ucsc args {
		-c {set continue $value}
	} {resultfile build dbname} 3 3 {
		download data from UCSC table dbname as a tsv file
	}
	if {[file isfile $resultfile]} {
		if {$continue} {
			putslog "'$resultfile' already exists: skipping download"
			return
		}
		error "The file '$resultfile' already exists..."
	}
	set resulttail [file tail $resultfile]
	set temp $resultfile.temp
	if {!$continue} {file delete -force $temp}
	file mkdir $temp
	# download documentation
	# ----------------
	cg_download_ucscinfo [gzroot $resultfile].info $build $dbname
	#
	# download data
	set single 1
	set sqlfile $dbname.sql
	if {![file exists $temp/$dbname.txt.gz]} {
		putslog "Downloading $dbname.txt.gz ....."
		wgetfile ftp://hgdownload.cse.ucsc.edu/goldenPath/$build/database/$dbname.txt.gz $temp/$dbname.txt.gz
	} else {
		putslog "Skipping download $dbname.txt.gz (already there)"
	}
	if {![file exists $temp/$dbname.txt.gz]} {
		set single 0
		wgetfiles ftp://hgdownload.cse.ucsc.edu/goldenPath/$build/database/chr*_$dbname.txt.gz $temp
		set files [ssort -natural [glob -nocomplain $temp/chr*_$dbname.txt.gz]]
		if {![llength $files]} {
			error "Could not download $dbname.txt.gz"
		}
		set sqlfile [file root [gzroot [file tail [lindex $files 0]]]].sql
	}
	if {![file exists $temp/$sqlfile]} {
		putslog "Downloading $dbname.sql ....."
		wgetfile ftp://hgdownload.cse.ucsc.edu/goldenPath/$build/database/$sqlfile $temp/$sqlfile
	} else {
		putslog "Skipping download $dbname.sql (already there)"
	}
	if {![file exists $temp/$sqlfile]} {
		error $errmsg
	}
	# header
	# ------
	putslog "Finding header...."
	if {[catch "open $temp/$sqlfile r" fileid]} {
		error "Could not open sql file - $fileid"
	}
	set line [gets $fileid]
	set header ""
	while {![eof $fileid]} {
		if {[lsearch $line "CREATE"] != -1} {
			set line [gets $fileid]
			set head [string trim [lindex $line 0] (`) ]
			while {![inlist {KEY PRIMARY UNIQUE )} $head]} {	
				if {$head eq ""} break
				lappend header $head
				if {[eof $fileid]} break
				set line [gets $fileid]
				set head [string trim [lindex $line 0] (`) ]
			}
		}
		set line [gets $fileid]
	}
	close $fileid
	set header [list_change $header {
		chromStart start chromEnd end
		genoName chrom genoStart start genoEnd end
		txStart start txEnd end
		tName chrom tStart start tEnd end
		repName name swScore score
	}]
	# write file
	# ----------
	set file_outid [open $temp/u_$dbname.tsv w]
	puts $file_outid [join $header \t]
	close $file_outid
	if {$single} {
		putslog "Gunzipping $dbname.txt.gz ...."
		if {[catch "exec zcat -f $temp/$dbname.txt.gz >> $temp/u_$dbname.tsv" errmsg]} {
			error "Gunzip failed - $errmsg"
		}
	} else {
		foreach file $files {
			putslog "Gunzipping $file ...."
			if {[catch "exec zcat -f $file >> $temp/u_$dbname.tsv" errmsg]} {
				error "Gunzip failed - $errmsg"
			}
		}
	}
	putslog "Sorting $resultfile ...."
	set fields [list_common {chrom chromosome begin start end stop name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts} $header]
	lappend fields {*}[list_lremove $header $fields]
	putslog "Writing $temp/$resulttail"
	exec cg select -s - -f $fields $temp/u_$dbname.tsv {*}[compresspipe $resulttail 12] > $temp/$resulttail
	# move to result
	putslog "move results to $resultfile and $resultfile.info"
	compress $temp/$resulttail $resultfile
	file delete -force $temp
}
