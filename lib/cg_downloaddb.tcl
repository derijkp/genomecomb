package require Extral

##############################################################
# 
#  This is a script that downloads a number of databases
#  
#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# parts by Annelies Cassiers (VIB - University of Antwerp)
#
# cg_downloaddb way is deprecated (old code organization), still kept for compat ...
# all to cg_download_<name> resultfile ...
#
##############################################################

proc downloaddb {path build dbname} {
	file mkdir $path
	set filename $path/$build/ucsc_${build}_$dbname.tsv
	set temp $filename.temp
	file delete -force $temp
	file mkdir $temp
	if {[file isfile $filename]} {
		puts "The file '$filename' already exists..."
		puts "Skipping the download..."
		puts "----------------------------------------------------"
		return
	}
	catch {file mkdir [file dir $filename]}
	set single 1
	set sqlfile $dbname.sql
	if {![file exists $temp/$dbname.txt.gz]} {
		puts "Downloading $dbname.txt.gz ....."
		wgetfile ftp://hgdownload.cse.ucsc.edu/goldenPath/$build/database/$dbname.txt.gz $temp/$dbname.txt.gz
	} else {
		puts "Skipping download $dbname.txt.gz (already there)"
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
		puts "Downloading $dbname.sql ....."
		wgetfile ftp://hgdownload.cse.ucsc.edu/goldenPath/$build/database/$sqlfile $temp/$sqlfile
	} else {
		puts "Skipping download $dbname.sql (already there)"
	}
	if {![file exists $temp/$sqlfile]} {
		error $errmsg
	}
	# header
	# ------
	puts "Finding header...."
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
		puts "Gunzipping $dbname.txt.gz ...."
		if {[catch "exec zcat -f $temp/$dbname.txt.gz >> $temp/u_$dbname.tsv" errmsg]} {
			puts "Gunzip failed - $errmsg"
		}
	} else {
		foreach file $files {
			puts "Gunzipping $file ...."
			if {[catch "exec zcat -f $file >> $temp/u_$dbname.tsv" errmsg]} {
				puts "Gunzip failed - $errmsg"
			}
		}
	}
	puts "Sorting $filename ...."
	set fields [list_common {chrom chromosome begin start end stop name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts} $header]
	lappend fields {*}[list_lremove $header $fields]
	cg select -s - -f $fields $temp/u_$dbname.tsv $temp/su_$dbname.tsv
	file rename -force $temp/su_$dbname.tsv $filename
	puts "----------------------------------------------------"
}

proc downloaddbinfo {path build dbname} {
	set temp $path/tmp/$build
	file mkdir $temp
	# do documentation
	# ----------------
	if {![file exists $temp/$dbname.html]} {
		puts "Downloading $dbname.html ....."
		wgetfile http://genome.ucsc.edu/cgi-bin/hgTrackUi?db=$build&g=$dbname $temp/$dbname.html
	} else {
		puts "Skipping download $dbname.html (already there)"
	}
	puts "Making reg_${build}_$dbname.info"
	set c [file_read $temp/$dbname.html]
	if {[regexp {HGERROR-START} $c]} {
		puts "Skipping reg_${build}_$dbname.info: Could not download info"
		return
	}
	set title $dbname
	regexp {<TITLE>([^<]+)</TITLE>} $c t title
	regexp {<B style='font-family:serif; font-size:200%;'>([^<]+)</B>} $c t title
#	regsub {^.*<H2>Description</H2>} $c {} c
	if {![regsub {^.*<A NAME='TRACK_HTML'></A>} $c {} c]} {
		set c {}
	} else {
		regsub {</td><td nowrap>.*$} $c {} c
		regsub {</H2>\n+$} $c {</H2>} c
		set c [string_change $c [list <P> {} </P> {} <H2> "== " </H2> " =="]]
	}
	file_write $path/$build/reg_${build}_$dbname.info "= $dbname: [string trim $title] =\n\n[string trim $c] \n\n== Category ==\nAnnotation\n"
}

proc cg_downloaddbinfo {args} {
	if {([llength $args] < 2)} {
		error "format is: $::base resultdir build database ?...?\n - downloads databases from ucsc, 1000g, ... and converts to useful format"
	}
	foreach {path build dbname} $args break
	set dbnames [lrange $args 2 end]
	puts "----------------------------------------------------"
	file mkdir $path
	foreach dbname $dbnames {
			downloaddbinfo $path $build $dbname
	}
}

proc cg_downloaddb {args} {
	if {([llength $args] < 2)} {
		error "format is: $::base resultdir build database ?...?\n - downloads databases from ucsc, 1000g, ... and converts to useful format"
	}
	foreach {path build dbname} $args break
	set dbnames [lrange $args 2 end]
	puts "----------------------------------------------------"
	file mkdir $path
	foreach dbname $dbnames {
		if {$dbname eq "1000g"} {
			cg_download_1000g $path $build
		} elseif {$dbname eq "1000glow"} {
			cg_download_1000glow $path $build
		} elseif {$dbname eq "1000g3"} {
			set resultfile $path/$build/var_${build}_1000g3.tsv
			cg_download_1000g3 $resultfile $build
		} elseif {[regexp {snp.*} $dbname]} {
			set resultfile $path/$build/var_${build}_$dbname.tsv
			cg_download_dbsnp $resultfile $build $dbname
		} elseif {$dbname eq "evs"} {
			set resultfile $path/$build/var_${build}_evs.tsv
			cg_download_evs $resultfile $build [lindex $args 3]
		} elseif {$dbname eq "exac"} {
			set resultfile $path/$build/var_${build}_exac.tsv
			cg_download_exac $resultfile $build [lindex $args 3]
		} elseif {$dbname eq "phenotype"} {
			set resultfile $path/${build}/geneannot_${build}_phenotype.tsv
			cg_download_phenotype $resultfile $build [lindex $args 3]
		} else {
			downloaddb $path $build $dbname
			downloaddbinfo $path $build $dbname
		}
	}
}
