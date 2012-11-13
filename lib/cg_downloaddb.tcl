#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

package require Extral

##############################################################
# 
#  This is a script that downloads a number of databases
#  
#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# parts by Annelies Cassiers (VIB - University of Antwerp)
#
##############################################################

proc downloaddb {path build dbname} {
	file mkdir $path
	set filename $path/$build/ucsc_${build}_$dbname.tsv
	set temp $path/tmp/$build
	file mkdir $temp
	if {[file isfile $filename]} {
		puts "The file '$filename' already exists..."
		puts "Skipping the download..."
		puts "----------------------------------------------------"
		return
	}
	catch {file mkdir -force [file dir $filename]}
	set single 1
	set chromosomes {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 M X Y}
	if {![file exists $temp/$dbname.txt.gz]} {
		puts "Downloading $dbname.txt.gz ....."
		wgetfile ftp://hgdownload.cse.ucsc.edu/goldenPath/$build/database/$dbname.txt.gz $temp/$dbname.txt.gz
	} else {
		puts "Skipping download $dbname.txt.gz (already there)"
	}
	if {![file exists $temp/$dbname.txt.gz]} {
		set single 0
		foreach chr $chromosomes {
			if {![file exists $temp/chr${chr}_$dbname.txt.gz]} {
				puts "Downloading chr${chr}_$dbname.txt.gz ....."
				wgetfile ftp://hgdownload.cse.ucsc.edu/goldenPath/$build/database/chr${chr}_$dbname.txt.gz $temp/chr${chr}_$dbname.txt.gz
				if {![file exists $temp/chr${chr}_$dbname.txt.gz]} {
					if {$chr eq "1"} {
						puts "Could not download $dbname.txt.gz"
					} else {
						puts "Could not download chr${chr}_$dbname.txt.gz"
					}
					exit 1
				}
			}
		}
	}
	if {$single} {
		set sqlfile $dbname.sql
	} else {
		set sqlfile chr1_$dbname.sql
	}
	if {![file exists $temp/$sqlfile]} {
		puts "Downloading $dbname.sql ....."
		catch {exec wget --directory-prefix=$temp ftp://hgdownload.cse.ucsc.edu/goldenPath/$build/database/$sqlfile} errmsg
	} else {
		puts "Skipping download $dbname.sql (already there)"
	}
	if {![file exists $temp/$sqlfile]} {
		puts $errmsg
		exit 1
	}
	# header
	# ------
	puts "Finding header...."
	if {[catch "open $temp/$sqlfile r" fileid]} {
		puts "Could not open sql file - $fileid"
		exit 1
	}
	set line [gets $fileid]
	set header ""
	while {![eof $fileid]} {
		if {[lsearch $line "CREATE"] != -1} {
			set line [gets $fileid]
			set head [string trim [lindex $line 0] (`) ]
			while {![inlist {KEY PRIMARY UNIQUE )} $head]} {	
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
		foreach chr $chromosomes {
			puts "Gunzipping chr${chr}_$dbname.txt.gz ...."
			if {[catch "exec zcat -f $temp/chr${chr}_$dbname.txt.gz >> $temp/u_$dbname.tsv" errmsg]} {
				puts "Gunzip failed - $errmsg"
			}
		}
	}
	puts "Sorting $filename ...."
	if {[catch {tsv_basicfields $fields 3}]} {
		file copy $temp/u_$dbname.tsv $filename
		# file delete -force $temp/$dbname.txt $temp/$dbname.sql $temp/u_$dbname.tsv
	} else {
		set fields [list_common {chrom start end name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts} $header]
		lappend fields {*}[list_lremove $header $fields]
		cg select -q {$chrom ~ /^chr[0-9XYM][0-9]*$/} -s {chrom start end} -f $fields $temp/u_$dbname.tsv $filename
		# file delete -force $temp/$dbname.txt $temp/$dbname.sql $temp/u_$dbname.tsv
	}
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
		puts stderr "format is: $::base resultdir build database ?...?"
		puts stderr " - downloads databases from ucsc, 1000g, ... and converts to useful format"
		exit 1
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
		puts stderr "format is: $::base resultdir build database ?...?"
		puts stderr " - downloads databases from ucsc, 1000g, ... and converts to useful format"
		exit 1
	}
	foreach {path build dbname} $args break
	set dbnames [lrange $args 2 end]
	puts "----------------------------------------------------"
	file mkdir $path
	foreach dbname $dbnames {
		if {$dbname eq "1000g"} {
			downloaddb_1000g $path $build
		} elseif {$dbname eq "1000glow"} {
			downloaddb_1000glow $path $build
		} elseif {[regexp {snp.*} $dbname]} {
			downloaddb_dbsnp $path $build $dbname
		} else {
			downloaddb $path $build $dbname
			downloaddbinfo $path $build $dbname
		}
	}
}

if 0 {
	set path /media/663d83bb-851c-4dbb-8c03-e8815d28e483/refseq
	set nuild hg18
}

proc cg_calcsequencedgenome {args} {
	if {([llength $args] != 2)} {
		puts stderr "format is: $::base resultdir build"
		puts stderr " - calculate sequenced genome from genome ifas (regions without Ns)"
		exit 1
	}
	foreach {path build} $args break
	set file $path/$build/genome_$build.ifas
	set f [open $path/$build/genome_$build.ifas]
	file mkdir $path/$build/extra/
	set o [open $path/$build/extra/reg_${build}_sequencedgenome.tsv w]
	puts $o chromosome\tbegin\tend\tsize
	while {![eof $f]} {
		set name [gets $f]
		if {$name eq ""} continue
		if {![regexp {chromosome ([0-9A-Z]+)} $name temp chr]} {
			if {![regexp {chr([0-9A-Z]+)} $name temp chr]} {
				error "no chromosome found in line $name"
			}
		}
		putslog $name\n$chr
		set seq [gets $f]
		set indices [regexp -all -inline -indices {[^N]{1,}} $seq]
		putslog Writing
		list_foreach {pbegin pend} $indices break
		list_foreach {begin end} $indices {
			incr end
			if {$begin > [expr {$pend+1}]} {
				puts $o chr$chr\t$pbegin\t$pend
				set pbegin $begin
			}
			set pend $end
		}
		puts $o chr$chr\t$pbegin\t$pend
	}
	close $o
	close $f
}

if {[info exists argv0] && [file tail [info script]] eq [file tail $argv0]} {
	package require pkgtools
	set appdir [file dir [pkgtools::startdir]]
	lappend auto_path $appdir/lib
	append env(PATH) :[file dir [file dir $appdir]]/bin:$appdir/bin
	package require Extral
	set ::base $scriptname
	cg_downloaddb {*}$argv
}


