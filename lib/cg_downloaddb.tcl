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
	catch {file mkdir [file dir $filename]}
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
		foreach chr $chromosomes {
			puts "Gunzipping chr${chr}_$dbname.txt.gz ...."
			if {[catch "exec zcat -f $temp/chr${chr}_$dbname.txt.gz >> $temp/u_$dbname.tsv" errmsg]} {
				puts "Gunzip failed - $errmsg"
			}
		}
	}
	puts "Sorting $filename ...."
	set fields [list_common {chrom chromosome begin start end stop name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts} $header]
	lappend fields {*}[list_lremove $header $fields]
	cg select -s - -f $fields $temp/u_$dbname.tsv $filename.temp
	file rename -force $filename.temp $filename
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

proc downloaddb_evs {path build {url {}}} {
	if {$url eq ""} {
		set url http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.protein-hgvs-update.snps_indels.vcf.tar.gz
	}
	file mkdir ${path}/tmp/$build/evs
	cd ${path}/tmp/$build/evs
	exec -ignorestderr wget -c --tries=45 --directory-prefix=${path}/tmp/$build/evs $url
	exec tar xvzf [file tail $url]
	foreach file [glob *.vcf] {
		cg vcf2tsv $file [file root $file].tsv
	}
	set files [lsort -dict [glob ESP6500SI-V2-SSA137.updatedProteinHgvs.*.snps_indels.tsv]]
	if {[llength $files] != 24} {error "not enough files found"}
	cg cat {*}$files > var_${build}_evs.tsv.temp
#	cg select -f {chromosome begin end type 
#		{ea_freqp=lrange(vformat("%.3f",(100.0 @* $EA_AC) @/ lsum($EA_AC)),0,"end-1")} 
#		{aa_freqp=lrange(vformat("%.3f",(100.0 @* $AA_AC) @/ lsum($AA_AC)),0,"end-1")}
#		*
#	} var_${build}_evs.tsv.temp var_${build}_evs.tsv.temp2
	set f [open var_${build}_evs.tsv.temp]
	set header [tsv_open $f]
	set bposs [tsv_basicfields $header]
	set o [open var_${build}_evs.tsv.temp2 w]
	set nheader [list {*}[list_sub $header $bposs] ea_freqp aa_freqp ea_mfreqp aa_mfreqp {*}[list_sub $header -exclude $bposs]]
	puts $o [join $nheader \t]
	set poss [list_cor $header {alt EA_AC AA_AC GTS EA_GTC AA_GTC}]
	while {![eof $f]} {
		set line [split [gets $f] \t]
		if {![llength $line]} continue
		foreach {alt EA_AC AA_AC GTS EA_GTC AA_GTC} [list_sub $line $poss] break
		set EA_AC [split $EA_AC ,]
		set AA_AC [split $AA_AC ,]
		set GTS [split $GTS ,]
		set EA_GTC [split $EA_GTC ,]
		set AA_GTC [split $AA_GTC ,]
		set sum [lmath_sum $EA_AC]
		set ea_freqp {}
		foreach el [lrange $EA_AC 0 end-1] {
			lappend ea_freqp [trimformat %.3f [expr {(100.0*$el)/$sum}]]
		}
		set sum [lmath_sum $AA_AC]
		set aa_freqp {}
		foreach el [lrange $AA_AC 0 end-1] {
			lappend aa_freqp [trimformat %.3f [expr {(100.0*$el)/$sum}]]
		}
		set alt [split $alt ,]
		set apos 1
		set ea_mfreqp {}
		set aa_mfreqp {}
		set easum [lmath_sum $EA_GTC]
		set aasum [lmath_sum $AA_GTC]
		foreach a $alt {
			set pos [lsearch $GTS ${a}${a}]
			if {$pos == -1} {
				set pos [lsearch $GTS A${apos}A${apos}]
			}
			if {$pos == -1} {
				lappend ea_mfreqp .
				lappend aa_mfreqp .
			} else {
				lappend ea_mfreqp [trimformat %.3f [expr {100.0*[lindex $EA_GTC $pos]/$easum}]]
				lappend aa_mfreqp [trimformat %.3f [expr {100.0*[lindex $AA_GTC $pos]/$aasum}]]
			}
			incr apos
		}
		set result [list_sub $line $bposs]
		lappend result [join $ea_freqp ,] [join $aa_freqp ,] [join $ea_mfreqp ,] [join $aa_mfreqp ,]
		lappend result {*}[list_sub $line -exclude $bposs]
		puts $o [join $result \t]
	}
	close $o
	close $f
	file_write $path/$build/extra/var_${build}_evs.tsv.opt "fields\t{ea_freqp aa_freqp ea_mfreqp aa_mfreqp}\n"
file_write $path/$build/extra/var_${build}_evs.tsv.info [string trim {
Exome variant server variants

Data Usage
----------
We request that any use of data obtained from the NHLBI GO ESP Exome Variant Server be cited in publications.

Citation
--------
Exome Variant Server, NHLBI GO Exome Sequencing Project (ESP), Seattle, WA (URL: http://evs.gs.washington.edu/EVS/) [date (month, yr) accessed].

Annotation fields added
-----------------------
ea_freqp: frequency (as percentage!) of alt aleles in european american population
aa_freqp: frequency (as percentage!) of alt aleles in african american population
ea_mfreqp: frequency (as percentage!) of homozygous alt alele genotype in european american population
aa_mfreqp: frequency (as percentage!) of homozygous alt alele genotype in african american population

more information on http://evs.gs.washington.edu/EVS/
}]
	file rename var_${build}_evs.tsv.temp2 $path/$build/extra/var_${build}_evs.tsv
}

proc cg_downloaddb {args} {
	if {([llength $args] < 2)} {
		puts stderr "format is: $::base resultdir build database ?...?"
		puts stderr " - downloads databases from ucsc, 1000g, ... and converts to useful format"
		exit 1
	}
	foreach {path build dbname} $args break
	if {$dbname eq "evs"} {
		downloaddb_evs $path $build [lindex $args 3]
		return
	} elseif {$dbname eq "exac"} {
		downloaddb_exac $path $build [lindex $args 3]
		return
	}
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
	set build hg18
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
	set o [open $path/$build/extra/reg_${build}_sequencedgenome.tsv.temp w]
	puts $o chromosome\tbegin\tend\tsize
	while {![eof $f]} {
		set name [gets $f]
		if {$name eq ""} continue
		if {![regexp {chromosome ([0-9A-Z]+)} $name temp chr]} {
			if {![regexp {chr([0-9A-Z]+)} $name temp chr]} {
				set chr $name
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
	file rename $path/$build/extra/reg_${build}_sequencedgenome.tsv.temp $path/$build/extra/reg_${build}_sequencedgenome.tsv
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


