#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

package require Extral


##############################################################
# 
#  This is a script that downloads a number of databases
#  
#  Written by Annelies Cassiers
#  adapted and extended by Peter De Rijk
#
##############################################################

if 0 {

	package require Tclx
	signal -restart error SIGINT
	lappend auto_path /home/peter/dev/completegenomics/lib
	package require Extral
	cd /complgen/refseq/hg18
	set path /complgen/refseq/hg18
	set build hg18
	set dbname tRNAs
	set dbname simpleRepeat
	set dbname rmsk
	set dbname cytoBand

cg downloaddb /complgen/refseq/hg18 hg18 simpleRepeat microsat rmsk genomicSuperDups chainSelf
cg downloaddb /complgen/refseq/hg18 hg18 omimGene phastConsElements44way oreganno rnaGene tRNAs tfbsConsSites targetScanS evofold
cg downloaddb /complgen/refseq/hg18 hg18 cytoBand dgv gwasCatalog phastConsElements28wayPlacMammal phastConsElements28way snp130 wgRna
cg downloaddb /complgen/refseq/hg18 hg18 mirbase

}

proc downloaddb_mirbase {path build} {
	set filename $path/reg_${build}_mirbase.tsv
	set temp $path/tmp/$build
	if {[file isfile $filename]} {
		puts "The file '$filename' already exists..."
		puts "Skipping the changes..."
		puts "----------------------------------------------------"
		return
	}
	if {![file isfile $temp/mirbase.bed]} {
		if {![file isfile $temp/hsa.gff]} {
			puts "Downloading mirbase.gff....."
			catch {exec wget --tries=45 --directory-prefix=$temp/ ftp://mirbase.org/pub/mirbase/16/genomes/hsa.gff} errmsg
			if {![file exists $temp/hsa.gff]} {
				puts $errmsg
				exit 1
			}
		} else {
			puts "The file 'hsa.gff' already exists..."
			puts "Skipping the download..."
		}
		# Van gff formaat naar bed zodat het in liftOver kan en in Aregio
		puts "Making some changes..... "
		puts "From hsa.gff to mirbase.bed...."
		if {[catch {open $temp/hsa.gff r} fileid]} {
			puts "Could not open frequency file - $fileid"
			exit 1
		}
		if {[catch {open $temp/mirbase.bed w} fileid_out]} {
			puts "Could not open outputfile - $fileid_out"
			exit 1
		}	
		set in [gets $fileid]
		while {![eof $fileid]} {
			# commend lines overslaan
			if {[string match #* $in]} {
				set in [gets $fileid]
			} else {
				set line [split $in "\t"]
				set name [lindex $line 8]
				set new_name [string map {"; " ","} $name]
				set new_name [string map {";" ""} $new_name]
				set new_line "chr[lindex $line 0] [lindex $line 3] [lindex $line 4] $new_name [lindex $line 6]"
				puts $fileid_out [join $new_line \t]
				set in [gets $fileid]
			}
		}
		close $fileid
		close $fileid_out
	} else { 
		puts "The file 'mirbase.bed' already exists..."
		puts "Skipping the changes..."
	}
	puts "Use liftOver to make hg18 from hg19...."
	set dir [file dir [exec which liftOver]]
	if {[catch {exec liftOver $temp/mirbase.bed $dir/hg19ToHg18.over.chain $temp/hg18_mirbase.bed $temp/unmapped.bed} errmsg]} {
		puts "$errmsg"
	}	
	# temp file maken met header, deze dan tegen de tsv file kleven...
	if {[catch {open $temp/mirbase.tsv w} fileid_out]} {
		puts "Could not open tsv file - $fileid_out"
		exit 1
	}
	set new_header {chrom start end name strand}
	puts $fileid_out [join $new_header \t]
	close $fileid_out
	if {[catch "exec cat $temp/hg18_mirbase.bed >> $temp/mirbase.tsv" errmsg]} {
		puts "Pasting head to file failed - $errmsg"
	}	
	puts "Sorting $filename ...."
	cg select -q {$chrom ~ /^chr[0-9XYM][0-9]*$/} -s {chrom start end} $temp/mirbase.tsv $filename
	# file delete -force $temp/mirbase.bed $temp/mirbase.gff unmapped.bed
	puts "----------------------------------------------------"
}

proc downloaddb {path build dbname} {
	file mkdir $path
	set filename $path/ucsc_${build}_$dbname.tsv
	set temp $path/tmp/$build
	file mkdir $temp
	if {[file isfile $filename]} {
		puts "The file '$path/ucsc_${build}_$dbname.tsv' already exists..."
		puts "Skipping the download..."
		puts "----------------------------------------------------"
		return
	}
	set single 1
	set chromosomes {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 M X Y}
	if {![file exists $temp/$dbname.txt.gz]} {
		puts "Downloading $dbname.txt.gz ....."
		catch {exec wget --tries=45 --directory-prefix=$temp ftp://hgdownload.cse.ucsc.edu/goldenPath/$build/database/$dbname.txt.gz} errmsg
	} else {
		puts "Skipping download $dbname.txt.gz (already there)"
	}
	if {![file exists $temp/$dbname.txt.gz]} {
		set single 0
		foreach chr $chromosomes {
			if {![file exists $temp/chr${chr}_$dbname.txt.gz]} {
				puts "Downloading chr${chr}_$dbname.txt.gz ....."
				set e [catch {exec wget --tries=45 --directory-prefix=$temp ftp://hgdownload.cse.ucsc.edu/goldenPath/$build/database/chr${chr}_$dbname.txt.gz} errmsg]
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
	set fields [list_common {chrom start end name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts} $header]
	lappend fields {*}[list_lremove $header $fields]
	cg select -q {$chrom ~ /^chr[0-9XYM][0-9]*$/} -s {chrom start end} -f $fields $temp/u_$dbname.tsv $filename
	# file delete -force $temp/$dbname.txt $temp/$dbname.sql $temp/u_$dbname.tsv
	puts "----------------------------------------------------"
}

proc cg_downloaddb {args} {
	if {([llength $args] < 2)} {
		puts stderr "format is: $::base resultdir build database ?...?"
		puts stderr " - downloads databases from ucsc and converts to region format"
		exit 1
	}
	foreach {path build dbname} $args break
	set dbnames [lrange $args 2 end]
	puts "----------------------------------------------------"
	foreach dbname $dbnames {
		if {$dbname eq "mirbase"} {
			downloaddb_mirbase $path $build
		} elseif {$dbname eq "1000g"} {
			downloaddb_1000g $path $build
		} elseif {[regexp {snp.*} $dbname]} {
			downloaddb_dbsnp $path $build $dbname
		} else {
			downloaddb $path $build $dbname
		}
	}
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


