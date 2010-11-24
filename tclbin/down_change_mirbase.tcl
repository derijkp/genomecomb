#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

package require Extral


##############################################################
# 
#  This is a script that downloads mirbase DB
#  and makes some necessary changes
#  
#  Written by Annelies Cassiers
#
##############################################################

if {[llength $argv] < 2} {
	puts "Format is: directory (eg '/home/user'), build_version (eg 'hg18') "
	exit 1
}

set path [lindex $argv 0]
set build [lindex $argv 1]

if {![file isdirectory $path/Databases/$build]} {
	file mkdir $path/Databases/$build
}


#downloading with wget

if {![file isfile $path/Databases/$build/hsa.tsv]} {
	if {![file isfile $path/Databases/$build/hsa.bed]} {
		if {![file isfile $path/Databases/$build/hsa.gff]} {

			puts "----------------------------------------------------"
			puts "Downloading hsa.gff....."
			puts "This may take a few seconds...."
			if {[catch "exec wget --tries=45 --directory-prefix=$path/Databases/$build/ ftp://mirbase.org/pub/mirbase/16/genomes/hsa.gff" errmsg]} {
				puts "Downloading database succeeded"
			}	
			puts "----------------------------------------------------"

		} else {
			puts "----------------------------------------------------"
			puts "The file 'hsa.gff' already exists..."
			puts "Skipping the download..."
			puts "----------------------------------------------------"
		}

		# Van gff formaat naar bed zodat het in liftOver kan en in Aregio
	
		puts "----------------------------------------------------"
		puts "Making some changes..... "
		puts "From hsa.gff to hsa.bed...."
	

		if {[catch "open $path/Databases/$build/hsa.gff r" fileid]} {
			puts "Could not open frequency file - $fileid"
			exit 1
		}

		if {[catch "open $path/Databases/$build/hsa.bed w" fileid_out]} {
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
		puts "----------------------------------------------------"

	} else { 
		puts "----------------------------------------------------"
		puts "The file 'hsa.bed' already exists..."
		puts "Skipping the changes..."
		puts "----------------------------------------------------"
	}
	puts "----------------------------------------------------"
	puts "Use liftOver to make hg18 from hg19...."
	
	set temp1 [tempfile get]
	if {[catch "exec liftOver $path/Databases/$build/hsa.bed $path/Liftover/hg19ToHg18.over.chain $temp1 unmapped.bed" errmsg]} {
		puts "$errmsg"
	}	
	
	# temp file maken met header, deze dan tegen de tsv file kleven...
	
	if {[catch "open $path/Databases/$build/hsa.tsv w" fileid_out]} {
		puts "Could not open tsv file - $fileid_out"
		exit 1
	}
		 
	set new_header {chrom chromStart chromEnd name strand}
	puts $fileid_out [join $new_header \t]
	close $fileid_out
	if {[catch "exec cat $temp1 >> $path/Databases/$build/hsa.tsv" errmsg]} {
		puts "Pasting head to file failed - $errmsg"
	}	
	file delete -force $path/Databases/$build/hsa.bed $path/Databases/$build/hsa.gff unmapped.bed
	puts "----------------------------------------------------"

} else {
	puts "----------------------------------------------------"
	puts "The file 'hsa.tsv' already exists..."
	puts "Skipping the changes..."
	puts "----------------------------------------------------"
}

tempfile clean
return 0
