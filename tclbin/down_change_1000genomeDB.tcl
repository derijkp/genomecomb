#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

package require Extral


##############################################################
# 
#  This is a script that downloads a number of databases\
#  from 1000 genomes and makes some necessary changes
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

if {![file isdirectory $path/Databases/s_$build]} {
	file mkdir $path/Databases/s_$build
}


#downloading with wget

set gen_1000 { 	{CEU} \
	   	{CHB+JPT} \
		{YRI} \
}

set i 0

while {$i < [llength $gen_1000]} {
	if {![file isfile $path/Databases/s_$build/[lindex $gen_1000 $i].tsv]} {

		if {![file isfile $path/Databases/s_$build/[lindex $gen_1000 $i].frq]} {	

			if {![file isfile $path/Databases/s_$build/[lindex $gen_1000 $i].SRP000031.2010_03.genotypes.vcf]} {
			
				puts "----------------------------------------------------"
				puts "Downloading [lindex $gen_1000 $i].SRP000031.2010_03.genotypes.vcf....."
				puts "This may take a couple minutes...."
				if {[catch "exec wget --tries=45 --directory-prefix=$path/Databases/s_$build/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/release/2010_03/pilot1/[lindex $gen_1000 $i].SRP000031.2010_03.genotypes.vcf.gz " errmsg]} {
					puts "Downloading database succeeded"
				}	

				puts "Gunzipping [lindex $gen_1000 $i].SRP000031.2010_03.genotypes.vcf...."			
				puts "This may take a few minutes...."
				if {[catch "exec gunzip -f $path/Databases/s_$build/[lindex $gen_1000 $i].SRP000031.2010_03.genotypes.vcf.gz " errmsg]} {
					puts "Gunzip failed - $errmsg"
				}
				puts "----------------------------------------------------"

			} else {
				puts "----------------------------------------------------"
				puts "The file '[lindex $gen_1000 $i].SRP000031.2010_03.genotypes.vcf.gz' already exists..."
				puts "Skipping the download..."
				puts "----------------------------------------------------"
			}

			puts "----------------------------------------------------"
			puts "Making the first necessary changes to [lindex $gen_1000 $i].SRP000031.2010_03.genotypes.vcf creating [lindex $gen_1000 $i].frq"
			puts "Get some coffee, this may take an half hour..."
			puts "Log file will be $path/Databases/s_$build/[lindex $gen_1000 $i].log "
			if {[catch "exec vcftools1.2 --vcf $path/Databases/s_$build/[lindex $gen_1000 $i].SRP000031.2010_03.genotypes.vcf --freq --out $path/Databases/s_$build/[lindex $gen_1000 $i]" errmsg]} {
				puts "Change with vcftools failed - $errmsg"
			}
		
			puts "----------------------------------------------------"

		} else {
			puts "----------------------------------------------------"
			puts "The file '[lindex $gen_1000 $i].frq' already exists..."
			puts "Skipping the first changes with vcftools..."
			puts "----------------------------------------------------"
		}
		
		puts "----------------------------------------------------"
		puts "Making some other changes..... "
		puts "including going from hg19 to hg18 with liftOver..."
		puts "From [lindex $gen_1000 $i].frq to [lindex $gen_1000 $i].tsv"
		

		if {[catch "open $path/Databases/s_$build/[lindex $gen_1000 $i].frq r" fileid]} {
			puts "Could not open frequency file - $fileid"
			exit 1
		}
		set temp1 [tempfile get]
		if {[catch "open $temp1 w" fileid_out]} {
			puts "Could not open outputfile - $fileid_out"
			exit 1
		}	

		set in [gets $fileid]
		set headers [split $in "\t"]

		set column 0
		foreach head $headers {
			switch -glob $head {
				CHROM {set chrom $column}
				POS {set begin $column}
				N_ALLELES {set end $column}
				\{ALLELE:FREQ\} {set allel_freq $column}
			}
			incr column
		}

		set new_headers [lreplace $headers $chrom $chrom "#chrom"]
		set new_headers [lreplace $new_headers $begin $begin "chromStart"]
		set new_headers [lreplace $new_headers $end $end "chromEnd"]
		set new_headers [lreplace $new_headers $allel_freq $allel_freq "allel_freq"]

		puts $fileid_out [join $new_headers \t]

		set j 0
		set in [gets $fileid]
		while {![eof $fileid]} {
				incr j
				set line [split $in "\t"]

				# Print every 100000 lines
				if {[expr {$j%100000}]==0} {	
					puts "reading line $j"
				}

				set new_line [lreplace $line $end $end [lindex $line $begin]]
				set new_line [lreplace $new_line $chrom $chrom "chr[lindex $line $chrom]"]
				set new_line [lreplace $new_line $begin $begin [expr [lindex $line $begin] - 1]]

				# Putting the different frequencies in 1 column, depending on N_ALLELES
				if {[lindex $line $end] == 1} {
					puts $fileid_out [join $new_line \t]
					set in [gets $fileid]
				} elseif {[lindex $line $end] == 2} {
					set new_line [lreplace $new_line $allel_freq $allel_freq "[lindex $line $allel_freq],[lindex $line [expr $allel_freq + 1]]"]
					set new_line [lreplace $new_line [expr $allel_freq + 1] [expr $allel_freq + 1] ]
					puts $fileid_out [join $new_line \t]
					set in [gets $fileid]
				} elseif {[lindex $line $end] == 3} {
					set new_line [lreplace $new_line $allel_freq $allel_freq "[lindex $line $allel_freq],[lindex $line [expr $allel_freq + 1]],[lindex $line [expr $allel_freq + 2]]"]
					set new_line [lreplace $new_line [expr $allel_freq + 1] [expr $allel_freq + 2] ]
					puts $fileid_out [join $new_line \t]
					set in [gets $fileid]
				} elseif {[lindex $line $end] == 4} {
					set new_line [lreplace $new_line $allel_freq $allel_freq "[lindex $line $allel_freq],[lindex $line [expr $allel_freq + 1]],[lindex $line [expr $allel_freq + 2]],[lindex $line [expr $allel_freq + 3]]"]
					set new_line [lreplace $new_line [expr $allel_freq + 1] [expr $allel_freq + 3] ]
					puts $fileid_out [join $new_line \t]
					set in [gets $fileid]

				}

		}
		close $fileid
		close $fileid_out

		#using liftOver
		set temp2 [tempfile get]
		if {[catch "exec liftOver $temp1 $path/Liftover/hg19ToHg18.over.chain $temp2 unmapped.bed" errmsg]} {
			puts "$errmsg"
		}
		if {[catch "open $path/Databases/s_$build/[lindex $gen_1000 $i].tsv w" fileid_out]} {
			puts "Could not open outputfile - $fileid_out"
			exit 1
		}
		set new_headers [lreplace $new_headers $chrom $chrom "chrom"]
		puts $fileid_out [join $new_headers \t]
		close $fileid_out
		if {[catch "exec cat $temp2 >> $path/Databases/s_$build/[lindex $gen_1000 $i].tsv" errmsg]} {
			puts "Pasting head to file failed - $errmsg"
		}
	
		file delete -force $path/Databases/s_$build/[lindex $gen_1000 $i].SRP000031.2010_03.genotypes.vcf $path/Databases/s_$build/[lindex $gen_1000 $i].frq unmapped.bed $path/Databases/s_$build/[lindex $gen_1000 $i].log
		puts "----------------------------------------------------"

	} else {
		puts "----------------------------------------------------"
		puts "The file '[lindex $gen_1000 $i].tsv' already exists..."
		puts "Skipping the changes...."
		puts "----------------------------------------------------"
	}
	incr i
}





tempfile clean
return 0





