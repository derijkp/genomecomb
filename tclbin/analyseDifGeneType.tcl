#!/bin/sh
# the next line restarts using tclsh \
exec tclsh "$0" "$@"

package require Extral

##############################################################
# 
#  This is a script that analyses the difference between
#  the type called by knownGene and ensGene
#
#  Written by Annelies Cassiers
#
##############################################################

if {[llength $argv] < 1} {
	puts "Format is: filename to be analysed "
	exit 1
}


if {[catch "exec cg select -h [lindex $argv 0]" headers]} {
	puts "could not execute cg select - $headers"
	puts 1
}

set column 1
foreach head $headers {
	switch -glob -nocase $head {
		knownGene_type {set UCSC $column}
		ensGene_type {set ens $column}
		avsift_minimalScore {set sift $column}
		phastConsElements44way_name {set phastCon $column}
		oreganno_name {set oreganno $column}
		xRef {set dbsnp $column}
		exonCategory {set ex $column}
	}
	incr column
}

puts "preparing inputfile.."
puts "this may take a few minutes.."

set temp1 [tempfile get] 
# ,${ex}
if {[catch "exec cut -f${UCSC},${ens},${sift},${phastCon},${oreganno},${dbsnp} [lindex $argv 0] > $temp1 " errmsg]} {
	puts "Something went wrong while making new file - $errmsg"
	exit 1
}

if {[catch "open $temp1 r" fileid]} {
	puts "Could not open new file - $fileid"
	exit 1
}

set column 0
set in [gets $fileid]
set headers [split $in "\t"]
foreach head $headers {
	switch -glob -nocase $head {
		knownGene_type {set UCSC $column}
		ensGene_type {set ens $column}
		avsift_minimalScore {set sift $column}
		phastConsElements44way_name {set phastCon $column}
		oreganno_name {set oreganno $column}
		xRef {set dbsnp $column}
		exonCategory {set ex $column}
	}
	incr column
}



# read real first line 
set in [gets $fileid]
set uni 0
set alone 0
set empty 0
set exon 0
set splice 0
set intron 0
set UTR5 0
set UTR3 0
set ncRNA 0
set upstream 0
set downstream 0
set intergen 0
set more 0
set sift_empty 0
set sift_under 0
set sift_up 0
set phast_empty 0
set phast_full 0
set or_empty 0
set or_full 0
set db_empty 0
set db_full 0
set ex_utr 0
set ex_ex 0

set i 1
while {![eof $fileid]} {
	# Print every 500000 lines
	if {[expr {$i%500000}]==0} {	
		puts "reading line $i"
	}
	set line [split $in "\t"]

	########
	# venn
	########
	
	if {[lindex $line $UCSC] == " "} {
		incr empty
	} elseif {[lindex $line $UCSC] == [lindex $line $ens]} {
		incr uni
	} else {
		incr alone
	}

	#######
	# bar
	#######
	#if you want to know were the variants lie in ensGene replace the 'UCSC' with 'ens'

	if {[lindex $line $UCSC] == " "} {
		#do nothing
	} elseif {[lindex $line $UCSC] == "exonic"} {
		incr exon
	} elseif {[lindex $line $UCSC] == "splicing"} {
		incr splice
	} elseif {[lindex $line $UCSC] == "intronic"} {
		incr intron
	} elseif {[lindex $line $UCSC] == "UTR5"} {
		incr UTR5
	} elseif {[lindex $line $UCSC] == "UTR3"} {
		incr UTR3
	} elseif {[lindex $line $UCSC] == "ncRNA"} {
		incr ncRNA
	} elseif {[lindex $line $UCSC] == "upstream"} {
		incr upstream
	} elseif {[lindex $line $UCSC] == "downstream"} {
		incr downstream
	} elseif {[lindex $line $UCSC] == "intergenic"} {
		incr intergen
	} else {
		incr more
	} 
	
	#######
	# sift
	#######

	if {[lindex $line $sift] == " "} {
		incr sift_empty
	} elseif {[expr ([lindex $line $sift] <= 0.05)]} {
		incr sift_under
	} elseif {[expr ([lindex $line $sift] > 0.05)]} {
		incr sift_up
	}
	
	######
	# phastCon
	######
	if {[lindex $line $phastCon] ==	""} {
		incr phast_empty
	} else {
		incr phast_full
	}

	######
	# oreganno
	######
	if {[lindex $line $oreganno] ==	""} {
		incr or_empty
	} else {
		incr or_full
	}

	

	######
	# dbsnp
	######
	if {[lindex $line $dbsnp] ==	""} {
		incr db_empty
	} else {
		incr db_full
	}
	
	######
	# utr cg
	######
	#switch [lindex $line $ex] {
	#	EXON {incr ex_ex}
	#	UTR {incr ex_utr}
	#}


	set in [gets $fileid]
	incr i

}

close $fileid

##################
#
# venn diagram
#
##################


puts "There are $empty variants were there was no annotation."
puts "There are $uni variants with the same annotation in UCSC and ensemble genes."
puts "There are $alone variants with different annotation in UCSC and ensemle genes."
set A [expr ($uni.0 + $alone.0) / 1000000]
set B [expr ($uni.0 + $alone.0) / 1000000]
set C [expr $uni.0 / 1000000]
set D [expr 100 - ( $alone.0 / ($uni.0 + $alone.0) *100)]
puts "That means that $D % of the variants have the same annotation in both databases"
puts " Use the following url in your browser: "
puts " http://chart.apis.google.com/chart?cht=v&chd=t:$A,$B,0,$C&chs=250x300&chco=ADDE63,63C6DE&chtt=Overlap+gene+databases|$D+%&chdl=UCSC|ensemble"

#################
#
# bar diagram
#
#################

puts " "
puts "When following the UCSC annotation: "
puts "There are $exon variants located in an exon"
puts "There are $splice variants located in a splice site"
puts "There are $intron variants located in an intron"
puts "There are $UTR5 variants located in an UTR5 region"
puts "There are $UTR3 variants located in an UTR3 region"
puts "There are $ncRNA variants located in a ncRNA region"
puts "There are $upstream variants located in an upstream region"
puts "There are $downstream variants located in a downstream region"
puts "There are $intergen variants located between genes (intergenic)"
puts "There are $more variants located in more then 1 of the above regions (eg upstream;downstream)"

set other [expr $exon + $splice + $UTR5 + $UTR3 + $ncRNA + $upstream + $downstream + $more]
set max 0
if {$intron > $intergen} {
	set max $intron
} else {
	set max $intergen
}

if {$max < $other} {
	set max $other
}

puts " Use the following urls in your browser: "
puts " http://chart.apis.google.com/chart?cht=bvs&chd=t:$intron,$intergen,$other,$empty&chds=0,$max&chs=500x300&chf=b0,lg,0,FFE7C6,0,76A4FB,1&chtt=Position+of+the+variant|based+on+UCSC+annotation|[file tail [file rootname [lindex $argv 0]]]&chxt=x,y,y&chxl=0:|intron|intergenic|other|unknown|2:|Number+of+variants|&chxp=2,50&chbh=a,25&chma=20,20,20,30|80,20&chxr=0,0,100|1,0,$max&chm=N,000000,0,-1,10"
set list_max "$exon $splice $UTR5 $UTR3 $ncRNA $upstream $downstream $more"
set max [lindex [lsort -integer -decreasing $list_max ] 0]
puts " "
puts " http://chart.apis.google.com/chart?cht=bvs&chd=t:$exon,$splice,$UTR5,$UTR3,$ncRNA,$upstream,$downstream,$more&chds=0,$max&chs=700x300&chf=b0,lg,0,FFE7C6,0,76A4FB,1&chtt=Position+in+'other'+of+the+variant|based+on+UCSC+annotation|[file tail [file rootname [lindex $argv 0]]]&chxt=x,y,y&chxl=0:|exon|spliceSite|UTR5|UTR3|ncRNA|upstream|downstream|more|2:|Number+of+variants|&chxp=2,50&chxr=0,0,$max&chbh=a,25&chma=20,20,20,30|80,20&chxr=0,0,100|1,0,$max&chm=N,000000,0,-1,10"

##############
#
# sift
#
##############



set sift_all [expr $sift_up + $sift_under]
set pro [expr ($sift_under.0 / ($sift_empty.0 + $sift_all.0 ) *100)]
puts " "
puts "There are $sift_empty variants annotated without a sift score"
#puts "There are $sift_all variants annotated with a sift score"
puts "There are $sift_up variants predicted by sift score to be benign (>0.05)"
puts "There are $sift_under variants predicted to be malignant (<= 0.05)"
puts "This is $pro % of the variants that are predicted to be malignant"

#############
#
# phastCons
#
#############

set pro [expr  ($phast_full.0 / ($phast_empty.0 + $phast_full.0) *100)]
puts " "
#puts "There are $phast_full variants located in a conserved region annotated by phastCons "
#puts "There are $phast_empty variants not located in such a region"
puts " $pro % of the variants are located in a conserved region annotated by phastCons"

#############
#
# oreganno
#
#############
set pro [expr  ( $or_full.0 / ($or_empty.0 + $or_full.0) *100)]
puts " "
#puts "There are $or_full variants located in a regulatory element annotated by ORegAnno "
#puts "There are $or_empty variants not located in such an element"
puts " $pro % of the variants are located in a regulatory element annotated by ORegAnno "


#############
#
# dbsnp
#
#############
set pro [expr  ( $db_full.0 / ($db_empty.0 + $db_full.0) *100)]
puts " "
#puts "There are $db_full variants known in dbsnp database "
#puts "There are $db_empty variants are not known in dbsnp database"
puts " $pro % of the variants are known in dbsnp database"

#############
#
# UTR annotation from complete genomics
# 
#############

#puts " "
#puts "There are $ex_utr variants annotated UTR by complete genomics"
#puts "There are $ex_ex variants annotated EXON by complete genomics"


tempfile clean
exit 0