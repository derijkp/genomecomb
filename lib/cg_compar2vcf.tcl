proc cg_compar2vcf {args} {
	if {([llength $args] < 0) || ([llength $args] > 2)} {
		errorformat exportvcf
		exit 1
	}
	if {[llength $args] > 0} {
		set filename [lindex $args 0]
		set f [gzopen $filename]
		set samples [cg select -n $filename] 
	} else {
		set f stdin
		set samples [cg select -n stdin] 
	}
	if {[llength $args] > 1} {
		set outfile [lindex $args 1]
		set o [open $outfile w]
	} else {
		set o stdout
	}
	
	set exportVAT 1
	
	###get header
	set header [tsv_open $f]
	###basic fields
	set bfield_pos [tsv_basicfields $header 6]
	set bfields [list_sub $header $bfield_pos]
	###get unique sample fields (i.e. part before 1st dash of fields containing one or multiple dashes)
	set index [list_find -regexp $header {-}]
	set sfields [list_remdup [list_regsub {(^[^-]+)-.+$} [list_sub $header $index] {\1}]]
	
	################
	###write header
	################
	set vcf_info {##fileformat=VCFv4.0}
	
	if {$exportVAT} {
		set info {
			{##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">}
			{##INFO=<ID=HM2,Number=0,Type=Flag,Description="HapMap2 membership">}
			{##INFO=<ID=HM3,Number=0,Type=Flag,Description="HapMap3 membership">}
			{##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">}
			{##INFO=<ID=AC,Number=1,Type=Integer,Description="total number of alternate alleles in called genotypes">}
			{##INFO=<ID=AN,Number=1,Type=Integer,Description="total number of alleles in called genotypes">} 
			{##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">}
			{##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth from MOSAIK BAM">}
		}
			
		foreach i $info {
			append vcf_info \n$i
		}
		
		set column_header [join {#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT} \t]
		
		
		foreach s $samples {
			append column_header \t$s
		}

		puts $o $vcf_info
		puts $o $column_header
	}	
	
	#####################################
	###add data
	#####################################
	
	set allele_pos [list_find -regexp $header {alleleSeq}]
	set cov_pos [list_find -regexp $header {^coverage-}]
	set sequenced_alleleSeq_coverage_pos [list_find -regexp $header {sequenced-|alleleSeq|^coverage-}]
	set sequenced_pos [list_find -regexp $header {sequenced-}]
	
	if {$exportVAT} {
		while {1} {
			set line [split [gets $f] \t]
			if {![eof $f]} {
				
				#skip line if variant is same as previous line  -- temp solution to deal with splitted compar files
				set var [list_sub $line [lrange $bfield_pos 0 3]]
				
				if { ![info exists old_var] || $var ne $old_var} {
					if {[list_sub $line [lindex $bfield_pos 3]] eq "snp"} {
						#set general info
						set chrom [list_sub $line [lindex $bfield_pos 0 ]]
						set pos [list_sub $line [lindex $bfield_pos 2 ]]
						set ref [list_sub $line [lindex $bfield_pos 4 ]]
						
						#get (alternative) alleles
						set alleles_list [list_sub $line $allele_pos]
						set alt_alleles_list [list_remove $alleles_list $ref -]
						set alt_alleles [list_remdup $alt_alleles_list]
						set alleles [list $ref {*}$alt_alleles]
						
						#prepare INFO
						#DP = total coverage, including depth for samples where no call could be made, i.e. with GT ./.
						#AC = count of occurrence of different alt alleles in sample genotypes in order listed
						#AN = total number of alleles called -> # of samples that are not unsequenced * 2
						#AA = ancestral allele -> base found in ancestor, e.g. chimp	
						set total_cov 0
						foreach i $cov_pos {
							set cov [lindex $line $i]
							if {[regexp {[0-9]+} $cov] } {set total_cov [expr $total_cov + $cov] } 
						}
						
						set alt_counts {}
						foreach a $alt_alleles {
							 lappend alt_counts [llength [list_find $alt_alleles_list $a]]
						}
						set alt_counts [join $alt_counts ","]
						
						set total_allele_num 0
						foreach i $sequenced_pos {
							if {[lindex $line $i] != "u"} {
								set total_allele_num [expr $total_allele_num + 2]
							}
						}
							
						set AA AA=.
						set AC AC=$alt_counts
						set AN AN=$total_allele_num
						set DP DP=$total_cov
						set INFO [join [list $AA $AC $AN $DP] ";"]
						
						
						#write general part new line
						set new_line {} 
						#lappend \n
						
						lappend new_line $chrom
						lappend new_line $pos
						lappend new_line .
						lappend new_line $ref
						lappend new_line [join $alt_alleles ","]
						lappend new_line .
						lappend new_line PASS
						lappend new_line $INFO
						lappend new_line GT:DP
						
						#write sample specific part new line
						foreach {s a1p a2p cp} $sequenced_alleleSeq_coverage_pos {
							if {[lindex $line $s] == "u"} {
								set GT ./.
							} else {
								if {[set a1 [list_cor $alleles [lindex $line $a1p]]] == -1 } { set a1 "."}
								if {[set a2 [list_cor $alleles [lindex $line $a2p]]] == -1 } { set a2 "."}
								set GT $a1/$a2
							}
							if {![regexp {[0-9]+} [set c [lindex $line $cp]]] } {set c 0} 
							lappend new_line $GT:$c
						}
						puts $o [join $new_line \t]
					}
				}
					
				set old_var $var 
				
			} else {
				break
			}
		} 
	}
	
	if {$o ne "stdout"} {catch {close $o}}
	if {$f ne "stdin"} {catch {close $f}}
	
}	






