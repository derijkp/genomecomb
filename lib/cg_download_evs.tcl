proc cg_download_evs {args} {
	cg_options download_phenotype args {
	} {resultfile url}
	set resultfile [file_absolute $resultfile]
	if {$url eq ""} {
		set url http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.protein-hgvs-update.snps_indels.vcf.tar.gz
		set evsbuild hg19
		if {$build ne "hg19"} {error "for build other than hg19, the url must be given"}
	}
	set version {}
	regexp {^[^.]*} [file tail $url] version
	set tempdir $resultfile.temp
	set tempresult $tempdir/result.tsv
	file mkdir $tempdir
	set keeppwd [pwd]
	cd $tempdir
	if {![file exists evs.tar.gz]} {
		wgetfile $url evs.tar.gz
	}
	putslog "Unpacking evs file"
	exec tar xvzf evs.tar.gz
	putslog "Converting to tsv"
	foreach file [glob *.vcf] {
		cg vcf2tsv $file [file root $file].tsv
	}
	set files [lsort -dict [glob ESP6500SI-*.snps_indels.tsv]]
	if {[llength $files] != 24} {error "not enough files found"}
	putslog "Concatenating files"
	cg cat {*}$files > $tempresult.temp
#	cg select -f {chromosome begin end type 
#		{ea_freqp=lrange(vformat("%.3f",(100.0 @* $EA_AC) @/ lsum($EA_AC)),0,"end-1")} 
#		{aa_freqp=lrange(vformat("%.3f",(100.0 @* $AA_AC) @/ lsum($AA_AC)),0,"end-1")}
#		*
#	} $tempresult.temp $tempresult.temp2
	putslog "Create results"
	set f [open $tempresult.temp]
	set header [tsv_open $f]
	set bposs [tsv_basicfields $header]
	set o [open $tempresult.temp2 w]
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
	file_write [gzroot $resultfile].opt "fields\t{ea_freqp aa_freqp ea_mfreqp aa_mfreqp}\n"
	#
file_write [gzroot $resultfile].info [subst [string trim {
Exome variant server
====================

Download info
-------------
dbname	evs
version	$version
citation	Exome Variant Server, NHLBI GO Exome Sequencing Project (ESP), Seattle, WA (URL: http://evs.gs.washington.edu/EVS/) \[[lindex [timestamp] 0]\].
website	http://evs.gs.washington.edu/EVS/
source	$url
time	[timestamp]

Data Usage
----------
We request that any use of data obtained from the NHLBI GO ESP Exome Variant Server be cited in publications.

Citation
--------
Exome Variant Server, NHLBI GO Exome Sequencing Project (ESP), Seattle, WA (URL: http://evs.gs.washington.edu/EVS/) \[date (month, yr) accessed\].

Annotation fields added
-----------------------
ea_freqp: frequency (as percentage!) of alt aleles in european american population
aa_freqp: frequency (as percentage!) of alt aleles in african american population
ea_mfreqp: frequency (as percentage!) of homozygous alt alele genotype in european american population
aa_mfreqp: frequency (as percentage!) of homozygous alt alele genotype in african american population

more information on http://evs.gs.washington.edu/EVS/
}]]
	#
	putslog "move results to $resultfile"
	compress $tempresult.temp2 $resultfile
	cd $keeppwd
	file delete -force $tempdir
}
