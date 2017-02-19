#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_download_1000g3 {args} {
	set url ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
	cg_options download_1000g3 args {
	} {resultfile build url} 2 3 {
		download data from 1000 genomes final release
	}
	if {$build ne "hg19"} {
		error "only build hg19 supported"
	}
	set dir [file tail [file dir $url]]
	set tail [file tail $url]
	set base [file root [gzroot $tail]]
	if {[file exists $resultfile]} {
		putslog "skipping file $resultfile: exists"
		return
	}
	set tempdir $resultfile.temp
	file mkdir $tempdir
	# info
	regsub {/[^/]+$} $url {} urldir
	set error [catch {
		wgetfile $urldir/README_phase3_callset_20150220 $tempdir/[file tail $infourl]
	} msg]
	if {$error} {
		wgetfiles $urldir/README_* $tempdir
	}
	set readme [lindex [glob -nocomplain $tempdir/README_phase3_callset_20150220 $tempdir/README_*callset* $tempdir/README_*] 0]
	set o [open $tempdir/result.info w]
	puts $o dbname\t1000g3
	puts $o "version\t$dir"
	puts $o "citation\tThe 1000 Genomes Project Consortium. 2015. A Global Reference for Human Genetic Variation. Nature 526 (7571): 68\u201374. doi:10.1038/nature15393."
	puts $o "website\thttp://1000genomes.org"
	puts $o "source\t$url"
	puts $o "time\t[timestamp]"
	puts $o ""
	close $o
	exec cat $readme >> $tempdir/result.info
	file rename -force $tempdir/result.info $resultfile.info
	
	# data
	wgetfile $url $tempdir/$tail
	set fields {chromosome begin end type ref alt {freq=vformat("%.3f",($allelecount @* 1.0) @/ $totalallelecount)} quality filter totalallelecount AMR_AF EAS_AF AFR_AF EUR_AF SAS_AF}
	set error [catch {
		exec cg vcf2tsv $tempdir/$tail $tempdir/$base.tsv.temp
	} msg]
	if {$error} {
		set msg [split $msg \n]
		set len [llength [list_find -regexp $msg {not described in header|broken pipe|child process exited abnormally}]]
		if {$len ne [llength $msg]} {
			error "error converting $tempdir/$base.vcf.gz:\n$msg"
		}
	}
	file rename -force $tempdir/$base.tsv.temp $tempdir/$base.tsv
	cg select -s - -f $fields $tempdir/$base.tsv $tempdir/result.tsv.temp
	file rename -force $tempdir/result.tsv.temp $resultfile
	file delete -force $resultfile.temp
}

proc cg_download_1000glow {args} {
	set url ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz
	cg_options download_1000glow args {
	} {resultfile build url} 2 3 {
		download data from 1000 genomes phase1 (low coverage) release
	}
	if {$build ne "hg19"} {
		error "only build hg19 supported"
	}
	if {[file exists $resultfile]} {
		putslog "skipping file $resultfile: exists"
		return
	}
	set tempdir $resultfile.temp
	file mkdir $tempdir
	# catch {exec wget --tries=45 --directory-prefix=$tempdir/ ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20101123/interim_phase1_release/ALL.wgs.phase1.projectConsensus.snps.sites.vcf.gz >@stdout 2>@stderr} errmsg
	# catch {exec wget -c --tries=45 --directory-prefix=$tempdir/ ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20110521/ALL.wgs.merged_beagle_mach.20101123.snps_indels_svs.sites.vcf.gz >@stdout 2>@stderr} errmsg
	# catch {exec wget -c --tries=45 --directory-prefix=$tempdir/ ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20110521/ALL.wgs.phase1_release_v2.20101123.snps_indels_sv.sites.vcf.gz >@stdout 2>@stderr} errmsg
	catch {exec wget -c --tries=45 --directory-prefix=$tempdir/ $url >@stdout 2>@stderr} errmsg
	cg vcf2sft $tempdir/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz $tempdir/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.tsv.temp
	file rename -force $tempdir/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.tsv.temp $tempdir/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.tsv
	cg select -s - -f {chromosome begin end type ref alt {freq=format("%.2f",double($allelecount)/$totalallelecount)} quality filter totalallelecount AMR_AF ASN_AF AFR_AF EUR_AF} $tempdir/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.tsv $resultfile.temp
	file rename -force $resultfile.temp $resultfile
	file delete -force $tempdir
}

proc cg_download_1000g {args} {
	cg_options download_1000g args {
	} {resultfilebase build url} 2 3 {
		download data from 1000 genomes pilot
	}
	if {$build ne "hg18"} {
		error "only build hg18 supported"
	}
	set tempdir $resultfilebase.temp
	file mkdir $tempdir
	set populations {CEU CHB+JPT YRI}
	foreach pop $populations {
		regsub -all {\+} $pop x rpop
		set resultfile $resultfilebase$rpop.tsv
		if {[file exists $resultfile]} {
			puts "$resultfile exists: skipping"
			continue
		}
		puts "Making $resultfile"
		if {![file exists $tempdir/$pop.SRP000031.2010_03.sites.vcf.gz]} {
			puts "downloading $pop.SRP000031.2010_03.sites.vcf.gz"
			wgetfile ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/release/2010_03/pilot1/$pop.SRP000031.2010_03.sites.vcf.gz $tempdir/$pop.SRP000031.2010_03.sites.vcf.gz
		}
#		if {![file exists $tempdir/$pop.SRP000031.2010_03.indels.sites.vcf.gz]} {
#			puts "downloading $pop.SRP000031.2010_03.indels.sites.vcf.gz"
#			wgetfile ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/release/2010_03/pilot1/indels/$pop.SRP000031.2010_03.indels.sites.vcf.gz $tempdir/$pop.SRP000031.2010_03.indels.sites.vcf.gz
#		}
		puts "Changing format"
		catch {close $f} ; catch {close $o}
		set f [gzopen $tempdir/$pop.SRP000031.2010_03.sites.vcf.gz]
		set o [open $tempdir/var_${build}_1000g$pop.tsv w]
		while {![eof $f]} {
			set line [gets $f]
			if {[string index $line 0] ne "#"} break
		}
		set line [split $line \t]
		puts $o [join {chrom start end type ref alt freq id} \t]
		set num 0; set next 100000
		set type snp
		while {![eof $f]} {
			incr num; if {$num >= $next} {puts $num; incr next 100000}
			if {![llength $line]} {
				set line [split [gets $f] \t]
				continue
			}
			foreach {chrom end id ref alt} $line break
			set start [expr {$end-1}]
			set temp [lindex $line end]
			regexp {AC=([0-9]+)} $temp t ac
			regexp {AN=([0-9]+)} $temp t an
			set freq [format %.3f [expr {$ac/double($an)}]]
			puts $o "chr$chrom\t$start\t$end\t$type\t$ref\t$alt\t$freq\t$id"
			set line [split [gets $f] \t]
		}
		close $o
		gzclose $f
		puts "Sorting $resultfile"
		cg select -s {chrom start end type alt} $tempdir/var_${build}_1000g$pop.tsv $tempdir/result.tsv.temp
		file rename -force $tempdir/result.tsv.temp $resultfile
	}
	file delete -force $tempdir
}
