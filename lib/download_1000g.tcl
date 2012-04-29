#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc downloaddb_1000glow {path build} {
	if {$build ne "hg19"} {
		error "only build hg19 supported"
	}
	set resultfile $path/$build/var_${build}_1000glow.tsv
	if {[file exists $resultfile]} {
		puts stderr "skipping file $resultfile: exists"
		return
	}
	set tempdir $path/tmp/$build
	# catch {exec wget --tries=45 --directory-prefix=$tempdir/ ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20101123/interim_phase1_release/ALL.wgs.phase1.projectConsensus.snps.sites.vcf.gz >@stdout 2>@stderr} errmsg
	# catch {exec wget -c --tries=45 --directory-prefix=$tempdir/ ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20110521/ALL.wgs.merged_beagle_mach.20101123.snps_indels_svs.sites.vcf.gz >@stdout 2>@stderr} errmsg
	catch {exec wget -c --tries=45 --directory-prefix=$tempdir/ ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/release/20110521/ALL.wgs.phase1_release_v2.20101123.snps_indels_sv.sites.vcf.gz >@stdout 2>@stderr} errmsg

	cg vcf2sft $tempdir/ALL.wgs.phase1_release_v2.20101123.snps_indels_sv.sites.vcf.gz $tempdir/ALL.wgs.phase1_release_v2.20101123.snps_indels_sv.sites.tsv.temp
	file rename $tempdir/ALL.wgs.phase1_release_v2.20101123.snps_indels_sv.sites.tsv.temp $tempdir/ALL.wgs.phase1_release_v2.20101123.snps_indels_sv.sites.tsv
	cg select -s {chromosome begin end type alt} -f {chromosome begin end type ref alt {freq=$allelecount/$totalallelecount} quality filter} $tempdir/ALL.wgs.phase1_release_v2.20101123.snps_indels_sv.sites.tsv $resultfile
}

proc downloaddb_1000g {path build} {
	if {$build ne "hg18"} {
		error "only build hg18 supported"
	}
	set tempdir $path/tmp/$build
	set populations {CEU CHB+JPT YRI}
	foreach pop $populations {
		regsub -all {\+} $pop x rpop
		set resultfile $path/$build/var_${build}_1000g$rpop.tsv
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
		close $f
		puts "Sorting $resultfile"
		cg select -s {chrom start end type alt} $tempdir/var_${build}_1000g$pop.tsv $resultfile.temp
		file rename $resultfile.temp $resultfile
	}
}
