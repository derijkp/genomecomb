if 0 {

	package require Tclx
	signal -restart error SIGINT
	lappend auto_path /home/peter/dev/completegenomics/lib
	package require Extral
	cd /complgen/refseq/hg18
	set path /complgen/refseq/hg18
	set build hg18
	set pop CEU

cg downloaddb /complgen/refseq/hg18 hg18 1000g

}

proc downloaddb_1000g {path build} {
	if {$build ne "hg18"} {
		error "only build hg18 supported"
	}
	set tempdir $path/tmp/$build
	set populations {CEU CHB+JPT YRI}
	foreach pop $populations {
		set resultfile $path/reg_${build}_1000g$pop.tsv
		if {[file exists $resultfile]} {
			puts "$resultfile exists: skipping"
			continue
		}
		puts "Making $resultfile"
		if {![file exists $tempdir/$pop.SRP000031.2010_03.sites.vcf.gz]} {
			puts "downloading $pop.SRP000031.2010_03.sites.vcf.gz"
			catch {exec wget --tries=45 --directory-prefix=$tempdir/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/release/2010_03/pilot1/$pop.SRP000031.2010_03.sites.vcf.gz} errmsg
			if {![file exists $tempdir/$pop.SRP000031.2010_03.sites.vcf.gz]} {
				puts $errmsg
				exit 1
			}
		}
		puts "Changing format"
		catch {close $f} ; catch {close $o}
		set f [rzopen $tempdir/$pop.SRP000031.2010_03.sites.vcf.gz]
		set o [open $tempdir/reg_${build}_1000g$pop.tsv w]
		while {![eof $f]} {
			set line [gets $f]
			if {[string index $line 0] ne "#"} break
		}
		set line [split $line \t]
		puts $o [join {chrom start end freq ref alt id} \t]
		set num 0; set next 100000
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
			puts $o "chr$chrom\t$start\t$end\t$freq\t$ref\t$alt\t$id"
			set line [split [gets $f] \t]
		}
		close $o
		close $f
		puts "Sorting $resultfile"
		cg select -s {chrom start end} $tempdir/reg_${build}_1000g$pop.tsv $resultfile.temp
		file rename $resultfile.temp $resultfile
	}
}
