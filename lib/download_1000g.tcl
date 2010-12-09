if 0 {

	package require Tclx
	signal -restart error SIGINT
	lappend auto_path /home/peter/dev/completegenomics/lib
	package require Extral
	cd /complgen/refseq/hg18test
	set path /complgen/refseq/hg18test
	set build hg18
	set pop CEU

cg downloaddb /complgen/refseq/hg18test hg18 1000g

}


proc downloaddb_1000g {path build} {
	set tempdir $path/tmp/$build
	set populations {CEU CHB+JPT YRI}
	foreach pop $populations {
		set resultfile $path/1000g_$pop.tsv
		if {[file exists $resultfile]} {
			puts "$resultfile exists: skipping"
			continue
		}
		puts "Making $resultfile"
		if {![file exists $tempdir/$pop.SRP000031.2010_03.genotypes.vcf]} {
			if {![file exists $tempdir/$pop.SRP000031.2010_03.genotypes.vcf.gz]} {
				puts "downloading $pop.SRP000031.2010_03.genotypes.vcf.gz"
				catch {exec wget --tries=45 --directory-prefix=$tempdir/ ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/release/2010_03/pilot1/$pop.SRP000031.2010_03.genotypes.vcf.gz} errmsg
				if {![file exists $tempdir/$pop.SRP000031.2010_03.genotypes.vcf.gz]} {
					puts $errmsg
					exit 1
				}
			}
			puts "Unzipping $pop.SRP000031.2010_03.genotypes.vcf.gz"
			exec gunzip $tempdir/$pop.SRP000031.2010_03.genotypes.vcf.gz
		}
		if {![file exists $tempdir/$pop.frq]} {
			puts "Creating frequency files (This can take a long time)"
			exec vcftools1.2 --vcf $tempdir/$pop.SRP000031.2010_03.genotypes.vcf --freq --out $tempdir/$pop >@stdout 2>@ stdout
			if {![file exists $tempdir/$pop.frq]} {
				exit 1
			}
		}
		puts "Changing format"
		catch {close $f} ; catch {close $o}
		set f [open $tempdir/$pop.frq]
		set o [open $tempdir/1000g_$pop.tsv w]
		gets $f
		puts $o [join {chrom start end freq altallele numalleles frequencies} \t]
		set num 0; set next 100000
		while {![eof $f]} {
			incr num; if {$num >= $next} {puts $num; incr next 100000}
			set line [split [gets $f] \t]
			if {![llength $line]} continue
			foreach {chrom end numalleles} $line break
			set start [expr {$end-1}]
			set frequencies [lrange $line 4 end]
			set temp {}
			foreach freq $frequencies {
				lappend temp [split $freq :]
			}
			set temp [lsort -real -decreasing -index 1 $temp]
			set altallele [lindex $temp end 0]
			set freq [lindex $temp end 1]
			if {$freq == 1 || $freq == 0} continue
			if {$freq > 0.5} {
				set freq [expr {1 - $freq}]
				set altallele ?
			}
			set freq [format %.3f $freq]
			puts $o "chr$chrom\t$start\t$end\t$freq\t$altallele\t$numalleles\t$frequencies"
		}
		close $o
		close $f
		puts "Sorting $resultfile"
		cg select -s {chrom start end} $tempdir/1000g_$pop.tsv $resultfile.temp
		file rename $resultfile.temp $resultfile
	}
}
