#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_download_genome {args} {
	set chromosomes {}
	set alt 0
	cg_options download_genome args {
		-chromosomes {set chromosomes $value}
		-alt {set alt $value}
	} {result build chromosomes} 2 3
	set keepdir [pwd]
	set result [file_absolute $result]
	file mkdir $result.temp
	cd $result.temp
	set tail [file tail $result]
	putslog "Downloading $build genome"
	if {[llength $chromosomes]} {
		set files {}
		foreach chr $chromosomes {
			lappend files chr$chr.fa.gz
			if {[file exists chr$chr.fa.gz]} {
				puts "Skipping chromosome $chr: chr$chr.fa.gz already there"
				continue
			}
			putslog "Downloading chromosome chr$chr.fa.gz"
			wgetfile ftp://hgdownload.cse.ucsc.edu/goldenPath/$build/chromosomes/chr$chr.fa.gz
		}
		putslog "Converting and indexing"
		exec zcat {*}$files | cg genome_indexfasta $tail
		file rename -force {*}[glob $tail*] ..
	} elseif {!$alt && ![catch {
		exec wget --tries=45 ftp://hgdownload.cse.ucsc.edu/goldenPath/$build/bigZips/analysisSet/$build.analysisSet.chroms.tar.gz
	} msg]} {
		exec tar xvzf $build.analysisSet.chroms.tar.gz
		set files [ssort -natural [glob $build.analysisSet.chroms/*.fa]]
		if {!$alt} {
			set files [list_sub $files -exclude [list_find -regexp $files _hap]]
		}
		putslog "Converting and indexing"
		exec cat {*}$files | cg genome_indexfasta $tail
		file rename -force {*}[glob $tail*] ..
	} elseif {![catch {
		wgetfiles ftp://hgdownload.cse.ucsc.edu/goldenPath/$build/chromosomes/chr*.fa.gz goldenPath-$build-chromosomes
		if {![llength [glob goldenPath-$build-chromosomes/*.fa.gz]]} {error "Could not download $build chromosomes"}
	} msg]} {
		set files [ssort -natural [glob goldenPath-$build-chromosomes/*.fa.gz]]
		if {!$alt} {
			set files [list_sub $files -exclude [list_find -regexp $files _hap]]
		}
		putslog "Converting and indexing"
		exec zcat {*}$files | cg genome_indexfasta $tail
		file rename -force {*}[glob $tail*] ..
	} else {
		exec wget --tries=45 ftp://hgdownload.cse.ucsc.edu/goldenPath/$build/bigZips/$build.fa.gz >@ stdout  2>@ stderr
		cg_fas2ifas $build.fa.gz $tail
		file rename -force $tail $tail.index ..
	}
	cd ..
	catch {file delete -force $result.temp}
	#
	# make samtools index
	exec samtools faidx $result
	#
	set rfile [file dir $result]/reg_[file root $tail].tsv
	putslog "Making $rfile"
	set data [file_read $result.index]
	set o [open $rfile.temp w]
	puts $o chromosome\tbegin\tend
	list_foreach {chromosome begin len} [lrange [split [string trim $data] \n] 0 end-1] {
		puts $o $chromosome\t0\t$len
	}
	close $o
	file rename -force $rfile.temp $rfile
	cd $keepdir
}
