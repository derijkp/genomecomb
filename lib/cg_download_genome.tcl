#
# Copyright (c) by Peter De Rijk (VIB - University of Antwerp)
# See the file "license.txt" for information on usage and redistribution of
# this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc cg_download_genome {args} {
	set chromosomes {}
	set alt 0
	set url {}
	cg_options download_genome args {
		-chromosomes {set chromosomes $value}
		-alt {set alt $value}
		-url {set url $value}
	} {result build chromosomes} 2 3
	set keepdir [pwd]
	set result [file_absolute $result]
	file mkdir $result.temp
	cd $result.temp
	set tail [file tail $result]
	putslog "Downloading $build genome"
	set source ""
	if {$url ne ""} {
		set source $url
		wgetfile $source
		cg_fas2ifas [file tail $source] $tail
		file rename -force -- $tail $tail.index ..
	} elseif {[llength $chromosomes]} {
		set source ftp://hgdownload.cse.ucsc.edu/goldenPath/$build/chromosomes/chr*.fa.gz
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
		file rename -force -- {*}[glob $tail*] ..
	} elseif {!$alt && ![catch {
		set source ftp://hgdownload.cse.ucsc.edu/goldenPath/$build/bigZips/analysisSet/$build.analysisSet.chroms.tar.gz
		wgetfile $source 
		if {![file exists $build.analysisSet.chroms.tar.gz]} {error "Could not download $build chromosomes"}
	} msg]} {
		exec tar xvzf $build.analysisSet.chroms.tar.gz
		set files [bsort [glob $build.analysisSet.chroms/*.fa]]
		if {!$alt} {
			set files [list_sub $files -exclude [list_find -regexp $files _hap]]
		}
		putslog "Converting and indexing"
		exec cat {*}$files | cg genome_indexfasta $tail
		file rename -force -- {*}[glob $tail*] ..
		file delete $build.analysisSet.chroms.tar.gz
	} elseif {![catch {
		set source ftp://hgdownload.cse.ucsc.edu/goldenPath/$build/chromosomes/chr*.fa.gz
		wgetfiles $source goldenPath-$build-chromosomes
		if {![llength [glob goldenPath-$build-chromosomes/*.fa.gz]]} {error "Could not download $build chromosomes"}
	} msg]} {
		set files [bsort [glob goldenPath-$build-chromosomes/*.fa.gz]]
		if {!$alt} {
			set files [list_sub $files -exclude [list_find -regexp $files _hap]]
		}
		putslog "Converting and indexing"
		exec zcat {*}$files | cg genome_indexfasta $tail
		file rename -force -- {*}[glob $tail*] ..
	} else {
		set source ftp://hgdownload.cse.ucsc.edu/goldenPath/$build/bigZips/$build.fa.gz
		wgetfile $source
		cg_fas2ifas $build.fa.gz $tail
		file rename -force -- $tail $tail.index ..
	}
	puts "genome downloaded from $source"
	file_write $result.info [subst [deindent {
		= Genome (build $build) =
		
		== Download info ==
		dbname	genome
		version	[timestamp]
		license	free
		source	$source
		time	[timestamp]
		
		== Description ==
		The genome sequence downloaded from source, and chromosome sorted according 
		to a natural sort (for use in genomecomb)
		
		== Category ==
		Genome
	}]]
	cd ..
	catch {file delete -force $result.temp}
	#
	# make samtools index
	exec samtools faidx $result
	#
	set rfile [file dir $result]/reg_[file root $tail].tsv
	putslog "Making $rfile"
	set data [file_read $result.fai]
	set o [open $rfile.temp w]
	puts $o chromosome\tbegin\tend
	list_foreach {chromosome len} [lrange [split [string trim $data] \n] 0 end] {
		puts $o $chromosome\t0\t$len
	}
	close $o
	file rename -force -- $rfile.temp $rfile
	cd $keepdir
}
