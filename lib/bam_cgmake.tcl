target bam_index {(.*).bam.bai} -deps {$target1.bam} -code {
	exec samtools index $dep >@ stdout 2>@ stderr
}

target bam_idxstats {(.*).bam.idxstats} -deps {$target1.bam} -code {
	exec samtools idxstats $dep > $target 2>@ stderr
}

proc bam2coverage {bamfile coverageprefix} {
	set base $coverageprefix
	if {[file exists $base.FINISHED]} {
		puts "coverage: $base already done"
		return
	}
	puts "coverage: making $base"
	exec samtools depth $bamfile | cg bcol make -p pos -c contig -t su $base-coverage coverage
	set files [glob -nocomplain $base-coverage-*.bcol]
	foreach file $files {
		# samtools depth is 1 based, this hack will correct this (without having to remap the entire coverage file)
		set c [split [string trim [file_read $file]] \n]
		set line [lindex $c end-1]
		set temp [lindex $line 0]
		lset line 0 [incr temp -1]
		lset c end-1 [join $line \t]
		set line [lindex $c end]
		set temp [lindex $line 0]
		lset line 0 [incr temp -1]
		lset c end [join $line \t]
		write_file $file [join $c \n]
	}
	write_file $base.FINISHED ""
}


target bam2coverage {(.*)/([^/]*).coverage-/coverage-([^/-]*).bcol} -deps {$target1/*.bam} -code {
}


# Added for autoloading lib with cgmake_lib
proc cgmakelib_bam {args} {
	return targetdir
}
