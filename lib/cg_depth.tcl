proc cg_depth {args} {
	set regionfile {}
	set max 1000000
	set opts {}
	set Q 0
	set q 0
	set opts {}
	cg_options depth args {
		-max {set max $value}
		-q {set q $value}
		-Q {set Q $value}
		-all {
			if {[true $value]} {
				lappend opts -a
			}
		}
	} {bamfile} 1 1 {
		run samtools depth on a bam/sam file, outputting tsv format
	}
	puts [join {chromosome end depth} \t]
	catch_exec samtools depth {*}$opts -d$max -Q $Q -q $q $bamfile >@ stdout 2>@ stderr
}
