proc cg_calcsequencedgenome {args} {
	cg_options calcsequencedgenome args {
	} {genomefile resultfile} 2 2 {
		calculate sequenced genome from genome ifas (regions without Ns)
	}
	set f [open $genomefile]
	file mkdir [file dir $resultfile]
	set tempresult [filetemp $resultfile]
	set o [open $tempresult w]
	puts $o chromosome\tbegin\tend\tsize
	while {![eof $f]} {
		set name [gets $f]
		if {$name eq ""} continue
		if {![regexp {chromosome ([0-9A-Z]+)} $name temp chr]} {
			if {![regexp {chr([0-9A-Z]+)} $name temp chr]} {
				set chr $name
			}
		}
		putslog $name\n$chr
		set seq [gets $f]
		set indices [regexp -all -inline -indices {[^N]{1,}} $seq]
		putslog Writing
		list_foreach {pbegin pend} $indices break
		list_foreach {begin end} $indices {
			incr end
			if {$begin > [expr {$pend+1}]} {
				puts $o chr$chr\t$pbegin\t$pend
				set pbegin $begin
			}
			set pend $end
		}
		puts $o chr$chr\t$pbegin\t$pend
	}
	close $o
	close $f
	file rename -force $tempresult $resultfile
}
