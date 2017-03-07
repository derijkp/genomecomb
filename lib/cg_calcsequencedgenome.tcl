proc cg_calcsequencedgenome {args} {
	cg_options calcsequencedgenome args {
	} {genomefile} 1 1 {
		calculate sequenced genome from genome ifas (regions without Ns)
	}
	set f [open $genomefile]
	set o stdout
	puts $o chromosome\tbegin\tend\tsize
	while {![eof $f]} {
		set name [gets $f]
		if {$name eq ""} continue
		set chr [chr_clip [lindex [string range $name 1 end] 0]]
		putslog $name\n$chr
		set seq [gets $f]
		set indices [regexp -all -inline -indices {[^N]{1,}} $seq]
		putslog Writing
		list_foreach {pbegin pend} $indices break
		list_foreach {begin end} $indices {
			incr end
			if {$begin > [expr {$pend+1}]} {
				puts $o chr$chr\t$pbegin\t$pend\t[expr {$end-$begin}]
				set pbegin $begin
			}
			set pend $end
		}
		puts $o chr$chr\t$pbegin\t$pend\t[expr {$end-$begin}]
	}
	close $f
}
