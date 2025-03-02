proc find_umis {fastq target {adaptorseq GATCGGAAGAGCACACGTCTGAACTCCAGTCAC} {barcodesize 8} {umisize 12} {ASlimit 0}} {
	set reffile [tempfile].fa
	set ref [sc_barcodes_ref $reffile $adaptorseq]
	# set resultdir [file dir $target]
	# set sam [tempfile].umi.sam
	set temptarget $target.temp[gzext $target]
	set sam $temptarget.umi_sam
	if {[file ext $fastq] in ".bam .cram .sam"} {
		set usefastq [tempfile].fastq.gz
		catch_exec samtools fastq -T "RG,CB,QT,MI,MM,ML,Mm,Ml" $fastq | gzip > $usefastq
		set ubams 1
	} else {
		set usefastq $fastq
		set ubams 0
	}
	# -Y In SAM output, use soft clipping for supplementary alignments
	# -t threads
	# -k k-mer size
	# -w window size
	# -n Discard chains consisting of <INT number of minimizers
	# -m Discard chains with chaining score <INT (def 40)
	# -s Minimal peak DP alignment score to output (def 40)
	# -A Matching score [2]
	# -B Mismatching penalty [4]
	# -O Gap open penalty
	# -E Gap extension penalty
	catch_exec minimap2 -Y -a --secondary=no -x map-ont -Y -t 4 -n 1 -k 5 -w 1 -m 10 -s 20 -A 2 -B 4 -O 2 -E 2 $ref $usefastq > $sam 2>@ stderr
	if {$ubams} {file delete $usefastq}

	catch {close $ff} ; catch {close $f} ; catch {close $o}
	unset -nocomplain a
	unset -nocomplain ba
	set a() 0
	set ba() 0
	set ff [gzopen $fastq]
	set ffline [gets $ff]
	set ffid [lindex [string range $ffline 1 end] 0]
	set f [open "| cg sam2tsv -fields AS $sam"]
	set o [wgzopen $temptarget]
	set header [tsv_open $f]
	set poss [list_cor $header {chromosome begin end strand qname qstart qend cigar seq supplementary mapquality AS}]
	set qnamepos [lsearch $header qname]
	set mqpos [lsearch $header mapquality]
	set strandpos [lsearch $header strand]
	set qstartpos [lsearch $header qstart]
	set num 0
	set todo {}
	set prevqname {}
	while 1 {
		if {[gets $f nline] == -1} break
		set nline [split $nline \t]
		if {![expr [incr num]%10000]} {puts $num}
		set nqname [lindex $nline $qnamepos]
		if {$nqname ne $prevqname && [llength $todo]} {
# if {$prevqname eq "f091ff67-a4ac-4963-b0ae-64d9f62d6493"} error
			# if more than one hit for adaptorseq, select "best" one
			if {[llength $todo] > 1} {
				# remove hits with lower mapping quality
				set qs [list_subindex $todo $mqpos]
				# set limitq [expr {[lmath_max $qs]-3}]
				set limitq [lmath_max $qs]
				set pos 0
				set keep {}
				foreach q $qs {
					if {$q >= $limitq} {lappend keep $pos}
					incr pos
				}
				if {[llength $keep] < [llength $todo]} {
					set todo [list_sub $todo $keep]
				}
			}
			if {[llength $todo] > 1} {
				# if still multiple left with similar mapquality, pick the outer one
				# (for e.g. when adaptor/read1 also included in adaptors added later, after 10x)
				if {[lindex $todo 0 $strandpos] eq "+"} {
					set line [lindex [lsort -integer -index $qstartpos $todo] 0]
				} else {
					set line [lindex [lsort -integer -index $qstartpos $todo] end]
				}
			} else {
				set line [lindex $todo 0]
			}
			foreach {chromosome begin end strand qname qstart qend cigar seq supplementary mapquality AS} [list_sub $line $poss] break
# putsvars ffid qname
# error
			if {$qname ne $ffid} {
				error "error in finding UMIs: id $qname from alignment does not match fastq id (from $fastq)"
			}
			set barcode {}
			if {$chromosome ne "*" && $AS >= $ASlimit} {
				# add UMI if match found with adaptor
				set start $qend
				if {[regexp H $cigar]} {
					error "hardclipped sequence in line: [list set line $line]"
				}
				if {$barcodesize > 0} {
					set barcode [string range $seq $start [expr {$start+$barcodesize-1}]]
					incr start $barcodesize
				}
				set umi [string range $seq $start [expr {$start+$umisize-1}]]
				incr a($umi)
				incr ba($barcode)
				if {[string length $umi]} {
					if {[string length $barcode] > 0} {
						set ffline @${barcode}_$umi\#[string range $ffline 1 end]\ CR:Z:$barcode\ MI:Z:$umi
					} else {
						set ffline @$umi\#[string range $ffline 1 end]\ MI:Z:$umi
					}
				}
			} else {
				incr a()
				incr ba()
			}
			puts $o $ffline
			set ffline [gets $ff]
			puts $o $ffline
			set ffline [gets $ff]
			puts $o $ffline
			set ffline [gets $ff]
			puts $o $ffline
			set ffline [gets $ff]
			set ffid [lindex [string range $ffline 1 end] 0]
			set todo [list]
		}
		# if {$supplementary} continue
		set prevqname $nqname
		lappend todo $nline
	}

	gzclose $o
	gzclose $f
	file rename -force $temptarget $target
	#
	if {[info exists sumresultfile]} {
		set o [wgzopen $sumresultfile.temp w]
		puts $o umi\tcount
		foreach umi [array names a] {
			puts $o $umi\t$a($umi)
		}
		gzclose $o
		cg select -overwrite 1 -s count $sumresultfile.temp $sumresultfile.temp2
		file rename -force $sumresultfile.temp2 $sumresultfile
		file delete $sumresultfile.temp
		set o [wgzopen $sumresultfile.tempb w]
		puts $o barcode\tcount
		foreach barcode [array names ba] {
			puts $o $barcode\t$ba($barcode)
		}
		gzclose $o
		cg select -overwrite 1 -s count $sumresultfile.tempb $sumresultfile.tempb2
		file rename -force $sumresultfile.tempb2 [file_root [gzroot $sumresultfile]]_bc[file ext [gzroot $sumresultfile]][gzext $sumresultfile]
		file delete $sumresultfile.tempb
	}

#	file delete $sam
}

proc add_umis_job args {
	upvar job_logdir job_logdir
	global appdir
	set cmdline [clean_cmdline cg add_umis {*}$args]
	set refseq {}
	# = cannot use P7 because a sizeable proportion of reads misses it
	# = Read2 (followed by barcode, then umi, then P7)
	set adaptorseq GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
	set barcodesize 8
	set umisize 12
	set ASlimit 20
	# set maxfastqdistr {} # todo
	set skips {}
	cg_options add_umis args {
		-adaptorseq {
			set adaptorseq $value
		}
		-barcodesize {
			set barcodesize $value
		}
		-umisize {
			set umisize $value
		}
		-ASlimit {
			set ASlimit $value
		}
		-skip {
			lappend skips -skip $value
		}
	} {fastqdir resultdir} 1 2 {
	}
	# logfile
	set fastqdir [file_absolute $fastqdir]
	if {![info exists resultdir]} {set resultdir [file dir $fastqdir]}
	set resultdir [file_absolute $resultdir]
	set sample [file tail $resultdir]
	# set workdir [shadow_workdir $resultdir]
	# set workdir [workdir $resultdir]
	set cleanupfiles {}
	if {![file isdir $fastqdir]} {
		error "$fastqdir is not a directory"
	}
	mkdir $resultdir
	job_logfile $resultdir/sc_barcodes_[file tail $resultdir] $resultdir $cmdline \
		{*}[versions minimap2]
	set fastqs [gzfiles $fastqdir/*.fq $fastqdir/*.fastq $fastqdir/*.bam $fastqdir/*.cram $fastqdir/*.sam]
	foreach fastq $fastqs {
		set root [file root [gzroot [file tail $fastq]]]
		set target $resultdir/$root-umi.fastq.gz
		job find_umis-$sample-[file tail $fastq] {*}$skips \
		  -deps {
			$fastq
		} -targets {
			$target
		} -vars {
			 adaptorseq fastq barcodesize umisize ASlimit
		} -procs {
			find_umis
		} -code {
			find_umis $fastq $target $adaptorseq $barcodesize $umisize $ASlimit
		}
	}
}

