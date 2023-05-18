proc sc_barcodes_ref {ref {adaptorseq CTACACGACGCTCTTCCGATCT}} {
	# (complement is AGATCGGAAGAGCGTCGTGTAG)
	set bc NNNNNNNNNNNNNNNNNNNNNNNNNN
	set post TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
	# set refseq $adaptorseq$bc$post
	set refseq $adaptorseq
	set begin [string length $adaptorseq]
	file_write $ref >adapter\n$refseq\n
	catch_exec minimap2 -x map-ont -k 5 -w 1 -d $ref.map-ont $ref
	return $ref.map-ont
}

proc find_barcodes {fastq resultfile sumresultfile adaptorseq {barcodesize 16} {umisize 12}} {
	# version 2 umisize = 10, v3 umi is 12
	set reffile [tempfile].fa
	set ref [sc_barcodes_ref $reffile $adaptorseq]
	set sam [tempfile].sam
	set sam $resultfile.sam
	catch_exec minimap2 -a --secondary=no -x map-ont -t 4 -n 1 -m 1 -k 5 -w 1 -s 20 $ref $fastq > $sam 2>@ stderr
	# cg sam2tsv $sam | cg select -g chromosome
	# exec ~/dev/genomecomb/bin/sc_getbarcodes adapter $begin 16 10 < $sam > temp
	catch {close $f} ; catch {close $o}
	unset -nocomplain a
	set f [open "|cg sam2tsv $sam"]
	set o [wgzopen $resultfile w]
	puts $o [join {id barcode umi start strand polyA} \t]
	set header [tsv_open $f]
	set poss [list_cor $header {chromosome begin end strand qname qstart qend cigar seq supplementary}]
	set num 0
	while 1 {
		if {[gets $f line] == -1} break
		set line [split $line \t]
		if {![expr [incr num]%10000]} {puts $num}
		foreach {chromosome begin end strand qname qstart qend cigar seq supplementary} [list_sub $line $poss] break
		if {$supplementary} continue
		if {$chromosome eq "*"} {
			puts $o [join [list $qname {} {} {} {} 0] \t]
			continue
		}
		set start $qend
		if {[regexp H $cigar]} {
			error "hardclipped sequence in line: [list set line $line]"
		}
		set barcode [string range $seq $start [expr {$start+$barcodesize-1}]]
		set umi [string range $seq [expr {$start+$barcodesize}] [expr {$start+$barcodesize+$umisize-1}]]
		# check Ts
		set post [string range $seq [expr {$start+$barcodesize+$umisize}] [expr {$start+$barcodesize+$umisize+14}]]
		set polya [regexp -all T $post]
#		if {$polya < 1} {
#			error "not enough Ts in line: [list set line $line]"
#		}
		if {![info exists a($barcode)]} {
			set a($barcode) 1
		} else {
			incr a($barcode)
		}
		puts $o [join [list $qname $barcode $umi $start $strand $polya] \t]
	}
	gzclose $o
	close $f
	#
	set o [wgzopen $sumresultfile.temp w]
	puts $o barcode\tcount
	foreach barcode [array names a] {
		puts $o $barcode\t$a($barcode)
	}
	gzclose $o
	cg select -overwrite 1 -s count $sumresultfile.temp $sumresultfile.temp2
	file rename -force $sumresultfile.temp2 $sumresultfile
	file delete $sumresultfile.temp
}

proc ali_ident_matrix {} {
return {distance
gap:0
indel:1
profile gap:0
    A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z
A   0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  
B   1  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  
C   1  1  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  
D   1  1  1  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  
E   1  1  1  1  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  
F   1  1  1  1  1  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  
G   1  1  1  1  1  1  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  
H   1  1  1  1  1  1  1  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  
I   1  1  1  1  1  1  1  1  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  
J   1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  
K   1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  
L   1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  
M   1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  1  1  1  1  1  1  1  1  
N   1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  1  1  1  1  1  1  1  
O   1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  1  1  1  1  1  1  
P   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  1  1  1  1  1  
Q   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  1  1  1  1  
R   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  1  1  1  
S   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  1  1  
T   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  1  
U   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  1  
V   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  1  
W   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  1  
X   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  1  
Y   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  1  
Z   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  }
}

proc sc_barcodes_job args {
	# putslog [list iso_isoquant_job {*}$args]
	upvar job_logdir job_logdir
	global appdir
	set cmdline [clean_cmdline cg sc_barcodes {*}$args]
	set whitelist {}
	set refseq {}
	set adaptorseq CTACACGACGCTCTTCCGATCT
	set barcodesize 16
	set umisize 12
	cg_options sc_barcodes args {
		-whitelist {
			set whitelist $value
		}
		-adaptorseq {
			set adaptorseq $value
		}
		-barcodesize {
			set barcodesize	$value
		}
		-umisize {
			set umisize $value
		}
	} {fastqdir resultdir} 1 2 {
		cg sc_barcodes finds (cell)barcodes and umis from 10x single cell sequencing in nanopore reads.
		I will create the results in resultdir (if not specified the parent dir of the fastqdir):
		* dir barcodeinfo will contain a tsv file (for each fastq file in fastq) with info on barcode and umi per read
		* dir bcfastq will contains fastq files where cellbarcode and umi has been added in the
		      read name (as @<cellbarcode>_<umi>#originalname) as well as in the info "fields" CB CR and MI
		* files
	}
	# logfile
	set fastqdir [file_absolute $fastqdir]
	if {![info exists resultdir]} {set resultdir [file dir $fastqdir]}
	set resultdir [file_absolute $resultdir]
	set sampledir [file dir $resultdir]
	set sample [file tail $sampledir]
	# set workdir [shadow_workdir $resultdir]
	# set workdir [workdir $resultdir]
	set cleanupfiles {}
	if {![file isdir $fastqdir]} {
		error "$fastqdir is not a directory"
	}
	job_logfile $resultdir/sc_barcodes_[file tail $resultdir] $resultdir $cmdline \
		{*}[versions minimap2]
	set fastqs [gzfiles $fastqdir/*.fq $fastqdir/*.fastq]
	set summaries {}
	mkdir $resultdir/barcodes
	foreach fastq $fastqs {
		set root [file root [gzroot [file tail $fastq]]]
		set target $resultdir/barcodes/$root.barcodes.tsv
		set target2 $resultdir/barcodes/$root.summary_barcodes.tsv
		lappend summaries $target2
		job find_barcodes-$sample-[file tail $fastq] -deps {
			$fastq
		} -targets {
			$target $target2
		} -vars {
			 adaptorseq fastq barcodesize umisize
		} -procs {
			find_barcodes
		} -code {
			find_barcodes $fastq $target $target2 $adaptorseq $barcodesize $umisize
		}
	}
	# merge barcodes and find cell barcodes
	set summaries [jobglob $resultdir/barcodes/*.summary_barcodes.tsv]
	set barcodefiles [jobglob $resultdir/barcodes/*.barcodes.tsv]
	job merge_barcodes-$sample -deps $summaries -targets {
		$resultdir/mergedbarcodes.tsv.zst
		$resultdir/barcode2celbarcode.tsv
		$resultdir/reads_per_cell.tsv
	} -vars {
		summaries resultdir whitelist
	} -procs {
		ali_ident_matrix
	} -code {

		set mergedbarcodesfile $resultdir/mergedbarcodes.tsv.zst
		set tempfile [tempfile].zst
		cg zcat {*}$barcodefiles | cg select -s {barcode * umi *} | cg zst > $tempfile
		catch {close $o} ; catch {close $f}
		set o [wgzopen $mergedbarcodesfile.temp.zst]
		puts $o barcode\tcount\tumicount
		set f [gzopen $tempfile]
		set header [tsv_open $f]
		if {$header ne "id barcode umi start strand polyA"} {error "wrong file $file?"}
		set read [gets $f line]
		foreach {id prevbarcode prevumi} [split $line \t] break
		set count 1 ; set umicount 1
		while {1} {
			foreach {id barcode umi} [split $line \t] break
			if {$barcode ne $prevbarcode} {
				if {$count > 0 && $barcode ne ""} {
					puts $o $prevbarcode\t$count\t$umicount
				}
				set prevbarcode $barcode
				set prevumi $umi
				set count 1 ; set umicount 1
			} elseif {$umi ne $prevumi} {
				incr umicount
				incr count
			} else {
				incr count
			}
			if {$read == -1} break
			set read [gets $f line]
			# putsvars barcode count umicount
		}
		gzclose $o ; gzclose $f

#		unset -nocomplain a
#		foreach file $summaries {
#			set f [gzopen $file]
#			set header [tsv_open $f]
#			if {$header ne "barcode count"} {error "wrong file $file?"}
#			while {[gets $f line] != -1} {
#				foreach {barcode count} [split $line \t] break
#				if {![info exists a($barcode)]} {
#					set a($barcode) $count
#				} else {
#					incr a($barcode) $count
#				}
#			}
#		}
#		# barcode counts list
#		set mergedbarcodesfile $resultdir/mergedbarcodes.tsv
#		set o [wgzopen $mergedbarcodesfile.temp]
#		puts $o barcode\tcount
#		foreach barcode [array names a] {
#			puts $o $barcode\t$a($barcode)
#		}
#		gzclose $o
		cg select -overwrite 1 -s -count $mergedbarcodesfile.temp.zst $mergedbarcodesfile.temp2.zst
		file rename -force $mergedbarcodesfile.temp2.zst $mergedbarcodesfile
		file delete $mergedbarcodesfile.temp.zst
		#
		# load whitelist
		catch {close $f}
		set f [gzopen $whitelist]
		unset -nocomplain wa
		while 1 {
			if {[gets $f barcode] == -1} break
			if {![info exists a($barcode)]} continue
			if {$a($barcode) <= 2} continue
			set wa($barcode) $a($barcode)
		}
		close $f
		#
		# match to whitelist barcodes
		# start from highest count
		# prepare matrix and list for ali_cost
		package require BioTcl
		set matrix [tempfile].mat
		file_write $matrix [ali_ident_matrix]
		set list {}
		foreach barcode [array names a] {
			if {$a($barcode) <= 2} continue
			lappend list $barcode $barcode
		}
		#

		set barcode2celbarcodefile $resultdir/barcode2celbarcode.tsv
		set f [gzopen $mergedbarcodesfile]
		set header [tsv_open $f]
		unset -nocomplain a
		while {[gets $f line] != -1} {
			foreach {barcode count} [split $line \t] break
			set a($barcode) $count
		}
		
		catch {close $f} ; catch {close $o}
		set f [gzopen $mergedbarcodesfile]
		set o [open $barcode2celbarcodefile.temp w]
		puts $o barcode\tcellbarcode\tcount\twhitelistmatch
		set header [tsv_open $f]
		# make matches with max 1 diff
		# even if a match is also in whitelist we take it as an error in the higher count barcode:
		# this is likely to come across, while it is highly unlikely to have 2 1-diff barcodes in actual cells
		# unset a($match) to indicate this one is already taken (assign to highest)
		set max {}
		# min size of barcode (start)
		set cutoff 200
		set num 0
		# maximum number of barcodes
		set maxnum 100000
		while 1 {
			puts [incr num].$line
			if {[gets $f line] == -1} break
			foreach {barcode count} [split $line \t] break
			if {![info exists wa($barcode)]} {
				# not in whitelist, but do we really want to skip these?
				continue
			}
			if {![info exists a($barcode)]} {
				# already used (as a match)
				continue
			}
			if {$cutoff eq ""} {
				if cutoff not set directly, calc as fraction from max
				set max $a($barcode)
				set cutoff [expr {$max*0.05}]
			}
			if {$count < $cutoff} break
			if {$num > $maxnum} break
			time {set matches [ali_cost -cutoff 2 -matrix $matrix $barcode $list]}
			foreach {match temp} $matches {
				if {![info exists a($match)]} continue
				if {$match eq $barcode} {
					set wmatch 2
				} elseif {[info exists wa($match)]} {
					set wmatch 1
				} else {
					set wmatch 0
				}
				puts $o $match\t$barcode\t$a($match)\t$wmatch
				unset a($match)
			}
		}
		file rename -force $barcode2celbarcodefile.temp $barcode2celbarcodefile
		# cg select -g cellbarcode barcode2celbarcode.tsv	| cg select -g all
		# cg select -g all -gc 'sum(count)' barcode2celbarcode.tsv
		cg select -g cellbarcode -gc {sum(count)} $barcode2celbarcodefile | cg select -f {cellbarcode count=$sum_count} -s -sum_count > $resultdir/reads_per_cell.tsv.temp
		file rename -force $resultdir/reads_per_cell.tsv.temp $resultdir/reads_per_cell.tsv
		unset -nocomplain a ; unset -nocomplain wa
	}

	# apply found barcodes to make barcodeinfo and bcfastq
	mkdir $resultdir/bcfastq
	mkdir $resultdir/barcodeinfo
	# add cell barcodes to per fastq files
	# set fastqs [gzfiles $fastqdir/*.fq $fastqdir/*.fastq]
	set infofiles {}
	set newfastqs {}
	foreach fastq $fastqs {
		set root [file root [gzroot [file tail $fastq]]]
		set dep2 $resultdir/barcodes/$root.barcodes.tsv
		set target $resultdir/bcfastq/$root-sc.fastq.gz
		set target2 $resultdir/barcodeinfo/$root.barcodes.tsv
		lappend bcfastqs $target
		lappend infofiles $target2
		job make_bcfastq-$sample-[file tail $fastq] -deps {
			$fastq $dep2 $resultdir/barcode2celbarcode.tsv
		} -targets {
			$target $target2
		} -vars {
			ref fastq resultdir
		} -procs {
		} -code {

			#
			# load barcodematching
			unset -nocomplain a
			catch {close $f}
			set f [open $resultdir/barcode2celbarcode.tsv]
			set header [tsv_open $f]
			if {$header ne {barcode cellbarcode count whitelistmatch}} {
				error "wrong header for file $resultdir/barcode2celbarcode.tsv"
			}
			while 1 {
				if {[gets $f line] == -1} break
				foreach {barcode cellbarcode count} [split $line \t] break
				# if {[isint $count]} {incr ca($cellbarcode) $count}
				set a($barcode) $cellbarcode
			}
			close $f
			#
			# process fastq
			catch {close $ff} ; catch {close $fb}
			catch {close $fqo} ; catch {close $fio}
			set ff [gzopen $fastq]
			set fb [gzopen $dep2]
			set header [tsv_open $fb]
			if {$header ne {id barcode umi start strand polyA}} {error "wrong header for file $dep2"}
			set fqo [wgzopen $target]
			set fio [wgzopen $target2]
			puts $fio [join {id cellbarcode oribarcode umi start strand polyA} \t]
			while 1 {
				set fbread [gets $fb line]
				set fqread [gets $ff fqname]
				if {$fqread == -1} {
					if {$fbread != -1} {
						error "mismatch between $fastq and $dep2"
					}
					break
				}
				set fqread [gets $ff fqseq]
				if {$fqread == -1} {error "incomplete error in $fastq"}
				set fqread [gets $ff fqtemp]
				if {$fqread == -1} {error "incomplete error in $fastq"}
				set fqread [gets $ff fqqual]
				if {$fqread == -1} {error "incomplete error in $fastq"}
				foreach {id barcode umi start strand polyA} [split $line \t] break
				set name [lindex $fqname 0]
				if {$barcode eq ""} {
					set oname $name
					set cellbarcode ""
				} elseif {[info exists a($barcode)]} {
					set cellbarcode $a($barcode)
					set oname "@${cellbarcode}_${umi}\#[string range $name 1 end] CB:Z:$cellbarcode CR:Z:$barcode MI:Z:$umi"
				} else {
					set oname @${barcode}_${umi}\#[string range $name 1 end]
					set celbarcode $barcode
					set oname "@${celbarcode}_${umi}\#[string range $name 1 end] CB:Z:$cellbarcode CR:Z:$barcode MI:Z:$umi"
				}
				puts $fqo $oname
				puts $fqo $fqseq
				puts $fqo $fqtemp
				puts $fqo $fqqual
				puts $fio [join [list $id $cellbarcode $umi $barcode $start $strand $polyA] \t]
			}

			close $ff ; close $fb
			close $fqo ; close $fio
		}
	}
}

proc cg_sc_barcodes args {
	set args [job_init {*}$args]
	sc_barcodes_job {*}$args
	job_wait
}
