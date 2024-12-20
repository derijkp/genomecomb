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
	# set sam $resultfile.sam
	if {[file ext $fastq] in ".bam .cram .sam"} {
		set usefastq [tempfile].fastq.gz
		catch_exec samtools fastq -T "RG,CB,QT,MI,MM,ML,Mm,Ml" $fastq | gzip > $usefastq
		set ubams 1
	} else {
		set usefastq $fastq
		set ubams 0
	}
	catch_exec minimap2 -Y -a --secondary=no -x map-ont -t 4 -n 1 -m 1 -k 5 -w 1 -s 20 $ref $usefastq > $sam 2>@ stderr
	if {$ubams} {file delete $usefastq}
	# cg sam2tsv $sam | cg select -g chromosome
	# exec ~/dev/genomecomb/bin/sc_getbarcodes adapter $begin 16 10 < $sam > temp
	catch {close $f} ; catch {close $o}
	unset -nocomplain a
	set f [open "|cg sam2tsv $sam"]
	set o [wgzopen $resultfile]
	puts $o [join {id barcode umi start strand polyA} \t]
	set header [tsv_open $f]
	set poss [list_cor $header {chromosome begin end strand qname qstart qend cigar seq supplementary}]
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
			# if more than one hit for adapter, select "best" one
			if {[llength $todo] > 1} {
				# remove hits with lower mapping quality
				set qs [list_subindex $todo $mqpos]
				set limitq [expr {[lmath_max $qs]-3}]
				set num 0
				set keep {}
				foreach q $qs {
					if {$q >= $limitq} {lappend keep $num}
					incr num
				}
				if {[llength $keep] < [llength $todo]} {
					set todo [list_sub $todo $keep]
				}
			}
			if {[llength $todo] > 1} {
				# if still multiple left with similar mapquality, pick the inner one
				# (for e.g. when adapter/read1 also included in adapters added later, after 10x)
				# for now ignoring that we can have multiple hits on different strands
				# set starts [list_subindex $todo $qstartpos]
				if {[lindex $todo 0 $strandpos] eq "+"} {
					set line [lindex [lsort -integer -index $qstartpos $todo] end]
				} else {
					set line [lindex [lsort -integer -index $qstartpos $todo] 0]
				}
			} else {
				set line [lindex $todo 0]
			}
			foreach {chromosome begin end strand qname qstart qend cigar seq supplementary} [list_sub $line $poss] break
			if {$chromosome eq "*"} {
				puts $o [join [list $qname {} {} {} {} 0] \t]
			} else {
				set start $qend
				if {[regexp H $cigar]} {
					error "hardclipped sequence in line: [list set line $line]"
				}
				set barcode [string range $seq $start [expr {$start+$barcodesize-1}]]
				set umi [string range $seq [expr {$start+$barcodesize}] [expr {$start+$barcodesize+$umisize-1}]]
				# check Ts
				set post [string range $seq [expr {$start+$barcodesize+$umisize}] [expr {$start+$barcodesize+$umisize+14}]]
				set polya [regexp -all T $post]
				# if {$polya < 1} {
				# 	error "not enough Ts in line: [list set line $line]"
				# }
				if {![info exists a($barcode)]} {
					set a($barcode) 1
				} else {
					incr a($barcode)
				}
				puts $o [join [list $qname $barcode $umi $start $strand $polya] \t]
			}
			set todo [list]
		}
		# if {$supplementary} continue
		set prevqname $nqname
		lappend todo $nline
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
	file delete $sam
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
	upvar job_logdir job_logdir
	global appdir
	set cmdline [clean_cmdline cg sc_barcodes {*}$args]
	set whitelist {}
	set refseq {}
	# = end truseq universal adapter
	set adaptorseq CTACACGACGCTCTTCCGATCT
	set barcodesize 16
	set umisize 12
	# cutoff at 2 means that the actual maximum cost/difference is 1
	set costcutoff 2
	set bcparts 50
	set maxcells 100000
	set mincells 1000
	set cutoff 20
	set orphanprotection 1
	set skips {}
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
		-maxcost {
			set costcutoff [expr {$value+1}]
		}
		-cutoff {
			set cutoff $value
		}
		-umisize {
			set umisize $value
		}
		-bcparts {
			set bcparts $value
		}
		-orphanprotection {
			set orphanprotection $value
		}
		-skip {
			lappend skips -skip $value
		}
	} {fastqdir resultdir} 1 2 {
		cg sc_barcodes finds (cell)barcodes and umis from 10x single cell sequencing in nanopore reads.
		I will create the results in resultdir (if not specified the parent dir of the fastqdir):
		* dir barcodeinfo will contain a tsv file (for each fastq file in fastq) with info on barcode and umi per read
		* dir bcfastq will contains fastq files where cellbarcode and umi has been added in the
		      read name (as @<cellbarcode>_<umi>#originalname) as well as in the info "fields" CB CR and MI
		* files
	}
	if {$whitelist ne ""} {
		if {![file exists $whitelist]} {
			if {$whitelist in "10Xv3 v3"} {
				set whitelist $::genomecombdir/whitelists/3M-february-2018.txt.gz
			} elseif {$whitelist in "10Xv4 v4"} {
				set whitelist $::genomecombdir/whitelists/3M-3pgex-may-2023.txt.gz
			} elseif {$whitelist in "10Xp5v3 p5v3"} {
				set whitelist $::genomecombdir/whitelists/3M-5pgex-jan-2023.txt.gz
			} elseif {$whitelist in "10Xv2 v2"} {
				set whitelist $::genomecombdir/whitelists/737K-august-2016.txt.gz
			} else {
				error "given sc_whitelist file \"$whitelist\" does not exist, must be an existing file or one of: v4, p5v3, v3, v2"
			}
		}
		set usewhitelist 1
	} else {
		set usewhitelist 0
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
	job_logfile $resultdir/sc_barcodes_[file tail $resultdir] $resultdir $cmdline \
		{*}[versions minimap2]
	set fastqs [gzfiles $fastqdir/*.fq $fastqdir/*.fastq $fastqdir/*.bam $fastqdir/*.cram $fastqdir/*.sam]
	set summaries {}
	shadow_mkdir $resultdir/barcodes
	job_cleanup_add_shadow $resultdir/barcodes
	foreach fastq $fastqs {
		set root [file root [gzroot [file tail $fastq]]]
		set target $resultdir/barcodes/$root.barcodes.tsv.zst
		set target2 $resultdir/barcodes/$root.summary_barcodes.tsv
		lappend summaries $target2
		job find_barcodes-$sample-[file tail $fastq] {*}$skips -skip [list \
			$resultdir/mergedbarcodes.tsv.zst \
			$resultdir/barcode2celbarcode.tsv \
			$resultdir/reads_per_cell_raw.tsv \
			$resultdir/umis_per_cell_raw.tsv \
			$resultdir/bcfastq/$root-sc.fastq.gz \
			$resultdir/barcodeinfo/$root.barcodes.tsv.zst \
		] -deps {
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
	set barcodefiles [jobglob $resultdir/barcodes/*.barcodes.tsv.zst]
	set mergedbarcodesfile $resultdir/mergedbarcodes.tsv.zst
	job merge_barcodes-$sample {*}$skips -deps $barcodefiles -targets {
		$resultdir/mergedbarcodes.tsv.zst
		$resultdir/barcode_cutoff-info.tsv
	} -vars {
		barcodefiles resultdir whitelist barcodesize threads mergedbarcodesfile
		mincells maxcells cutoff
	} -code {

		#
		# load whitelist
		if {$whitelist ne ""} {
			catch {close $f}
			set f [gzopen $whitelist]
			unset -nocomplain wa
			while 1 {
				if {[gets $f barcode] == -1} break
				set wa($barcode) 1
			}
			close $f
			# llength [array names wa]
		}

		#
		# make mergedbarcodes
		set tempfile $resultdir/mergedbarcodes-prefile.zst
		# set tempfile [tempfile].zst
		cg zcat {*}$barcodefiles | cg select -s {barcode * umi *} | cg zst > $tempfile
		catch {close $o} ; catch {close $f}
		set o [wgzopen $mergedbarcodesfile.temp.zst]
		puts $o barcode\tcount\tumicount\twhitelist
		set f [gzopen $tempfile]
		set header [tsv_open $f]
		if {$header ne "id barcode umi start strand polyA"} {error "wrong file $file?"}
		set read [gets $f line]
		foreach {id prevbarcode prevumi} [split $line \t] break
		set count 1 ; set umicount 1
		while {1} {
			foreach {id barcode umi} [split $line \t] break
			if {$barcode ne $prevbarcode} {
				incr a($prevbarcode) $count
				if {$count > 0 && $barcode ne ""} {
					if {$whitelist eq ""} {
						set white -1
					} elseif {[info exists wa($prevbarcode)]} {
						set white 1
					} else {
						set white 0
					}
					puts $o $prevbarcode\t$count\t$umicount\t$white
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
		cg select -overwrite 1 -s -count $mergedbarcodesfile.temp.zst $mergedbarcodesfile.temp2.zst
		file rename -force $mergedbarcodesfile.temp2.zst $mergedbarcodesfile
		file delete $mergedbarcodesfile.temp.zst
		#
		# barcode_cutoff-info.tsv
		catch {close $f}
		set f [gzopen $mergedbarcodesfile]
		set header [tsv_open $f]
		set line [gets $f]
		foreach {barcode count umicount white} [split $line \t] break
		if {$barcode eq ""} {
			set line [gets $f]
			foreach {barcode count umicount white} [split $line \t] break
		}
		set max $count
		if {$cutoff ne ""} {
			#default cutoff is set at 20, not ""
			set ucutoff $cutoff
		} else {
			# if cutoff set to "", calc as fraction from max
			set ucutoff [expr {$max*0.05}]
			# go to at least 100
			if {$ucutoff > 100} {
				set ucutoff 100
			}
		}
		set numlines 1
		if {$white} {set whitenum 1} else {set whitenum 0}
		set minlen [expr {$barcodesize-2}]
		while 1 {
			if {[gets $f line] == -1} break
			foreach {barcode count umicount white} [split $line \t] break
			if {$barcode eq "" || $count <= 2 || [string length $barcode] <= $minlen} continue
			if {$white} {incr whitenum $white}
			if {$whitenum > $maxcells} break
			if {$count < $ucutoff && $whitenum > $mincells} break
			incr numlines
		}
		gzclose $f
		file_write $resultdir/barcode_cutoff-info.tsv [join [list \
			key\tvalue \
			max\t$max \
			cutoff\t$cutoff \
			ucutoff\t$ucutoff \
			numlines\t$numlines \
			whitenum\t$whitenum \
		] \n]
	}
	shadow_mkdir $resultdir/barcodes.temp
	job_cleanup_add_shadow $resultdir/barcodes.temp
	set barcode_matches {}
	for {set part 1} {$part <= $bcparts} {incr part} {
		set target $resultdir/barcodes.temp/barcode_matches-$part.tsv.zst
		lappend barcode_matches $target
		job sc_findbcmatches-barcode_matches-$sample-$part {*}$skips -skip [list \
			$resultdir/barcode2celbarcode.tsv \
			$resultdir/reads_per_cell_raw.tsv \
			$resultdir/umis_per_cell_raw.tsv \
		] \
		-deps {
			$resultdir/mergedbarcodes.tsv.zst
			$resultdir/barcode_cutoff-info.tsv
		} -targets {
			$target
		} -vars {
			part costcutoff bcparts barcodesize resultdir usewhitelist
		} -procs {
			ali_ident_matrix
		} -code {
			# get barcode list
			unset -nocomplain wa
			set f [gzopen $resultdir/mergedbarcodes.tsv.zst]
			set header [tsv_open $f]
			set minlen [expr {$barcodesize-2}]
			set list {}
			while 1 {
				if {[gets $f line] == -1} break
				foreach {barcode count umicount whitelist} [split $line \t] break
				if {$barcode eq "" || $count <= 2 || [string length $barcode] <= $minlen} continue
				if {$whitelist == 1} {set wa($barcode) 1}
				lappend list $barcode
			}
			gzclose $f
			# prepare matrix and list for ali_cost
			package require BioTcl
			set matrix [tempfile].mat
			file_write $matrix [ali_ident_matrix]
			# get part to test
			array set infoa [file_read $resultdir/barcode_cutoff-info.tsv]
			set partsize [expr {$infoa(numlines)/$bcparts}]
			set pos [expr {($part-1)*$partsize}]
			set todo [lrange $list $pos [expr {$pos+$partsize-1}]]
			# get and write matches
			set o [wgzopen $target.temp.zst]
			puts $o barcode\tmatches
			set num 0
			foreach barcode $todo {
				puts [incr num].$barcode
				if {$usewhitelist && ![info exists wa($barcode)]} continue
				set matches [ali_cost -names 0 -cutoff $costcutoff -matrix $matrix $barcode $list]
				puts $o $barcode\t$matches
			}
			gzclose $o
			file rename $target.temp.zst $target
		}
	}
	job barcode2celbarcode-$sample {*}$skips -deps [list \
		$mergedbarcodesfile \
		$resultdir/barcode_cutoff-info.tsv \
		{*}$barcode_matches \
	] -targets {
		$resultdir/barcode2celbarcode.tsv
	} -vars {
		barcode_matches resultdir whitelist barcodesize threads
		maxcells mincells mergedbarcodesfile cutoff usewhitelist orphanprotection
	} -code {

		# load counts in a (counts) and ua (umicounts)
		set f [gzopen $mergedbarcodesfile]
		set header [tsv_open $f]
		unset -nocomplain a
		unset -nocomplain ua
		unset -nocomplain wa
		while {[gets $f line] != -1} {
			foreach {barcode count umicount whitelist} [split $line \t] break
			set a($barcode) $count
			set ua($barcode) $umicount
			if {$whitelist == 1} {set wa($barcode) 1}
		}
		gzclose $f

		# matches with max 1 diff already done in barcodes.temp/barcode_matches-$part.tsv.zst
		# even if a match is also in whitelist we take it as an error barcode in the higher count barcode:
		# this is likely to come across, while it is highly unlikely to have 2 1-diff barcodes in actual cells
		# unset a($match) to indicate this one is already taken (assign to highest)
		catch {close $f} ; catch {close $o} ; catch {close $or} ; catch {close $fm}
		set o [open $resultdir/barcode2celbarcode.tsv.temp w]
		puts $o barcode\tcellbarcode\tcount\tumicount\twhitelistmatch
		# min size of barcode (start)
		set num 0
		# cutoff info
		array set infoa [file_read $resultdir/barcode_cutoff-info.tsv]
		set max $infoa(max)
		set or [open $resultdir/rejected_barcode2celbarcode.tsv.temp w]
		foreach barcode_match $barcode_matches {
			set fm [gzopen $barcode_match]
			set header [tsv_open $fm]
			while 1 {
				if {[gets $fm line] == -1} break
				foreach {barcode matches} [split $line \t] break
				if {![info exists a($barcode)]} {
					# already used (as a match)
					continue
				}
				puts "[incr num].$barcode ($a($barcode))"
				set count $a($barcode)
				if {$count < $cutoff && $num > $mincells} break
				if {$num > $maxcells} break
				set orphan 0
				if {$orphanprotection} {
					set testmax [expr {1.50*$count}]
					foreach {match temp} $matches {
						if {![info exists a($match)]} continue
						if {$a($match) > $testmax} {
							set orphan 1
							break
						}
						# highest should come first anyway
						break 
					}
				}
				foreach {match temp} $matches {
					if {![info exists a($match)]} continue
					if {![info exists wa($barcode)]} {
						if {$match eq $barcode} {
							set wmatch -1
						} elseif {[info exists wa($match)]} {
							set wmatch -2
						} else {
							set wmatch -3
						}
					} else {
						if {$match eq $barcode} {
							set wmatch 2
						} elseif {[info exists wa($match)]} {
							set wmatch 1
						} else {
							set wmatch 0
						}
					}
					if {!$orphan} {
						puts $o $match\t$barcode\t$a($match)\t$ua($match)\t$wmatch
					} else {
						puts $or $match\t$barcode\t$a($match)\t$ua($match)\t$wmatch
					}
					unset a($match)
					unset ua($match)
				}
			}
			gzclose $fm
		}
		close $or
		close $o

		file rename -force $resultdir/barcode2celbarcode.tsv.temp $resultdir/barcode2celbarcode.tsv
		file rename -force $resultdir/rejected_barcode2celbarcode.tsv.temp $resultdir/rejected_barcode2celbarcode.tsv
		# cg select -g cellbarcode barcode2celbarcode.tsv	| cg select -g all
		# cg select -g all -gc 'sum(count)' barcode2celbarcode.tsv
		shadow_delete $resultdir/barcodes.temp
	}
	job reads_per_cell-$sample -deps {
		$resultdir/barcode2celbarcode.tsv
	} -targets {
		$resultdir/reads_per_cell_raw.tsv
		$resultdir/umis_per_cell_raw.tsv
	} -vars {
		resultdir
	} -code {
		# readcounts
		cg select -g cellbarcode -gc {sum(count)} $resultdir/barcode2celbarcode.tsv | cg select -f {cellbarcode count=$sum_count} -s -sum_count > $resultdir/reads_per_cell_raw.tsv.temp
		file rename -force $resultdir/reads_per_cell_raw.tsv.temp $resultdir/reads_per_cell_raw.tsv
		# umicounts
		cg select -g cellbarcode -gc {sum(umicount)} $resultdir/barcode2celbarcode.tsv | cg select -f {cellbarcode count=$sum_umicount} -s -sum_umicount > $resultdir/umis_per_cell_raw.tsv.temp
		file rename -force $resultdir/umis_per_cell_raw.tsv.temp $resultdir/umis_per_cell_raw.tsv
		unset -nocomplain a ; unset -nocomplain wa
		# foreach barcode_match $barcode_matches {
		# file delete $barcode_match
		# }
	}

	# apply found barcodes to make barcodeinfo and bcfastq
	shadow_mkdir $resultdir/bcfastq
	job_cleanup_add_shadow $resultdir/bcfastq
	mkdir $resultdir/barcodeinfo
	# add cell barcodes to per fastq files
	# set fastqs [gzfiles $fastqdir/*.fq $fastqdir/*.fastq]
	set infofiles {}
	set newfastqs {}
	foreach fastq $fastqs {
		set root [file root [gzroot [file tail $fastq]]]
		set dep2 $resultdir/barcodes/$root.barcodes.tsv.zst
		set target $resultdir/bcfastq/$root-sc.fastq.gz
		set target2 $resultdir/barcodeinfo/$root.barcodes.tsv.zst
		lappend bcfastqs $target
		lappend infofiles $target2
		job make_bcfastq-$sample-[file tail $fastq] {*}$skips -deps {
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
			if {[lrange $header 0 2] ne {barcode cellbarcode count}} {
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
			catch {gzclose $ff} ; catch {gzclose $fb}
			catch {gzclose $fqo} ; catch {gzclose $fio}
			if {[file ext $fastq] in ".bam .cram .sam"} {
				set ff [open [list | samtools fastq -T "RG,CB,QT,MI,MM,ML,Mm,Ml" $fastq]]
				set ubams 1
			} else {
				set ff [gzopen $fastq]
				set ubams 0
			}
			set fb [gzopen $dep2]
			set header [tsv_open $fb]
			if {$header ne {id barcode umi start strand polyA}} {error "wrong header for file $dep2"}
			set fqo [wgzopen $target]
			set fio [wgzopen $target2]
			puts $fio [join {id cellbarcode umi oribarcode start strand polyA} \t]
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
					set cellbarcode $barcode
					set oname "@${cellbarcode}_${umi}\#[string range $name 1 end] CB:Z:$cellbarcode CR:Z:$barcode MI:Z:$umi"
				}
				puts $fqo $oname
				puts $fqo $fqseq
				puts $fqo $fqtemp
				puts $fqo $fqqual
				puts $fio [join [list $id $cellbarcode $umi $barcode $start $strand $polyA] \t]
			}

			if {$ubams} {
				gzclosesamtools $ff
			} else {
				gzclose $ff
			}
			gzclose $fb
			gzclose $fqo ; gzclose $fio
		}
	}
}

proc cg_sc_barcodes args {
	set args [job_init {*}$args]
	sc_barcodes_job {*}$args
	job_wait
}
