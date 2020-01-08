proc version_gridss {} {
	set version ?
	catch {execjar gridss --version} version
	return $version
}

proc sv_gridss_job {args} {
	upvar job_logdir job_logdir
	set cmdline "[list cd [pwd]] \; [list cg sv_gridss {*}$args]"
	set refseq {}
	set opts {}
	set split 1
	set threads 2
	set cleanup 0
	set regmincoverage 3
	set resultfiles 0
	set skips {}
	set resultfile {}
	cg_options sv_gridss args {
		-refseq {
			set refseq $value
		}
		-split {
			set split $value
		}
		-threads - -t {
			set threads $value
		}
		-cleanup {
			set cleanup $value
		}
		-resultfiles {
			set resultfiles $value
		}
		-exome {
			# notused
		}
		-skip {
			lappend skips -skip $value
		}
		default {
			if {[regexp {^-..} $key]} {set key -$key}
			lappend opts $key $value
		}
	} {bamfile resultfile} 1 2
	set bamfile [file_absolute $bamfile]
	set refseq [refseq $refseq]
	if {$resultfile eq ""} {
		set root gridss-[file_rootname $bamfile]
		set resultfile [file dir $bamfile]/sv-$root.tsv.zst
	} else {
		set root [file_rootname $resultfile]
	}
	set resultanalysisinfo [gzroot $resultfile].analysisinfo
	set destdir [file dir $resultfile]
	# resultfiles
	set resultlist [list $resultfile $resultanalysisinfo]
	if {$resultfiles} {
		return $resultlist
	}
	# logfile
	job_logfile $destdir/sv_gridss_[file tail $resultfile] $destdir $cmdline \
		{*}[versions gridss gnusort8 zst os]
	# start
	set bwarefseq [refseq_bwa_job $refseq]
	## Produce gridss sv calls
	set bamfileindex $bamfile.[indexext $bamfile]
	set workdir $resultfile.gridssrun
	file mkdir $resultfile.gridssrun
	set vcffile $resultfile.gridssrun/results.vcf
	set blacklist [glob -nocomplain [file dir $refseq]/reg_*_wgEncodeDacMapabilityConsensusExcludable.tsv*]
	if {[llength $blacklist]} {
		lappend opts BLACKLIST=\"[lindex $blacklist 0]\"
	}
	job sv_gridss-$root.vcf {*}$skips -mem [expr {8+1*$threads}]G -cores $threads \
	-skip [list $resultfile $resultfile.analysisinfo] \
	-deps {
		$bamfile $bwarefseq $bamfileindex $bwarefseq.fai
	} -targets {
		$vcffile $vcffile.analysisinfo
	} -vars {
		resultfile workdir bamfile vcffile gridss opts bwarefseq threads root
	} -code {
		analysisinfo_write $dep $vcffile sample $root varcaller gridss varcaller_version [version gridss] varcaller_cg_version [version genomecomb]
		# -Dreference_fasta is only required for CRAM input files
		# -Dgridss.gridss.output_to_temp_file=true allows GRIDSS to continue where it left off without data errors due to truncated files
		# -Dsamjdk.create_index=true is required for multi-threaded operation
		# -Dsamjdk.use_async_io allow for async read/write on background threads which improves BAM I/O performancce
		set jar [findjar gridss GRIDSS]
		set java [findjava $jar]
		set mem [expr {8+1*$threads}]G
		set error [catch {
			exec $java -Xms1G -Xmx$mem -XX:ParallelGCThreads=1 \
				-Dreference_fasta=$bwarefseq \
				-Dsamjdk.create_index=true \
				-Dsamjdk.use_async_io_read_samtools=true \
				-Dsamjdk.use_async_io_write_samtools=true \
				-Dsamjdk.use_async_io_write_tribble=true \
				-Dgridss.gridss.output_to_temp_file=true \
			-cp $jar gridss.CallVariants \
				TMP_DIR=[scratchdir] \
				WORKING_DIR=$workdir \
				REFERENCE_SEQUENCE=$bwarefseq \
				INPUT=$bamfile \
				OUTPUT=$vcffile.temp.vcf \
				WORKER_THREADS=$threads \
				ASSEMBLY=$workdir/assembly.bam \
				{*}$opts >@ stdout 2>@ stderr
		} msg opt]
		if {$error} {
			if {$::errorCode ne "NONE"} {
				if {![regexp {gridss.CallVariants done. Elapsed time} $msg]} {
					dict unset opt -level
					return -options $opt $msg
				}
			}
		}
		file rename $vcffile.temp.vcf $vcffile
	}
	# 
	job sv-gridss-vcf2tsv-$root {*}$skips -deps {
		$resultfile.gridssrun/results.vcf
	} -targets {
		$resultfile $resultanalysisinfo
	} -vars {
		sample split resultfile
	} -code {
		analysisinfo_write $dep $target
		cg vcf2tsv -split $split -removefields {name filter AN AC AF AA ExcessHet InbreedingCoeff MLEAC MLEAF NDA RPA RU STR} $resultfile.gridssrun/results.vcf $target.temp[gzext $target]
		file rename -force -- $target.temp[gzext $target] $target
	}
	# cleanup
	return $resultlist
}

proc cg_sv_gridss {args} {
	set args [job_init {*}$args]
	sv_gridss_job {*}$args
	job_wait
}

if 0 {
#install.packages("stringr")
#source("https://bioconductor.org/biocLite.R")
#biocLite("VariantAnnotation")
#install.packages("devtools")
#library(devtools)
#install_github("PapenfussLab/StructuralVariantAnnotation")
library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(stringr)
#' Simple SV type classifier
simpleEventType <- function(gr) {
  return(ifelse(seqnames(gr) != seqnames(partner(gr)), "ITX", # inter-chromosomosal
          ifelse(gr$insLen >= abs(gr$svLen) * 0.7, "INS", # TODO: improve classification of complex events
           ifelse(strand(gr) == strand(partner(gr)), "INV",
            ifelse(xor(start(gr) < start(partner(gr)), strand(gr) == "-"), "DEL",
             "DUP")))))
}
# using the example in the GRIDSS /example directory
vcf <- readVcf("/data/genomecomb.testdata/tmp/sv-gridss/sv-gridss-dsbwa-ERR194147_30x_NA12878-chr21part.tsv.zst.gridssrun/results.vcf", "hg19")
gr <- breakpointRanges(vcf)
svtype <- simpleEventType(gr)
info(vcf)$SIMPLE_TYPE <- NA_character_
#info(vcf[gr$vcfId])$SIMPLE_TYPE <- svtype
info(vcf[gr$vcfId])$SVTYPE <- svtype
info(vcf[gr$vcfId])$SVLEN <- gr$svLen
writeVcf(vcf, "/data/genomecomb.testdata/tmp/sv-gridss/sv-gridss-dsbwa-ERR194147_30x_NA12878-chr21part.tsv.zst.gridssrun/results-annot.vcf")

}
