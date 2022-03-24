# code to download giab data
# and run benchmarks via their tools
# This is just an easy (but a bit more structured) way to get benchmark material for genomecomb
# It is not complete (will not download whole giab, etc.), nor generic 


# Genome in a Bottle based resources and tools
# Genome in a Bottle overview: https://www.nist.gov/programs-projects/genome-bottle
# giab FAQ: https://www.nist.gov/programs-projects/faqs-genome-bottle
#
# Best Practices for Benchmarking Germline Small Variant Calls in Human Genomes
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6699627/
# giab benchmarking tools: https://github.com/ga4gh/benchmarking-tools
# https://github.com/ga4gh/benchmarking-tools/blob/master/resources/high-confidence-sets/giab.md
# An open resource for accurately benchmarking small variant and reference calls
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6500473/
# ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release
# results precisionFDA in https://precision.fda.gov/challenges/truth/results
# https://github.com/ga4gh/benchmarking-tools/blob/master/resources/high-confidence-sets/giab.md

# data from precision FDA?: https://precision.fda.gov/challenges/10
# https://precision.fda.gov/challenges/truth/results
# https://precision.fda.gov/files/file-FpZG9Jj0xbJy0QXjB7yb6fpX-1 -> requires login
# -> from google: https://googlegenomics.readthedocs.io/en/staging/use_cases/discover_public_data/precision_fda.html
# -> links at gs://genomics-public-data/precision-fda

# also check
# Benchmarking workflows to assess performance and suitability of germline variant calling pipelines in clinical diagnostic assays
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03934-3

# giabsv
# ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/
# https://github.com/spiralgenetics/truvari
# https://github.com/spiralgenetics/truvari/wiki/bench

proc cg_giab_gettruth {args} {
	set version 3.3.2
	set ref hg38
	set basedir {}
	set d sge
	cg_options giab_gettruth args {
		-r - -ref {
			if {$value ne "hg38"} {error "Only hg38 supported for now"}
			set ref $value
		}
		-d {
			set d $value
		}
	} {version basedir} 0 2
	if {$version eq "hybrid"} {
		if {$basedir eq ""} {
			set basedir ~/public/giab/truth/truth_hg38_hybrid2017-1.0
		}
		puts stderr "Making $basedir"
		mkdir $basedir
		cd $basedir
		foreach url {
			https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/hybrid/README.md
			https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/hybrid/hg38.hybrid.bed.gz
			https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/hybrid/hg38.hybrid.bed.gz.tbi
			https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/hybrid/hg38.hybrid.vcf.gz
			https://s3.eu-central-1.amazonaws.com/platinum-genomes/2017-1.0/hg38/hybrid/hg38.hybrid.vcf.gz.tbi
		} {
			exec wget -c $url 2>@ stderr >@ stdout
		}
		cg vcf2tsv hg38.hybrid.vcf.gz | cg select -s - | cg zst > var_hg38.hybrid.tsv.zst
		cg bed2tsv hg38.hybrid.bed.gz | cg select -s - | cg zst > reg_hg38.hybrid.tsv.zst
		return
	}
	if {[regexp {sv(.*)} $version temp svversion]} {
		if {$svversion eq ""} {
			set svversion 0.6
		} elseif {$svversion ne 0.6} {
			error "only sv0.6 supported"
		}
		if {$basedir eq ""} {
			set basedir ~/public/giab/truth/truth_hg38_sv$svversion
		}
		puts stderr "Making $basedir"
		mkdir $basedir
		cd $basedir
		set baseurl ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v$svversion/
		foreach file {
			*.vcf.gz
			*.vcf.gz.tbi
			*.bed
		} {
			exec wget -c $baseurl/$file 2>@ stderr >@ stdout
		}
		foreach vcf [glob *.vcf.gz] {
			cg vcf2tsv $vcf | cg select -s - | cg zst > sv_hg38_[file root [gzroot $vcf]].tsv.zst
		}
		foreach bed [glob *.bed] {
			cg bed2tsv $bed | cg select -s - | cg zst > reg_hg38_[file root [gzroot $bed]].tsv.zst
		}
		return
	}
	if {$basedir eq ""} {
		set basedir ~/public/giab/truth_hg38_giab_v$version
	}
	puts stderr "Making $basedir"
	mkdir $basedir
	cd $basedir
	# set baseurl ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp
	set baseurl ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab
	catch {exec wget -c $baseurl/release/AshkenazimTrio/HG002_NA24385_son/NISTv${version}/README_NISTv${version}.txt}
	foreach {subdir url} [list \
		pilot/HG001_NA12878_hg38 $baseurl/release/NA12878_HG001/NISTv${version}/GRCh38/HG001_GRCh38_* \
		AshkenazimTrio/HG002_NA24385_son_hg38 $baseurl/release/AshkenazimTrio/HG002_NA24385_son/NISTv${version}/GRCh38/HG002_GRCh38_* \
		AshkenazimTrio/HG003_NA24149_father_hg38 $baseurl/release/AshkenazimTrio/HG003_NA24149_father/NISTv${version}/GRCh38/HG003_GRCh38_* \
		AshkenazimTrio/HG004_NA24143_mother_hg38 $baseurl/release/AshkenazimTrio/HG004_NA24143_mother/NISTv${version}/GRCh38/HG004_GRCh38_* \
	] {
		putsvars subdir url
		mkdir $basedir/$subdir
		cd $basedir/$subdir
		if {[catch {
			exec wget -c $url 2>@ stderr >@ stdout
		} msg]} {
			puts "skipping $url: $msg"
			continue
		}
		set file [glob *.vcf.gz]
		if {[file exists [file root [gzroot $file]].tsv.zst]} {
			puts "Skipping vcf file $file: done"
			continue
		} else {
			puts "Converting vcf file $file to [file root [gzroot $file]].tsv.zst"
			cg vcf2tsv $file [file root [gzroot $file]].temp.tsv.zst
			file rename -force [file root [gzroot $file]].temp.tsv.zst [file root [gzroot $file]].tsv.zst
		}
		set file [glob *.bed]
		if {[file exists [gzroot $file].tsv.zst]} {
			puts "Skipping bed file $file: done"
			continue
		} else {
			puts "Converting bed file $file to [gzroot $file].tsv.zst"
			cg bed2tsv $file [gzroot $file].temp.tsv.zst
			file rename -force [gzroot $file].temp.tsv.zst [gzroot $file].tsv.zst
		}
	}
}

proc giab_getdata_job {args} {
	upvar job_logdir job_logdir
	set version precisionfda_v2016_04
	set basedir {}
	set threads 1
	set todo {}
	set paralleldownload 0
	set parts 200
	set align {}
	set refseq {}
	cg_options giab_gettruth args {
		-paralleldownload {set paralleldownload $value}
		-todo {set todo $value}
		-parts {set parts $value}
		-align {set align $value}
		-refseq {set refseq [refseq $value]}
		-threads {set threads $value}
	} {version basedir} 0 2
	if {$basedir eq ""} {
		set basedir $::env(HOME)/public/giab/$version
	}
	job_logfile $basedir/giab_getdata_$version $basedir
	puts stderr "Making $basedir"
	mkdir $basedir
	if {[regexp ^precisionfda $version]} {
		if {$version ne "precisionfda_v2016_04"} {
			error "error downloading precision FDA data: only version precisionfda_v2016_04 supported"
		}
		set todo {
			HG002_NA24385_son {
				https://storage.googleapis.com/genomics-public-data/precision-fda/input/HG002-NA24385-pFDA_S2_L002_R1_001.fastq.gz
				https://storage.googleapis.com/genomics-public-data/precision-fda/input/HG002-NA24385-pFDA_S2_L002_R2_001.fastq.gz
			}
			HG001_NA12878 {
				https://storage.googleapis.com/genomics-public-data/precision-fda/input/HG001-NA12878-pFDA_S1_L001_R1_001.fastq.gz
				https://storage.googleapis.com/genomics-public-data/precision-fda/input/HG001-NA12878-pFDA_S1_L001_R2_001.fastq.gz
			}
		}
		set parts 200
		# precisionFDA Truth Challenge V2: Calling variants from short- and long-reads in difficult-to-map regions
		# https://www.biorxiv.org/content/10.1101/2020.11.13.380741v4
		# cannot find download without login .. ?
	} elseif {$version eq "platinum_genomes"} {
		# not tested in this form !
		# NA12878: Utah woman, parents are NA12891 and NA12892, genetic disease (CYP2D6 mutation)
		# pilot of giab (https://jimb.stanford.edu/giab)
		# also in platinum genomes
		# - (https://emea.illumina.com/platinumgenomes.html): 
		# - 17 member CEPH pedigree 1463 fully sequenced
		# - Eberle, MA et al. (2017) A reference data set of 5.4 million phased human variants validated by genetic inheritance from sequencing a three-generation 17-member pedigree. Genome Research 27: 157-164. doi:10.1101/gr.210500.116
		set pretodo {
			ERR194147 HG001_NA12878	30x
			ERR194146	NA12877	30x
			ERR194158	NA12889	30x
			ERR194159	NA12890	30x
			ERR194160	NA12891	30x
			ERR194161	NA12892	30x
		}
		set pretodo {
			ERR194147 HG001_NA12878	30x
			ERR194160	NA12891	30x
			ERR194161	NA12892	30x
		}
		set todo {}
		foreach {runid name var} $pretodo {
			lappend todo ${runid}_${var}_$name
			lappend todo [todo \
				ftp://ftp.sra.ebi.ac.uk/vol1/fastq/[string range $runid 0 5]/${runid}/${runid}_1.fastq.gz \
				ftp://ftp.sra.ebi.ac.uk/vol1/fastq/[string range $runid 0 5]/${runid}/${runid}_2.fastq.gz \
			]
		}
		set parts 200
	} elseif {$version eq "giab_ont_ultralong"} {
		puts stderr "Making $basedir"
		mkdir $basedir
		set todo {
			HG002_NA24385_son_guppy3.4.5 {
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/Ultralong_OxfordNanopore/guppy-V3.4.5/HG002_ONT-UL_GIAB_20200204.fastq.gz
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/Ultralong_OxfordNanopore/guppy-V3.4.5/HG002_ONT-UL_GIAB_20200204.sequencing_summary.txt.gz
			}
			HG002_NA24385_son {
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/GM24385_1.fastq.gz
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/GM24385_2.fastq.gz
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/GM24385_3.fastq.gz
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/GM24385_1.sequencing_summary.txt.gz
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/GM24385_2.sequencing_summary.txt.gz
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/GM24385_3.sequencing_summary.txt.gz
			}
			HG003_NA24149_father {
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/UCSC_Ultralong_OxfordNanopore_Promethion/GM24149_1.fastq.gz
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/UCSC_Ultralong_OxfordNanopore_Promethion/GM24149_2.fastq.gz
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/UCSC_Ultralong_OxfordNanopore_Promethion/GM24149_3.fastq.gz
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/UCSC_Ultralong_OxfordNanopore_Promethion/GM24149_1.sequencing_summary.txt.gz
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/UCSC_Ultralong_OxfordNanopore_Promethion/GM24149_2.sequencing_summary.txt.gz
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/UCSC_Ultralong_OxfordNanopore_Promethion/GM24149_3.sequencing_summary.txt.gz
			}
			HG004_NA24143_mother {
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/UCSC_Ultralong_OxfordNanopore_Promethion/GM24143_1.fastq.gz
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/UCSC_Ultralong_OxfordNanopore_Promethion/GM24143_2.fastq.gz
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/UCSC_Ultralong_OxfordNanopore_Promethion/GM24143_3.fastq.gz
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/UCSC_Ultralong_OxfordNanopore_Promethion/GM24143_1.sequencing_summary.txt.gz
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/UCSC_Ultralong_OxfordNanopore_Promethion/GM24143_2.sequencing_summary.txt.gz
				ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/UCSC_Ultralong_OxfordNanopore_Promethion/GM24143_3.sequencing_summary.txt.gz
			}
		}
		set parts 1000
	} elseif {![llength $todo]} {
		error "version $version not supported by cg giab_getdata; must be one of: platinum_genomes precisionfda_v2016_04"
	}
	foreach {sample urls} $todo {
		mkdir $basedir/$sample
		cd $basedir/$sample
		mkdir fastqsplit
		set files {}
		set endtargets {}
		foreach url $urls {
			set file [file tail $url]
			set ext [file extension [gzroot $file]]
			lappend files $file
			if {$ext ni ".fastq .fq"} continue
			for {set part 1} {$part <= $parts} {incr part} {
				lappend endtargets fastqsplit/p${part}_$file
			}
		}
		if {$paralleldownload} {
			job giab_getdata-download-$version-$sample -skip $endtargets -deps {
			} -targets $files -vars {
				basedir version sample urls parts
			} -code {
				mkdir tmp
				cd tmp
				foreach url $urls {
					exec wget -c $url >@ stdout 2>@ stderr
					lappend files tmp/[file tail $url]
				}
				cd ..
				file rename {*}$files .
				file delete tmp
			}
		} else {
			set donedownload 1
			foreach file $files {
				if {![file exists $file]} {set donedownload 0}
			}
			set done 1
			foreach file $endtargets {
				if {![file exists $file]} {set done 0}
			}
			if {![llength $endtargets]} {set done 0}
			if {!$done && !$donedownload} {
				mkdir tmp
				cd tmp
				set tempfiles {}
				foreach url $urls {
					exec wget -c $url >@ stdout 2>@ stderr
					lappend tempfiles tmp/[file tail $url]
				}
				cd ..
				file rename {*}$tempfiles .
				file delete tmp
			}
		}
		cd $basedir/$sample
		foreach file $files {
			set ext [file extension [gzroot $file]]
			if {$ext ni ".fastq .fq"} continue
			puts stderr "splitting $file"
			fastq_split_job -parts $parts -threads $threads $file fastqsplit/$file
		}
		if {$align ne ""} {
			set fastqs [jobglob fastqsplit/*.fastq.gz]
			if {[llength $fastqs]} {
				map_job -method $align -paired 0 -threads 8 \
					map-sminimap2-${sample}_hg38.bam \
					$refseq \
					$sample {*}$fastqs
				bam_index_job map-sminimap2-${sample}_hg38.bam
			}
		}
	}
}

proc cg_giab_getdata {args} {
	set args [job_init {*}$args]
	giab_getdata_job {*}$args
	job_wait
}

if 0 {
set output_prefix tmp/precisionFDA/hap.py-
set sample tmp/precisionFDA/samples/HG002_NA24385_son
set truthsample tmp/precisionFDA/samples/truth_HG002_NA24385_son
set refseq /complgen/refseq/hg38
}

proc version_hap.py {} {
	set happy [exec which hap.py]
	set happy [file_resolve $happy]
	set version [file tail [file dir [file dir $happy]]]
	regsub ^happy- $version {} version
	return $version
}

proc giab_benchmark_job {args} {
	set refseq {}
	set cmdline "[list cd [pwd]] \; [list cg giab_benchmark {*}$args]"
	cg_options giab_benchmark args {
		-r - -refseq {
			set refseq $value
		}
	} {output_prefix truthsample sample} 3
	set refseq [refseq $refseq]
	set args [list $sample {*}$args]
	set refseq [gatk_refseq [refseq $refseq]]
	set truthvcf [gzfile $truthsample/var-*.vcf]
	set confidenttsv [gzfile $truthsample/sreg-*.tsv]
	set vcffiles {}
	foreach file $args {
		if {[file isdir $file]} {
			lappend vcffiles {*}[gzfiles $file/var-*.vcf]
		} else {
			lappend vcffiles $file
		}
	}
	# logfile
	# -------
	set destdir [file dir $output_prefix]
	job_logfile $destdir/hap.py-[file tail $output_prefix] $destdir $cmdline \
		{*}[versions hap.py os]
	# run benchmarks
	foreach queryvcf $vcffiles {
		set target ${output_prefix}[file_rootname $queryvcf].summary.tsv
		job [job_relfile2name happy- $queryvcf] -deps {
			$refseq $confidenttsv $truthvcf $queryvcf
		} -targets {
			$target
		} -vars {
			refseq confidenttsv truthvcf queryvcf output_prefix
		} -code {
			set confidentbed [tempfile]
			cg tsv2bed $confidenttsv $confidentbed
			set workdir [workdir ${output_prefix}[file_rootname $queryvcf]]
			set workprefix $workdir/[file tail ${output_prefix}][file_rootname $queryvcf]
			set ::env(HGREF) $refseq
			exec hap.py $truthvcf $queryvcf \
				-f $confidentbed \
				-o $workprefix \
				-r $refseq
			cg csv2tsv $workprefix.summary.csv $workprefix.summary.tsv
			foreach file [glob $workdir/*.csv*] {
				cg csv2tsv $file | cg zst > [file root [gzroot $file]].tsv.zst
			}
			foreach file [glob $workdir/*] {
				file rename $file [file join [file dir ${output_prefix}] [file tail $file]]
			}
			file delete $workdir
		}
	}
}

proc cg_giab_benchmark {args} {
	set args [job_init {*}$args]
	giab_benchmark_job {*}$args
	job_wait
}

proc cg_giab {subcmd args} {
	if {[auto_load cg_giab_$subcmd]} {
		cg_giab_$subcmd {*}$args
	} else {
		set subcmds {}
		foreach cmd [array names ::auto_index cg_giab_*] {
			lappend subcmds [string range $cmd 8 end]
		}
		error "subcmd not supported, must be one of: $subcmds"
	}
		
#	cg_options benchmark args {
#		-a - -analyses {
#			set analyses $value
#		}
#		-r - -regionfile {
#			set regionfile $value
#		}
#		-refcurve {
#			set refcurve $value
#		}
#	} {cmd} 3 3
}
