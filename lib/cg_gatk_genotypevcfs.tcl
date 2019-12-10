proc cg_gatk_genotypevcfs args {
	upvar job_logdir job_logdir
	set distrreg chr
	set dbdir {}
	set usecombinegvcfs 0
	set batchsize 50
	set maxmem 10
	set vqsr {}
	set newqual true
	set gatkres {}
	set cmdline [list cg gatk_genotypevcfs {*}$args]
	set opts {}
	cg_options gatk_genotypevcfs args {
		-dbdir {
			set dbdir $value
		}
		-newqual {
			if {[true $value]} {
				set newqual true
			} else {
				set newqual false
			}
		}
		default {
			lappend opts $key $value
		}
	} {gvcf vcf} 2 2 {
		extract variants from a gvcf into a vcf (using GenotypeGVCFs)
	}
	set gvcf [file_absolute $gvcf]
	set vcf [file_absolute $vcf]
	set refseq [lindex [glob $dbdir/genome_*.ifas] 0]
	set gatkrefseq [gatk_refseq_job $refseq]
	if {[file extension $gvcf] eq ".gz" && ![file exists $gvcf.tbi]} {
		puts stderr "Making index for $gvcf"
		exec tabix -p vcf $gvcf
	}
	gatkexec [list -XX:ParallelGCThreads=1 -d64 -Xms${maxmem}g -Xmx${maxmem}g -Djava.io.tmpdir=[scratchdir]] GenotypeGVCFs \
		-R $gatkrefseq \
		-V $gvcf \
		-O $vcf.temp[gzext $vcf] \
		-G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
		--verbosity ERROR \
		-new-qual $newqual \
		--create-output-variant-index true \
		{*}$opts >@ stdout 2>@stderr
	file rename -- $vcf.temp[gzext $vcf] $vcf
}
