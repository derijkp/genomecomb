# make data for tests/process_small.tcl
# =====================================

# this is not complete
# assumes the test data will be in ~/genomecomb.smalltestdata
cd ~/genomecomb.smalltestdata/ori

# truthset
# --------
# https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/NIST_HG002_v5.0q_variant-benchmarksets_README.md

mkdir truth_hg38_v5.0q
cd truth_hg38_v5.0q

wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_smvar.benchmark.bed

wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_stvar.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_stvar.vcf.gz.tbi
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_GRCh38_v5.0q_stvar.benchmark.bed

cg vcf2tsv HG002_GRCh38_v5.0q_smvar.vcf.gz var-truth_HG002_hg38.tsv.zst
cg bed2tsv HG002_GRCh38_v5.0q_smvar.benchmark.bed sreg-truth_HG002_hg38.tsv.zst
cg vcf2tsv HG002_GRCh38_v5.0q_stvar.vcf.gz sv-svtruth_HG002_hg38.tsv.zst
cg bed2tsv HG002_GRCh38_v5.0q_stvar.benchmark.bed sreg-svtruth_HG002_hg38.tsv.zst

cg select -overwrite 1 -f '* {size=if($type eq "del",$end - $begin,length($alt))}' \
	-q 'region("chr21:41361957-41409417") and $size < 50' \
	var-truth_HG002_hg38.tsv.zst var-mx2-truth_HG002_hg38.tsv.zst
cg select -overwrite 1 -q 'region("chr21:41361957-41409417")' sreg-truth_HG002_hg38.tsv.zst sreg-mx2-truth_HG002_hg38.tsv.zst
cg select -overwrite 1 -f '* {size=if($type eq "del",$end - $begin,length($alt))}' \
	-q 'region("chr21:41361957-41409417") and $size >= 50' \
	sv-svtruth_HG002_hg38.tsv.zst sv-mx2-svtruth_HG002_hg38.tsv.zst
cg select -overwrite 1 -q 'region("chr21:41361957-41409417")' sreg-svtruth_HG002_hg38.tsv.zst sreg-mx2-svtruth_HG002_hg38.tsv.zst

cd ..

# ont data 
# --------
# from https://epi2me.nanoporetech.com/giab-2025.01/
mkdir ont_mx2
rm calls.sorted.bam.bai
CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt samtools view -h -b \
	https://42basepairs.com/download/s3/ont-open-data/giab_2025.01/basecalling/sup/HG002/PAW71238/calls.sorted.bam \
	chr21:41361957-41409417 \
	> ont_mx2/ont_HG002_mx2.bam
samtools index ont_mx2/ont_HG002_mx2.bam
samtools collate -l 1 -@ 4 --no-PG -O ont_mx2/ont_HG002_mx2.bam | samtools reset --output-fmt BAM,level=1 > ont_mx2/unaligned_ont_HG002_mx2.bam

rm calls.sorted.bam.bai
CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt samtools view -h -b \
	https://42basepairs.com/download/s3/ont-open-data/giab_2025.01/basecalling/sup/HG003/PAY87794/calls.sorted.bam \
	chr21:41361957-41409417 \
	> ont_mx2/ont_HG003_mx2.bam
samtools index ont_mx2/ont_HG003_mx2.bam
samtools collate -l 1 -@ 4 --no-PG -O ont_mx2/ont_HG003_mx2.bam | samtools reset --output-fmt BAM,level=1 > ont_mx2/unaligned_ont_HG003_mx2.bam

rm calls.sorted.bam.bai
CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt samtools view -h -b \
	https://42basepairs.com/download/s3/ont-open-data/giab_2025.01/basecalling/sup/HG004/PAY88428/calls.sorted.bam \
	chr21:41361957-41409417 \
	> ont_mx2/ont_HG004_mx2.bam
samtools index ont_mx2/ont_HG004_mx2.bam
samtools collate -l 1 -@ 4 --no-PG -O ont_mx2/ont_HG004_mx2.bam | samtools reset --output-fmt BAM,level=1 > ont_mx2/unaligned_ont_HG004_mx2.bam

# pacbio data 
# -----------
mkdir pacbio_mx2
CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt samtools view -h -b \
	https://downloads.pacbcloud.com/public/2026Q2/HG002-SPRQ-Nx/Use1/analysis/HG002.GRCh38.haplotagged.bam \
	chr21:41361957-41409417 \
	> pacbio_mx2/pacbio_HG002_mx2.bam
samtools index pacbio_mx2/pacbio_HG002_mx2.bam
samtools collate -l 1 -@ 4 --no-PG -O pacbio_mx2/pacbio_HG002_mx2.bam | samtools reset --output-fmt BAM,level=1 > pacbio_mx2/unaligned_pacbio_HG002_mx2.bam
