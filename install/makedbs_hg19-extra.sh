build=hg19
dest=/complgen/refseq/

mkdir -p ${dest}/${build}/extra
cd ${dest}/${build}/extra

# targetseq exome
wget http://tools.invitrogen.com/content/sfs/manuals/TargetSeq_exome_named_targets_hg19.bed
cg bed2sft TargetSeq_exome_named_targets_hg19.bed ureg_exome_targetseq.tsv
cg select -s 'chromosome begin end' ureg_exome_targetseq.tsv sreg_exome_targetseq.tsv
cg collapseoverlap -o reg_exome_targetseq.tsv sreg_exome_targetseq.tsv
rm TargetSeq_exome_named_targets_hg19.bed sreg_exome_targetseq.tsv ureg_exome_targetseq.tsv
# cg select -f '* {id=NR "-" $name}' reg_exome_targetseq.tsv | less
