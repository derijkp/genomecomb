use log::info;
use rand::seq::IteratorRandom;
use rust_htslib::bam::record::Cigar;
use rust_htslib::{bam, bam::Read, htslib};
use std::path::PathBuf;

pub fn extract(
    input: PathBuf,
    threads: usize,
    num_reads: usize,
) -> Result<(Vec<u64>, Vec<usize>), Box<dyn std::error::Error>> {
    let mut lengths = vec![];
    let mut exons = vec![];
    let mut bam = bam::Reader::from_path(&input)
        .expect("Error opening BAM/CRAM file.\nIs the input file correct?\n\n\n\n");
    bam.set_threads(threads)
        .expect("Failure setting decompression threads");
    for read in bam
        .rc_records()
        .choose_multiple(&mut rand::thread_rng(), num_reads)
        .into_iter()
        .map(|r| r.expect("Failure parsing Bam file"))
        .filter(|record| record.flags() & (htslib::BAM_FUNMAP | htslib::BAM_FSECONDARY) as u16 == 0)
    {
        lengths.push(read.seq_len() as u64);
        exons.push(get_exon_number(&read));
    }
    info!("Extracted read lengths and exon numbers from bam file");
    Ok((lengths, exons))
}

fn get_exon_number(record: &bam::Record) -> usize {
    let mut exon_count = 1;

    for op in record.cigar().iter() {
        if let Cigar::RefSkip(_len) = op {
            exon_count += 1;
        }
    }

    exon_count
}
