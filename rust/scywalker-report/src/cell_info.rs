use crate::knee::KneeMetrics;
use glob::glob;
use log::info;
use plotly::Plot;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

pub struct CellInfo {
    pub table: String,
    pub plot: String,
}

pub fn genes(
    directory: &str,
    metrics: KneeMetrics,
) -> Result<CellInfo, Box<dyn std::error::Error>> {
    let mut cell_info =
        String::from("<h2>Cells and genes metrics</h2><table class=\"styled-table\">");
    cell_info.push_str(&format!(
        "<tr><td>Total cells</td><td>{total_cells}</td></tr>",
        total_cells = metrics.num_cells,
    ));
    cell_info.push_str(&format!(
        "<tr><td>Cells passing filter</td><td>{good_cells}</td></tr>",
        good_cells = metrics.num_good_cells,
    ));
    cell_info.push_str(&format!(
        "<tr><td>% umis in cells passing filter</td><td>{percent_reads_in_good_cells:.2}%</td></tr>",
        percent_reads_in_good_cells = metrics.percent_reads_in_good_cells,
    ));
    cell_info.push_str(&format!(
        "<tr><td>Median UMI count</td><td>{median_umi_count}</td></tr>",
        median_umi_count = metrics.median_umi_count,
    ));
    let mut counts = glob(&format!(
        "{directory}/sc_gene_counts_filtered-isoquant_sc-*.tsv.zst"
    ))?;
    let mut genes_per_cell =
        get_genes_per_cell(counts.next().expect("Could not find gene_counts file")?)?;
    genes_per_cell.sort_unstable();
    let median_genes_per_cell = crate::utils::median(&genes_per_cell);
    cell_info.push_str(&format!(
        "<tr><td>Median genes per cell</td><td>{median_genes_per_cell}</td></tr>",
    ));

    cell_info.push_str("</table>");
    info!("Parsed cell and genes info");
    Ok(CellInfo {
        table: cell_info,
        plot: plot(genes_per_cell),
    })
}

fn get_genes_per_cell(path: PathBuf) -> Result<Vec<u32>, Box<dyn std::error::Error>> {
    let decoder = zstd::Decoder::new(std::fs::File::open(path)?)?;
    let reader = BufReader::new(decoder);
    let mut genes_per_cell = HashMap::new();
    for line in reader.lines().skip(1) {
        let line = line?;
        let mut line = line.split('\t');
        let cell_barcode = line
            .nth(6)
            .expect("Problem parsing gene counts file")
            .to_string();
        *genes_per_cell.entry(cell_barcode).or_insert(0) += 1;
    }
    let genes_per_cell = genes_per_cell.into_values().collect::<Vec<_>>();
    Ok(genes_per_cell)
}

fn plot(genes_per_cell: Vec<u32>) -> String {
    let mut plot = Plot::new();
    let hist = crate::utils::histogram_to_bar(genes_per_cell.to_vec());
    plot.add_trace(hist.name("Genes per cell"));
    let layout = crate::layout::specify_layout("Genes per cell", "Number of cells");
    plot.set_layout(layout);

    format!(
        "<div class=\"plot\"><h2>Genes per cell</h2>{}</div>",
        plot.to_inline_html(Some("genes_per_cell"))
    )
}
