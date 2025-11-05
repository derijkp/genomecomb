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

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;
    use tempfile::TempDir;

    fn create_test_gene_counts_file(temp_dir: &TempDir, filename: &str) -> PathBuf {
        let file_path = temp_dir.path().join(filename);
        let file = fs::File::create(&file_path).unwrap();
        let mut encoder = zstd::Encoder::new(file, 0).unwrap();
        
        // Write header
        writeln!(encoder, "col1\tcol2\tcol3\tcol4\tcol5\tcol6\tcell_barcode\tcol8").unwrap();
        
        // Cell1 has 3 genes
        writeln!(encoder, "x\tx\tx\tx\tx\tx\tCELL1\tx").unwrap();
        writeln!(encoder, "x\tx\tx\tx\tx\tx\tCELL1\tx").unwrap();
        writeln!(encoder, "x\tx\tx\tx\tx\tx\tCELL1\tx").unwrap();
        
        // Cell2 has 2 genes
        writeln!(encoder, "x\tx\tx\tx\tx\tx\tCELL2\tx").unwrap();
        writeln!(encoder, "x\tx\tx\tx\tx\tx\tCELL2\tx").unwrap();
        
        // Cell3 has 1 gene
        writeln!(encoder, "x\tx\tx\tx\tx\tx\tCELL3\tx").unwrap();
        
        encoder.finish().unwrap();
        file_path
    }

    #[test]
    fn test_get_genes_per_cell() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = create_test_gene_counts_file(&temp_dir, "test_counts.tsv.zst");
        
        let mut result = get_genes_per_cell(file_path).unwrap();
        result.sort_unstable();
        
        assert_eq!(result.len(), 3);
        assert_eq!(result, vec![1, 2, 3]);
    }

    #[test]
    fn test_plot_returns_valid_html() {
        let genes_per_cell = vec![1, 2, 3, 4, 5];
        let html = plot(genes_per_cell);
        
        assert!(html.contains("<div class=\"plot\">"));
        assert!(html.contains("<h2>Genes per cell</h2>"));
        assert!(html.contains("</div>"));
    }

    #[test]
    fn test_genes_creates_valid_table() {
        let temp_dir = TempDir::new().unwrap();
        create_test_gene_counts_file(
            &temp_dir,
            "sc_gene_counts_filtered-isoquant_sc-test.tsv.zst",
        );

        let metrics = KneeMetrics {
            num_cells: 100,
            num_good_cells: 80,
            median_umi_count: 1500.0,
            percent_reads_in_good_cells: 85.5,
        };

        let result = genes(temp_dir.path().to_str().unwrap(), metrics).unwrap();
        
        assert!(result.table.contains("<h2>Cells and genes metrics</h2>"));
        assert!(result.table.contains("Total cells"));
        assert!(result.table.contains("100"));
        assert!(result.table.contains("Cells passing filter"));
        assert!(result.table.contains("80"));
        assert!(result.table.contains("% umis in cells passing filter"));
        assert!(result.table.contains("85.50%"));
        assert!(result.table.contains("Median UMI count"));
        assert!(result.table.contains("1500"));
        assert!(result.table.contains("Median genes per cell"));
        assert!(result.table.contains("2")); // median of [1, 2, 3]
        
        assert!(result.plot.contains("<div class=\"plot\">"));
        assert!(result.plot.contains("<h2>Genes per cell</h2>"));
    }
}
