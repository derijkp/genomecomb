use log::info;
use plotly::common::{Marker, Mode};
use plotly::{Plot, Scatter};
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

pub struct KneePlot {
    pub html: String,
    pub metrics: KneeMetrics,
}

pub struct KneeMetrics {
    pub num_cells: usize,
    pub num_good_cells: usize,
    pub median_umi_count: f64,
    pub percent_reads_in_good_cells: f64,
}

pub fn knee_plot(directory: &str) -> Result<KneePlot, Box<dyn std::error::Error>> {
    let pattern = "sc_cellinfo_raw-isoquant_sc-sminimap2_splice-*.tsv.zst";
    let cellinfo_file = crate::utils::find_file(directory, pattern)
        .unwrap_or_else(|| panic!("Failed to find cellinfo file {pattern}"));
    let (good_umis, bad_umis) = get_umis(&cellinfo_file).expect("Failed to get umis");
    let good_umis_only = good_umis.iter().flatten().copied().collect::<Vec<u32>>();
    let bad_umis_only = bad_umis.iter().flatten().copied().collect::<Vec<u32>>();
    // for filtered cells
    let good_cells_count = good_umis_only.len();

    let median_umi_count = crate::utils::median(&good_umis_only);
    let good_umi_sum = good_umis_only.iter().sum::<u32>();

    // for unfiltered cells
    let total_cell_count = good_umis.len(); // this is the same as bad_umis.len()
    let total_umi_sum = bad_umis_only.iter().sum::<u32>() + good_umi_sum;

    // calculate the percentage of reads in good cells
    let percent_reads_in_good_cells = good_umi_sum as f64 / total_umi_sum as f64 * 100.0;

    // reduce the number of points to plot by only plotting a subset of the points, showing increasingly fewer points as the x value increases
    // except the first and last points - which are always shown
    // downsample the cell rank index to subset the points later
    let cell_rank_index = (0..total_cell_count)
        .filter(|i| {
            i < &20
                || i < &100 && i % 10 == 0
                || i < &1000 && i % 50 == 0
                || i % 500 == 0
                || i >= &(good_cells_count - 50)
        })
        .collect::<Vec<_>>();

    // from both the good and bad umis, get the umi counts
    // if the index is in the cell_rank_index
    let (good_index, good_value): (Vec<usize>, Vec<u32>) = good_umis
        .iter()
        .enumerate()
        .filter(|(i, c)| cell_rank_index.contains(i) && c.is_some())
        .map(|(i, c)| (i + 1, c.expect("Problem when iterating over good umis")))
        .collect::<Vec<_>>()
        .into_iter()
        .unzip();
    let (bad_index, bad_value): (Vec<usize>, Vec<u32>) = bad_umis
        .iter()
        .enumerate()
        .filter(|(i, c)| cell_rank_index.contains(i) && c.is_some())
        .map(|(i, c)| (i + 1, c.expect("Problem when iterating over bad umis")))
        .collect::<Vec<_>>()
        .into_iter()
        .unzip();

    let mut plot = Plot::new();
    plot.add_trace(
        Scatter::new(good_index, good_value)
            .mode(Mode::LinesMarkers)
            .marker(Marker::new().size(4).color("blue"))
            .name("UMIs per cell"),
    );
    plot.add_trace(
        Scatter::new(bad_index, bad_value)
            .mode(Mode::Markers)
            .marker(Marker::new().size(4).color("grey"))
            .name("UMIs per cell (empty)"),
    );
    let layout = crate::layout::specify_layout_loglog("Cell rank", "UMIs per cell");
    plot.set_layout(layout);
    let knee_html = plot.to_inline_html(Some("knee_plot"));
    let html = format!("<div class=\"plot\"><h2>Knee plot</h2>{knee_html}</div>");
    info!("Gathered knee plot information");
    Ok(KneePlot {
        html,
        metrics: KneeMetrics {
            num_cells: total_cell_count,
            num_good_cells: good_cells_count,
            median_umi_count,
            percent_reads_in_good_cells,
        },
    })
}

type UmiVec = Vec<Option<u32>>;

fn get_umis(file: &PathBuf) -> Result<(UmiVec, UmiVec), Box<dyn std::error::Error>> {
    info!("Reading file {:?}...", file);
    let decoder = zstd::Decoder::new(std::fs::File::open(file)?)
        .unwrap_or_else(|_| panic!("Failed to open file {:?}", file));
    let reader = BufReader::new(decoder);
    let mut good_umis = Vec::new();
    let mut bad_umis = Vec::new();
    // skip the header
    for line in reader.lines().skip(1) {
        let line = line?;
        let split = line.split('\t').collect::<Vec<_>>();
        if split[3] == "1" {
            good_umis.push(Some(split[2].parse()?));
            bad_umis.push(None);
        } else {
            bad_umis.push(Some(split[2].parse()?));
            good_umis.push(None);
        }
    }
    Ok((good_umis, bad_umis))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;
    use tempfile::TempDir;

    fn create_test_cellinfo_file(temp_dir: &TempDir, filename: &str) -> PathBuf {
        let file_path = temp_dir.path().join(filename);
        let file = fs::File::create(&file_path).unwrap();
        let mut encoder = zstd::Encoder::new(file, 0).unwrap();
        
        // Write header (4 columns minimum, index 2 is UMI count, index 3 is pass filter flag)
        writeln!(encoder, "col1\tcol2\tumi_count\tpass_filter").unwrap();
        
        // Good cells (pass_filter = 1)
        writeln!(encoder, "x\tx\t1000\t1").unwrap();
        writeln!(encoder, "x\tx\t900\t1").unwrap();
        writeln!(encoder, "x\tx\t800\t1").unwrap();
        
        // Bad cells (pass_filter = 0)
        writeln!(encoder, "x\tx\t100\t0").unwrap();
        writeln!(encoder, "x\tx\t50\t0").unwrap();
        writeln!(encoder, "x\tx\t25\t0").unwrap();
        
        encoder.finish().unwrap();
        file_path
    }

    #[test]
    fn test_get_umis() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = create_test_cellinfo_file(&temp_dir, "test_cellinfo.tsv.zst");
        
        let (good_umis, bad_umis) = get_umis(&file_path).unwrap();
        
        // Check total length
        assert_eq!(good_umis.len(), 6);
        assert_eq!(bad_umis.len(), 6);
        
        // Check good cells have values in good_umis and None in bad_umis
        assert_eq!(good_umis[0], Some(1000));
        assert_eq!(bad_umis[0], None);
        assert_eq!(good_umis[1], Some(900));
        assert_eq!(bad_umis[1], None);
        assert_eq!(good_umis[2], Some(800));
        assert_eq!(bad_umis[2], None);
        
        // Check bad cells have values in bad_umis and None in good_umis
        assert_eq!(good_umis[3], None);
        assert_eq!(bad_umis[3], Some(100));
        assert_eq!(good_umis[4], None);
        assert_eq!(bad_umis[4], Some(50));
        assert_eq!(good_umis[5], None);
        assert_eq!(bad_umis[5], Some(25));
    }

    #[test]
    fn test_knee_plot_structure() {
        let temp_dir = TempDir::new().unwrap();
        create_test_cellinfo_file(
            &temp_dir,
            "sc_cellinfo_raw-isoquant_sc-sminimap2_splice-test.tsv.zst",
        );
        
        let result = knee_plot(temp_dir.path().to_str().unwrap()).unwrap();
        
        // Check metrics
        assert_eq!(result.metrics.num_cells, 6);
        assert_eq!(result.metrics.num_good_cells, 3);
        assert_eq!(result.metrics.median_umi_count, 900.0);
        
        // Calculate expected percentage
        let good_sum = 1000 + 900 + 800;
        let total_sum = good_sum + 100 + 50 + 25;
        let expected_percent = (good_sum as f64 / total_sum as f64) * 100.0;
        assert!((result.metrics.percent_reads_in_good_cells - expected_percent).abs() < 0.01);
        
        // Check HTML output
        assert!(result.html.contains("<div class=\"plot\">"));
        assert!(result.html.contains("<h2>Knee plot</h2>"));
        assert!(result.html.contains("</div>"));
    }

    #[test]
    fn test_knee_metrics_calculation() {
        let temp_dir = TempDir::new().unwrap();
        create_test_cellinfo_file(
            &temp_dir,
            "sc_cellinfo_raw-isoquant_sc-sminimap2_splice-test.tsv.zst",
        );
        
        let result = knee_plot(temp_dir.path().to_str().unwrap()).unwrap();
        let metrics = result.metrics;
        
        // Verify all metrics are reasonable
        assert!(metrics.num_good_cells <= metrics.num_cells);
        assert!(metrics.median_umi_count > 0.0);
        assert!(metrics.percent_reads_in_good_cells >= 0.0);
        assert!(metrics.percent_reads_in_good_cells <= 100.0);
    }
}
