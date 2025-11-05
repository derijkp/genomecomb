use glob::glob;
use log::info;
use plotly::{Bar, Plot};
use std::{
    collections::HashMap,
    io::{BufRead, BufReader},
};

pub fn plot(directory: &str) -> Result<String, Box<dyn std::error::Error>> {
    let classifications = get_classifications(directory)?;

    let mut plot = Plot::new();
    plot.add_trace(Bar::new(
        classifications.keys().cloned().collect::<Vec<_>>(),
        classifications.values().cloned().collect::<Vec<_>>(),
    ));
    let layout = crate::layout::specify_layout("Classification", "Reads count");
    plot.set_layout(layout);
    let classification_html = plot.to_inline_html(Some("classification"));
    let html =
        format!("<div class=\"plot\"><h2>Read classification</h2>{classification_html}</div>");
    info!("Gathered read classification information");
    Ok(html)
}

fn get_classifications(
    directory: &str,
) -> Result<HashMap<String, i32>, Box<dyn std::error::Error>> {
    let mut pathglob = glob(&format!(
        "{directory}/read_assignments-isoquant_sc-sminimap2_splice-*.tsv.zst"
    ))?;
    let path = pathglob
        .next()
        .expect("Could not find read_assignments file")?;
    let decoder = zstd::Decoder::new(std::fs::File::open(path)?)?;
    let reader = BufReader::new(decoder);
    let mut read_assignment = HashMap::new();
    // skip the header
    for line in reader.lines().skip(1) {
        let line = line?;
        let mut line = line.split('\t');
        let classification = line
            .nth(17)
            .expect("Problem parsing read assignment file")
            .to_string();
        *read_assignment.entry(classification).or_insert(0) += 1;
    }
    Ok(read_assignment)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;
    use tempfile::TempDir;

    fn create_test_read_assignments_file(temp_dir: &TempDir, filename: &str) {
        let file_path = temp_dir.path().join(filename);
        let file = fs::File::create(&file_path).unwrap();
        let mut encoder = zstd::Encoder::new(file, 0).unwrap();

        // Write header (18 columns, index 17 is classification)
        writeln!(
            encoder,
            "c1\tc2\tc3\tc4\tc5\tc6\tc7\tc8\tc9\tc10\tc11\tc12\tc13\tc14\tc15\tc16\tc17\tclassification"
        )
        .unwrap();

        // Write some classifications
        for _ in 0..5 {
            writeln!(
                encoder,
                "x\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tmapped"
            )
            .unwrap();
        }
        for _ in 0..3 {
            writeln!(
                encoder,
                "x\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tunmapped"
            )
            .unwrap();
        }
        for _ in 0..2 {
            writeln!(
                encoder,
                "x\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tx\tambiguous"
            )
            .unwrap();
        }

        encoder.finish().unwrap();
    }

    #[test]
    fn test_get_classifications() {
        let temp_dir = TempDir::new().unwrap();
        create_test_read_assignments_file(
            &temp_dir,
            "read_assignments-isoquant_sc-sminimap2_splice-test.tsv.zst",
        );

        let result = get_classifications(temp_dir.path().to_str().unwrap()).unwrap();

        assert_eq!(result.len(), 3);
        assert_eq!(result.get("mapped"), Some(&5));
        assert_eq!(result.get("unmapped"), Some(&3));
        assert_eq!(result.get("ambiguous"), Some(&2));
    }

    #[test]
    fn test_plot_creates_valid_html() {
        let temp_dir = TempDir::new().unwrap();
        create_test_read_assignments_file(
            &temp_dir,
            "read_assignments-isoquant_sc-sminimap2_splice-test.tsv.zst",
        );

        let result = plot(temp_dir.path().to_str().unwrap()).unwrap();

        assert!(result.contains("<div class=\"plot\">"));
        assert!(result.contains("<h2>Read classification</h2>"));
        assert!(result.contains("</div>"));
    }

    #[test]
    fn test_get_classifications_empty_file() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir
            .path()
            .join("read_assignments-isoquant_sc-sminimap2_splice-test.tsv.zst");
        let file = fs::File::create(&file_path).unwrap();
        let mut encoder = zstd::Encoder::new(file, 0).unwrap();
        writeln!(
            encoder,
            "c1\tc2\tc3\tc4\tc5\tc6\tc7\tc8\tc9\tc10\tc11\tc12\tc13\tc14\tc15\tc16\tc17\tclassification"
        )
        .unwrap();
        encoder.finish().unwrap();

        let result = get_classifications(temp_dir.path().to_str().unwrap()).unwrap();

        assert_eq!(result.len(), 0);
    }
}
