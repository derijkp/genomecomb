use log::info;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub fn parse_cramino(directory: &str) -> Result<String, Box<dyn std::error::Error>> {
    // this function parses the cramino output and returns a string with the metrics, to a html table
    // would be much nicer if I could grab the json output from cramino but that is not yet implemented
    let path = Path::new(&directory).join("cramino-output.tsv");
    let cramino_output = std::fs::File::open(path)?;
    let reader = BufReader::new(cramino_output);
    let mut cramino_content =
        String::from("<h2>Alignment summary</h2><table class=\"styled-table\">");
    for line in reader.lines() {
        let line = line?;
        // ignore some fields
        if line.starts_with("Yield [Gb] (>25kb)")
            || line.starts_with("Mean coverage")
            || line.starts_with("N50")
            || line.starts_with("Creation time")
            || line.starts_with("Path")
        {
            continue;
        }
        if line.starts_with('#') || line.is_empty() {
            break;
        }
        let mut line = line.split('\t');
        let metric = line.next().expect("Could not parse line from cramino");
        let value = line.next().expect("Could not parse line from cramino");

        cramino_content.push_str(&format!("<tr><td>{}</td><td>{}</td></tr>", metric, value));
    }
    cramino_content.push_str("</table>");
    info!("Parsed cramino output");
    Ok(cramino_content)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;
    use tempfile::TempDir;

    #[test]
    fn test_parse_cramino_valid_file() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("cramino-output.tsv");
        let mut file = fs::File::create(&file_path).unwrap();
        writeln!(file, "Number of reads\t12345").unwrap();
        writeln!(file, "Total bases\t987654321").unwrap();
        writeln!(file, "Mean read length\t5432.1").unwrap();
        writeln!(file, "#End of file").unwrap();

        let result = parse_cramino(temp_dir.path().to_str().unwrap()).unwrap();
        
        assert!(result.contains("<h2>Alignment summary</h2>"));
        assert!(result.contains("<table class=\"styled-table\">"));
        assert!(result.contains("Number of reads"));
        assert!(result.contains("12345"));
        assert!(result.contains("Total bases"));
        assert!(result.contains("987654321"));
        assert!(result.contains("Mean read length"));
        assert!(result.contains("5432.1"));
        assert!(result.contains("</table>"));
    }

    #[test]
    fn test_parse_cramino_filtered_fields() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("cramino-output.tsv");
        let mut file = fs::File::create(&file_path).unwrap();
        writeln!(file, "Number of reads\t100").unwrap();
        writeln!(file, "Yield [Gb] (>25kb)\t50.5").unwrap();
        writeln!(file, "Mean coverage\t30x").unwrap();
        writeln!(file, "N50\t12345").unwrap();
        writeln!(file, "Creation time\t2024-01-01").unwrap();
        writeln!(file, "Path\t/some/path").unwrap();
        writeln!(file, "Valid metric\t999").unwrap();
        writeln!(file, "").unwrap();

        let result = parse_cramino(temp_dir.path().to_str().unwrap()).unwrap();
        
        // Should contain unfiltered metrics
        assert!(result.contains("Number of reads"));
        assert!(result.contains("100"));
        assert!(result.contains("Valid metric"));
        assert!(result.contains("999"));
        
        // Should NOT contain filtered metrics
        assert!(!result.contains("Yield [Gb] (>25kb)"));
        assert!(!result.contains("Mean coverage"));
        assert!(!result.contains("N50"));
        assert!(!result.contains("Creation time"));
        assert!(!result.contains("Path"));
    }

    #[test]
    fn test_parse_cramino_empty_file() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("cramino-output.tsv");
        fs::File::create(&file_path).unwrap();

        let result = parse_cramino(temp_dir.path().to_str().unwrap()).unwrap();
        
        assert!(result.contains("<h2>Alignment summary</h2>"));
        assert!(result.contains("<table class=\"styled-table\">"));
        assert!(result.contains("</table>"));
    }

    #[test]
    fn test_parse_cramino_missing_file() {
        let temp_dir = TempDir::new().unwrap();
        
        let result = parse_cramino(temp_dir.path().to_str().unwrap());
        
        assert!(result.is_err());
    }
}
