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
