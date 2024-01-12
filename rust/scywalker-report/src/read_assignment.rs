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
