use glob::glob;
use log::info;
use plotly::box_plot::BoxPoints;
use plotly::{BoxPlot, Plot};
pub struct Lengths {
    pub exon_vs_length: String,
    pub length: String,
}

pub fn plot(directory: &str) -> Result<Lengths, Box<dyn std::error::Error>> {
    let mut bam = glob(&format!("{directory}/map*am"))?;
    let path = bam.next().expect("Could not find bam file")?;
    // only returning 1000 reads from the bam file to make the plot lighter
    // I would prefer to use a violin plot, though
    // but those are currently not implemented in Rust plotly
    let (lengths, exons) = crate::parse_bam::extract(path, 4, 1000)?;
    let length_plot = length_plot(&lengths);
    // filter the data to remove outliers
    let zipped = std::iter::zip(lengths, exons)
        .filter(|(length, exons)| (exons < &10) && (length < &5000))
        .collect::<Vec<_>>();
    let (lengths, exons) = zipped.into_iter().unzip();
    let mut plot = Plot::new();
    plot.add_trace(
        BoxPlot::new_xy(exons, lengths)
            .box_points(BoxPoints::Outliers)
            .jitter(0.3)
            .name("Read length vs number of exons"),
    );
    let layout = crate::layout::specify_layout("Number of exons", "Read length");
    plot.set_layout(layout);
    let splice_html = plot.to_inline_html(Some("splice_plot"));
    info!("Gathered read length information");
    Ok(Lengths {
        exon_vs_length: format!(
            "<div class=\"plot\"><h2>Read length vs number of exons</h2>{splice_html}</div>"
        ),
        length: format!("<div class=\"plot\"><h2>Read lengths</h2>{length_plot}</div>"),
    })
}

fn length_plot(lengths: &[u64]) -> String {
    let mut plot = Plot::new();
    let hist = crate::utils::histogram_to_bar(lengths.to_vec());

    plot.add_trace(hist);

    let layout = crate::layout::specify_layout("Read length", "Number of reads");
    plot.set_layout(layout);
    plot.to_inline_html(Some("length_plot"))
}
