use std::path::PathBuf;

use ndarray::Array;
use ndarray_stats::{
    histogram::{strategies::Auto, GridBuilder},
    HistogramExt,
};

use glob::glob;
use num::Integer;
use plotly::Bar;
use serde::Serialize;

pub fn median(values: &[u32]) -> f64 {
    // values are assumed to be sorted
    let mid = values.len() / 2;
    if values.len() % 2 == 0 {
        (values[mid] + values[mid - 1]) as f64 / 2.0
    } else {
        values[mid] as f64
    }
}

pub fn histogram_to_bar<T>(values: Vec<T>) -> Box<plotly::Bar<T, usize>>
where
    T: Integer + Serialize + Clone + num_traits::FromPrimitive + std::fmt::Debug,
{
    let values = Array::from_shape_vec((values.len(), 1), values).unwrap();
    let grid = GridBuilder::<Auto<_>>::from_array(&values).unwrap().build();

    let histogram = values.histogram(grid);
    let mut bin_edges = vec![];
    let bins = &histogram.grid().projections()[0];
    for index in 0..bins.len() {
        let range = bins.index(index);
        bin_edges.push(range.start);
    }
    let hist_counts = histogram.counts().to_owned().into_raw_vec();
    Bar::new(bin_edges, hist_counts).name("Read length")
}

pub fn find_file(directory: &str, pattern: &str) -> Result<PathBuf, Box<dyn std::error::Error>> {
    let mut glob_pattern = PathBuf::from(directory);
    glob_pattern.push(pattern);
    let mut glob_iter = glob(glob_pattern.to_str().unwrap())?;
    let file_path = glob_iter
        .next()
        .unwrap_or_else(|| panic!("Could not find file {pattern}"))?;
    Ok(file_path)
}
