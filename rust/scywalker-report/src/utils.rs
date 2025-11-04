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
    let values = Array::from_shape_vec((values.len(), 1), values).expect("Failed to create array");
    let grid = GridBuilder::<Auto<_>>::from_array(&values)
        .expect("Problem when constructing grid for histogram")
        .build();

    let histogram = values.histogram(grid);
    let mut bin_edges = vec![];
    let bins = &histogram.grid().projections()[0];
    for index in 0..bins.len() {
        let range = bins.index(index);
        bin_edges.push(range.start);
    }
    let (hist_counts, _offset) = histogram.counts().to_owned().into_raw_vec_and_offset();
    Bar::new(bin_edges, hist_counts).name("Read length")
}

pub fn find_file(directory: &str, pattern: &str) -> Option<PathBuf> {
    let mut glob_pattern = PathBuf::from(directory);
    glob_pattern.push(pattern);
    let mut glob_iter = glob(
        glob_pattern
            .to_str()
            .expect("Failed to convert glob to str"),
    )
    .expect("Failed to read glob pattern");
    glob_iter
        .next()
        .map(|x| x.expect(format!("Failed to glob path {pattern}", pattern = pattern).as_str()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_median_odd_length() {
        let values = vec![1, 2, 3, 4, 5];
        assert_eq!(median(&values), 3.0);
    }

    #[test]
    fn test_median_even_length() {
        let values = vec![1, 2, 3, 4];
        assert_eq!(median(&values), 2.5);
    }

    #[test]
    fn test_median_single_element() {
        let values = vec![42];
        assert_eq!(median(&values), 42.0);
    }

    #[test]
    fn test_median_two_elements() {
        let values = vec![10, 20];
        assert_eq!(median(&values), 15.0);
    }

    #[test]
    fn test_median_large_values() {
        let values = vec![1000, 2000, 3000, 4000, 5000];
        assert_eq!(median(&values), 3000.0);
    }

    #[test]
    fn test_histogram_to_bar_creates_valid_bar() {
        let values = vec![1, 2, 2, 3, 3, 3, 4, 4, 4, 4];
        // Just verify it doesn't panic and creates a bar
        let _bar = histogram_to_bar(values);
        // If we get here without panicking, the test passes
    }

    #[test]
    fn test_histogram_to_bar_with_range() {
        // Test with a proper range of values
        let values = vec![1, 5, 10, 15, 20, 25, 30];
        let _bar = histogram_to_bar(values);
    }

    #[test]
    fn test_histogram_to_bar_larger_dataset() {
        // Test with a larger, more varied dataset
        let values = (1..100).collect::<Vec<u32>>();
        let _bar = histogram_to_bar(values);
    }
}
