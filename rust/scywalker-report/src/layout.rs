use plotly::{
    common::Title,
    layout::{Axis, AxisType, Margin},
    Layout,
};

pub fn specify_layout(xlabel: &str, ylabel: &str) -> Layout {
    Layout::new()
        .x_axis(Axis::new().title(Title::new(xlabel)))
        .y_axis(Axis::new().title(Title::new(ylabel)))
        .width(1200)
        .height(600)
        .margin(Margin::new().top(20).bottom(100).left(70).right(50))
}

pub fn specify_layout_loglog(xlabel: &str, ylabel: &str) -> Layout {
    Layout::new()
        .x_axis(Axis::new().title(Title::new(xlabel)).type_(AxisType::Log))
        .y_axis(Axis::new().title(Title::new(ylabel)).type_(AxisType::Log))
        .width(1200)
        .height(600)
        .margin(Margin::new().top(20).bottom(100).left(70).right(50))
}
