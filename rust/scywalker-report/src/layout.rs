use plotly::{
    layout::{Axis, AxisType, Margin},
    Layout,
};

pub fn specify_layout(xlabel: &str, ylabel: &str) -> Layout {
    Layout::new()
        .x_axis(Axis::new().title(xlabel))
        .y_axis(Axis::new().title(ylabel))
        .width(1200)
        .height(600)
        .margin(Margin::new().top(20).bottom(100).left(70).right(50))
}

pub fn specify_layout_loglog(xlabel: &str, ylabel: &str) -> Layout {
    Layout::new()
        .x_axis(Axis::new().title(xlabel).type_(AxisType::Log))
        .y_axis(Axis::new().title(ylabel).type_(AxisType::Log))
        .width(1200)
        .height(600)
        .margin(Margin::new().top(20).bottom(100).left(70).right(50))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_specify_layout() {
        let layout = specify_layout("X Label", "Y Label");
        // Verify the layout is created without panicking
        // The plotly Layout type doesn't expose getters, so we can't directly assert values
        // but we can verify it builds successfully
        let _ = layout;
    }

    #[test]
    fn test_specify_layout_loglog() {
        let layout = specify_layout_loglog("Log X", "Log Y");
        // Verify the layout is created without panicking
        let _ = layout;
    }

    #[test]
    fn test_layout_with_empty_labels() {
        let layout = specify_layout("", "");
        let _ = layout;
    }

    #[test]
    fn test_layout_loglog_with_special_chars() {
        let layout = specify_layout_loglog("X (units)", "Y [counts]");
        let _ = layout;
    }
}
