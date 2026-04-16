use ehpx::curtis;

fn main() {
    let max_stem: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(12);

    curtis::run_curtis(max_stem);
}
