use criterion::{Criterion, criterion_group, criterion_main};
use std::path::PathBuf;

fn bench_window_pi(c: &mut Criterion) {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../rsomics-fixtures/vcf-window-pi/bench_100k.vcf");
    if !fixture.exists() {
        eprintln!("fixture not found: {}", fixture.display());
        return;
    }
    c.bench_function("window_pi_100k", |b| {
        b.iter(|| rsomics_vcf_window_pi::compute_window_pi(&fixture, 10000).unwrap());
    });
}

criterion_group!(benches, bench_window_pi);
criterion_main!(benches);
