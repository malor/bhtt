use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

use bhtt::Histogram;

const HISTOGRAM_SIZES: [usize; 6] = [8, 16, 32, 64, 128, 256];

fn insert(c: &mut Criterion) {
    let dataset = utilities::Dataset::from_file("utilities/testdata/pings.txt").unwrap();

    let mut group = c.benchmark_group("update_histogram_of_size_X_10000_times");
    for size in HISTOGRAM_SIZES.iter() {
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, &size| {
            b.iter(|| {
                let mut h = Histogram::new(size);
                for v in dataset.values() {
                    h.insert(black_box(*v));
                }
            });
        });
    }
    group.finish();
}

fn from_iter(c: &mut Criterion) {
    let dataset = utilities::Dataset::from_file("utilities/testdata/pings.txt").unwrap();

    let mut group = c.benchmark_group("create_histogram_of_size_X_from_10000_values");
    for size in HISTOGRAM_SIZES.iter() {
        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, &size| {
            b.iter(|| black_box(Histogram::from_iter(size, dataset.values())));
        });
    }
    group.finish();
}

criterion_group!(benches, insert, from_iter);
criterion_main!(benches);
