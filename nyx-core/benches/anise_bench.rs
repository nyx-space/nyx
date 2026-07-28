use criterion::{criterion_group, criterion_main, Criterion};
use anise::prelude::Almanac;
use anise::constants::frames::*;
use hifitime::Epoch;

fn bench_transformations(c: &mut Criterion) {
    // Determine path based on running context (workspace root vs nyx-core dir)
    let path = if std::path::Path::new("data/01_planetary/de440s.bsp").exists() {
        "data/01_planetary/de440s.bsp"
    } else {
        "../data/01_planetary/de440s.bsp"
    };

    let ctx = Almanac::default().load(path).unwrap();
    let epoch = Epoch::from_gregorian_utc_at_midnight(2021, 1, 1);

    c.bench_function("single translation", |b| b.iter(|| {
        ctx.translate(EARTH_J2000, MOON_J2000, epoch, None).unwrap()
    }));

    c.bench_function("single rotation", |b| b.iter(|| {
        ctx.rotate(EARTH_J2000, MOON_J2000, epoch).unwrap()
    }));

    c.bench_function("full transformation", |b| b.iter(|| {
        ctx.transform(EARTH_J2000, MOON_J2000, epoch, None).unwrap()
    }));
}

criterion_group!(benches, bench_transformations);
criterion_main!(benches);
