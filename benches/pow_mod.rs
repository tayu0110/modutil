use criterion::{black_box, criterion_group, criterion_main, Criterion};
use modutil::pow_mod;

const EXP: u64 = 0xCCCCCCCCCCCCCCCC;
const BASE: u64 = 1000000000000000000;
const ITER: u64 = 1000;

#[inline(never)]
const fn pow_mod_naive(a: u64, mut exp: u64, m: u64) -> u64 {
    let (mut res, mut val) = (1u128, a as u128);
    while exp > 0 {
        if exp & 1 != 0 {
            res = res * val % m as u128;
        }
        val = val * val % m as u128;
        exp >>= 1;
    }
    res as u64
}

fn pow_mod_bench(c: &mut Criterion) {
    c.bench_function("modutil::pow_mod", |b| {
        b.iter(|| {
            for m in black_box(1..ITER) {
                for a in 0..m {
                    pow_mod(black_box(a), black_box(EXP), black_box(m));
                }
            }
            for m in black_box(BASE..BASE + ITER) {
                for a in BASE..m {
                    pow_mod(black_box(a), black_box(EXP), black_box(m));
                }
            }
        })
    });
}

fn pow_mod_naive_bench(c: &mut Criterion) {
    c.bench_function("pow_mod_naive", |b| {
        b.iter(|| {
            for m in black_box(1..ITER) {
                for a in 0..m {
                    pow_mod_naive(black_box(a), black_box(EXP), black_box(m));
                }
            }
            for m in black_box(BASE..BASE + ITER) {
                for a in BASE..m {
                    pow_mod_naive(black_box(a), black_box(EXP), black_box(m));
                }
            }
        })
    });
}

criterion_group!(benches, pow_mod_bench, pow_mod_naive_bench);
criterion_main!(benches);
