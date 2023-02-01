// Copyright (c) Facebook, Inc. and its affiliates.
//
// This source code is licensed under the MIT license found in the
// LICENSE file in the root directory of this source tree.

use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use math::fields::f64::BaseElement;
use rand_utils::{rand_array, rand_value};
use winter_crypto::{
    hashers::{Tip4_256, Tip4p_256, Tip5_320},
    ElementHasher, Hasher,
};

fn tip5_320(c: &mut Criterion) {
    c.bench_function("tip5_320 merge", |b| {
        b.iter_batched(
            || {
                [
                    Tip5_320::hash(&rand_value::<u64>().to_le_bytes()),
                    Tip5_320::hash(&rand_value::<u64>().to_le_bytes()),
                ]
            },
            |state| Tip5_320::merge(&state),
            BatchSize::SmallInput,
        )
    });

    c.bench_function("tip5_320 hash 10 elements", |b| {
        b.iter_batched(
            || rand_array::<BaseElement, 10>(),
            |state| Tip5_320::hash_elements(&state),
            BatchSize::SmallInput,
        )
    });
}

fn tip4_256(c: &mut Criterion) {
    c.bench_function("tip4_256 merge", |b| {
        b.iter_batched(
            || {
                [
                    Tip4_256::hash(&rand_value::<u64>().to_le_bytes()),
                    Tip4_256::hash(&rand_value::<u64>().to_le_bytes()),
                ]
            },
            |state| Tip4_256::merge(&state),
            BatchSize::SmallInput,
        )
    });

    c.bench_function("tip4_256 merge (4 digests with Jive)", |b| {
        b.iter_batched(
            || {
                [
                    Tip4_256::hash(&rand_value::<u64>().to_le_bytes()),
                    Tip4_256::hash(&rand_value::<u64>().to_le_bytes()),
                    Tip4_256::hash(&rand_value::<u64>().to_le_bytes()),
                    Tip4_256::hash(&rand_value::<u64>().to_le_bytes()),
                ]
            },
            |state| Tip4_256::merge_4_digests(&state),
            BatchSize::SmallInput,
        )
    });

    c.bench_function("tip4_256 hash 12 elements", |b| {
        b.iter_batched(
            || rand_array::<BaseElement, 12>(),
            |state| Tip4_256::hash_elements(&state),
            BatchSize::SmallInput,
        )
    });
}

fn tip4p_256(c: &mut Criterion) {
    c.bench_function("tip4_256p merge", |b| {
        b.iter_batched(
            || {
                [
                    Tip4p_256::hash(&rand_value::<u64>().to_le_bytes()),
                    Tip4p_256::hash(&rand_value::<u64>().to_le_bytes()),
                ]
            },
            |state| Tip4p_256::merge(&state),
            BatchSize::SmallInput,
        )
    });

    c.bench_function("tip4p_256 hash 8 elements", |b| {
        b.iter_batched(
            || rand_array::<BaseElement, 8>(),
            |state| Tip4p_256::hash_elements(&state),
            BatchSize::SmallInput,
        )
    });
}

criterion_group!(tip_hash_group, tip5_320, tip4_256, tip4p_256,);
criterion_main!(tip_hash_group);
