// Copyright (c) Facebook, Inc. and its affiliates.
//
// This source code is licensed under the MIT license found in the
// LICENSE file in the root directory of this source tree.

use super::super::mds::mds_f64_16x16::mds_multiply;
use super::{Digest, ElementHasher, Hasher};
use core::convert::TryInto;
use core::ops::Range;
use math::{fields::f64::BaseElement, FieldElement, StarkField};

mod digest;
pub use digest::ElementDigest;

#[cfg(test)]
mod tests;

// CONSTANTS
// ================================================================================================

/// Sponge state is set to 16 field elements or 128 bytes; 10 elements are reserved for rate and
/// the remaining 6 elements are reserved for capacity.
const STATE_WIDTH: usize = 16;

/// The rate portion of the state is located in elements 0 through 9.
const RATE_RANGE: Range<usize> = 0..10;
const RATE_WIDTH: usize = RATE_RANGE.end - RATE_RANGE.start;

const INPUT1_RANGE: Range<usize> = 0..5;
const INPUT2_RANGE: Range<usize> = 5..10;

/// The capacity portion of the state is located in elements 10 through 15.
const CAPACITY_RANGE: Range<usize> = 10..16;

/// The output of the hash function is a digest which consists of 5 field elements or 40 bytes.
///
/// The digest is returned from state elements 0 through 4 (the first four elements of the
/// rate portion).
const DIGEST_RANGE: Range<usize> = 0..5;
const DIGEST_SIZE: usize = DIGEST_RANGE.end - DIGEST_RANGE.start;

/// The number of rounds is set to 5 to target 160-bit security level
const NUM_ROUNDS: usize = 5;

// HASHER IMPLEMENTATION
// ================================================================================================

/// Implementation of [Hasher] trait for Tip5 hash function with 320-bit output.
pub struct Tip5_320();

impl Hasher for Tip5_320 {
    type Digest = ElementDigest;

    const COLLISION_RESISTANCE: u32 = 160;

    fn hash(bytes: &[u8]) -> Self::Digest {
        // compute the number of elements required to represent the string; we will be processing
        // the string in 7-byte chunks, thus the number of elements will be equal to the number
        // of such chunks (including a potential partial chunk at the end).
        let num_elements = if bytes.len() % 7 == 0 {
            bytes.len() / 7
        } else {
            bytes.len() / 7 + 1
        };

        // initialize state to all zeros, except for the first element of the capacity part, which
        // is set to the number of elements to be hashed. this is done so that adding zero elements
        // at the end of the list always results in a different hash.
        let mut state = [BaseElement::ZERO; STATE_WIDTH];
        state[CAPACITY_RANGE.start] = BaseElement::new(num_elements as u64);

        // break the string into 7-byte chunks, convert each chunk into a field element, and
        // absorb the element into the rate portion of the state. we use 7-byte chunks because
        // every 7-byte chunk is guaranteed to map to some field element.
        let mut i = 0;
        let mut buf = [0_u8; 8];
        for chunk in bytes.chunks(7) {
            if i < num_elements - 1 {
                buf[..7].copy_from_slice(chunk);
            } else {
                // if we are dealing with the last chunk, it may be smaller than 7 bytes long, so
                // we need to handle it slightly differently. we also append a byte with value 1
                // to the end of the string; this pads the string in such a way that adding
                // trailing zeros results in different hash
                let chunk_len = chunk.len();
                buf = [0_u8; 8];
                buf[..chunk_len].copy_from_slice(chunk);
                buf[chunk_len] = 1;
            }

            // convert the bytes into a field element and absorb it into the rate portion of the
            // state; if the rate is filled up, apply the Tip5 permutation and start absorbing
            // again from zero index.
            state[RATE_RANGE.start + i] += BaseElement::new(u64::from_le_bytes(buf));
            i += 1;
            if i % RATE_WIDTH == 0 {
                Self::apply_permutation(&mut state);
                i = 0;
            }
        }

        // if we absorbed some elements but didn't apply a permutation to them (would happen when
        // the number of elements is not a multiple of RATE_WIDTH), apply the Tip5 permutation.
        // we don't need to apply any extra padding because we injected total number of elements
        // in the input list into the capacity portion of the state during initialization.
        if i > 0 {
            Self::apply_permutation(&mut state);
        }

        // return the first 4 elements of the state as hash result
        ElementDigest::new(state[DIGEST_RANGE].try_into().unwrap())
    }

    fn merge(values: &[Self::Digest; 2]) -> Self::Digest {
        // initialize the state by copying the digest elements into the rate portion of the state
        // (8 total elements), and set the first capacity element to 8 (the number of elements to
        // be hashed).
        let mut state = [BaseElement::ZERO; STATE_WIDTH];
        state[RATE_RANGE].copy_from_slice(Self::Digest::digests_as_elements(values));
        state[CAPACITY_RANGE.start] = BaseElement::new(RATE_WIDTH as u64);

        // apply the Rescue permutation and return the first four elements of the state
        Self::apply_permutation(&mut state);
        ElementDigest::new(state[DIGEST_RANGE].try_into().unwrap())
    }

    fn merge_with_int(seed: Self::Digest, value: u64) -> Self::Digest {
        // initialize the state as follows:
        // - seed is copied into the first 4 elements of the rate portion of the state.
        // - if the value fits into a single field element, copy it into the fifth rate element
        //   and set the first capacity element to 5 (the number of elements to be hashed).
        // - if the value doesn't fit into a single field element, split it into two field
        //   elements, copy them into rate elements 5 and 6, and set the first capacity element
        //   to 6.
        let mut state = [BaseElement::ZERO; STATE_WIDTH];
        state[INPUT1_RANGE].copy_from_slice(seed.as_elements());
        state[INPUT2_RANGE.start] = BaseElement::new(value);
        if value < BaseElement::MODULUS {
            state[CAPACITY_RANGE.start] = BaseElement::new(DIGEST_SIZE as u64 + 1);
        } else {
            state[INPUT2_RANGE.start + 1] = BaseElement::new(value / BaseElement::MODULUS);
            state[CAPACITY_RANGE.start] = BaseElement::new(DIGEST_SIZE as u64 + 2);
        }

        // apply the Rescue permutation and return the first four elements of the state
        Self::apply_permutation(&mut state);
        ElementDigest::new(state[DIGEST_RANGE].try_into().unwrap())
    }
}

impl ElementHasher for Tip5_320 {
    type BaseField = BaseElement;

    fn hash_elements<E: FieldElement<BaseField = Self::BaseField>>(elements: &[E]) -> Self::Digest {
        // convert the elements into a list of base field elements
        let elements = E::as_base_elements(elements);

        // initialize state to all zeros, except for the last element of the capacity part, which
        // is set to the number of elements to be hashed. this is done so that adding zero elements
        // at the end of the list always results in a different hash.
        let mut state = [BaseElement::ZERO; STATE_WIDTH];
        state[CAPACITY_RANGE.start] = BaseElement::new(elements.len() as u64);

        // absorb elements into the state one by one until the rate portion of the state is filled
        // up; then apply the Tip5 permutation and start absorbing again; repeat until all
        // elements have been absorbed
        let mut i = 0;
        for &element in elements.iter() {
            state[RATE_RANGE.start + i] += element;
            i += 1;
            if i % RATE_WIDTH == 0 {
                Self::apply_permutation(&mut state);
                i = 0;
            }
        }

        // if we absorbed some elements but didn't apply a permutation to them (would happen when
        // the number of elements is not a multiple of RATE_WIDTH), apply the Tip5 permutation.
        // we don't need to apply any extra padding because we injected total number of elements
        // in the input list into the capacity portion of the state during initialization.
        if i > 0 {
            Self::apply_permutation(&mut state);
        }

        // return the first 4 elements of the state as hash result
        ElementDigest::new(state[DIGEST_RANGE].try_into().unwrap())
    }
}

// HASH FUNCTION IMPLEMENTATION
// ================================================================================================

impl Tip5_320 {
    // CONSTANTS
    // --------------------------------------------------------------------------------------------

    pub const NUM_ROUNDS: usize = NUM_ROUNDS;

    pub const STATE_WIDTH: usize = STATE_WIDTH;

    pub const RATE_RANGE: Range<usize> = RATE_RANGE;

    pub const CAPACITY_RANGE: Range<usize> = CAPACITY_RANGE;

    pub const DIGEST_RANGE: Range<usize> = DIGEST_RANGE;

    pub const LOOKUP: [u8; 256] = LOOKUP_TABLE;

    pub const MDS: [[BaseElement; STATE_WIDTH]; STATE_WIDTH] = MDS;

    pub const ARK: [[BaseElement; STATE_WIDTH]; NUM_ROUNDS] = ARK;

    // Tip5 PERMUTATION
    // --------------------------------------------------------------------------------------------

    /// Applies Tip5 permutation to the provided state.
    pub fn apply_permutation(state: &mut [BaseElement; STATE_WIDTH]) {
        for i in 0..NUM_ROUNDS {
            Self::apply_round(state, i);
        }
    }

    /// Tip5 round function.
    #[inline(always)]
    pub fn apply_round(state: &mut [BaseElement; STATE_WIDTH], round: usize) {
        Self::apply_sbox(state);
        Self::apply_mds(state);
        Self::add_constants(state, &ARK[round]);
    }

    // HELPER FUNCTIONS
    // --------------------------------------------------------------------------------------------

    #[inline(always)]
    fn apply_mds(state: &mut [BaseElement; STATE_WIDTH]) {
        mds_multiply(state)
    }

    #[inline(always)]
    fn add_constants(state: &mut [BaseElement; STATE_WIDTH], ark: &[BaseElement; STATE_WIDTH]) {
        state.iter_mut().zip(ark).for_each(|(s, &k)| *s += k);
    }

    #[inline(always)]
    fn apply_sbox(state: &mut [BaseElement; STATE_WIDTH]) {
        state[0] = apply_split_and_lookup(state[0]);
        state[1] = apply_split_and_lookup(state[1]);
        state[2] = apply_split_and_lookup(state[2]);
        state[3] = apply_split_and_lookup(state[3]);
        state[4] = state[4].exp7();
        state[5] = state[5].exp7();
        state[6] = state[6].exp7();
        state[7] = state[7].exp7();
        state[8] = state[8].exp7();
        state[9] = state[9].exp7();
        state[10] = state[10].exp7();
        state[11] = state[11].exp7();
        state[12] = state[12].exp7();
        state[13] = state[13].exp7();
        state[14] = state[14].exp7();
        state[15] = state[15].exp7();
    }
}

// HELPER METHODS
// ================================================================================================

#[inline]
fn apply_split_and_lookup(x: BaseElement) -> BaseElement {
    // x is already in Montgomery form, we can directly extract its bytes
    let mut bytes = x.inner().to_le_bytes();

    for i in 0..8 {
        bytes[i] = LOOKUP_TABLE[bytes[i] as usize];
    }

    BaseElement::from_mont(u64::from_le_bytes(bytes))
}

// LOOKUP
// ================================================================================================
/// Tip5 Lookup table for the S-Box S
pub const LOOKUP_TABLE: [u8; 256] = [
    0, 7, 26, 63, 124, 215, 85, 254, 214, 228, 45, 185, 140, 173, 33, 240, 29, 177, 176, 32, 8,
    110, 87, 202, 204, 99, 150, 106, 230, 14, 235, 128, 213, 239, 212, 138, 23, 130, 208, 6, 44,
    71, 93, 116, 146, 189, 251, 81, 199, 97, 38, 28, 73, 179, 95, 84, 152, 48, 35, 119, 49, 88,
    242, 3, 148, 169, 72, 120, 62, 161, 166, 83, 175, 191, 137, 19, 100, 129, 112, 55, 221, 102,
    218, 61, 151, 237, 68, 164, 17, 147, 46, 234, 203, 216, 22, 141, 65, 57, 123, 12, 244, 54, 219,
    231, 96, 77, 180, 154, 5, 253, 133, 165, 98, 195, 205, 134, 245, 30, 9, 188, 59, 142, 186, 197,
    181, 144, 92, 31, 224, 163, 111, 74, 58, 69, 113, 196, 67, 246, 225, 10, 121, 50, 60, 157, 90,
    122, 2, 250, 101, 75, 178, 159, 24, 36, 201, 11, 243, 132, 198, 190, 114, 233, 39, 52, 21, 209,
    108, 238, 91, 187, 18, 104, 194, 37, 153, 34, 200, 143, 126, 155, 236, 118, 64, 80, 172, 89,
    94, 193, 135, 183, 86, 107, 252, 13, 167, 206, 136, 220, 207, 103, 171, 160, 76, 182, 227, 217,
    158, 56, 174, 4, 66, 109, 139, 162, 184, 211, 249, 47, 125, 232, 117, 43, 16, 42, 127, 20, 241,
    25, 149, 105, 156, 51, 53, 168, 145, 247, 223, 79, 78, 226, 15, 222, 82, 115, 70, 210, 27, 41,
    1, 170, 40, 131, 192, 229, 248, 255,
];

// MDS
// ================================================================================================
/// Tip5 MDS matrix
const MDS: [[BaseElement; STATE_WIDTH]; STATE_WIDTH] = [
    [
        BaseElement::new(18446742626305572865),
        BaseElement::new(12103477672291081763),
        BaseElement::new(2060678330275253760),
        BaseElement::new(2379906257383880660),
        BaseElement::new(13229869363167232),
        BaseElement::new(12103947580051285184),
        BaseElement::new(3631871650260640512),
        BaseElement::new(18375880622847003104),
        BaseElement::new(18446743700047396865),
        BaseElement::new(7491684177549991903),
        BaseElement::new(16328081945490043393),
        BaseElement::new(11603770461046470701),
        BaseElement::new(18433515290973110273),
        BaseElement::new(5185371269165997376),
        BaseElement::new(14872856315882446081),
        BaseElement::new(4542937756286289441),
    ],
    [
        BaseElement::new(4542937756286289441),
        BaseElement::new(18446742626305572865),
        BaseElement::new(12103477672291081763),
        BaseElement::new(2060678330275253760),
        BaseElement::new(2379906257383880660),
        BaseElement::new(13229869363167232),
        BaseElement::new(12103947580051285184),
        BaseElement::new(3631871650260640512),
        BaseElement::new(18375880622847003104),
        BaseElement::new(18446743700047396865),
        BaseElement::new(7491684177549991903),
        BaseElement::new(16328081945490043393),
        BaseElement::new(11603770461046470701),
        BaseElement::new(18433515290973110273),
        BaseElement::new(5185371269165997376),
        BaseElement::new(14872856315882446081),
    ],
    [
        BaseElement::new(14872856315882446081),
        BaseElement::new(4542937756286289441),
        BaseElement::new(18446742626305572865),
        BaseElement::new(12103477672291081763),
        BaseElement::new(2060678330275253760),
        BaseElement::new(2379906257383880660),
        BaseElement::new(13229869363167232),
        BaseElement::new(12103947580051285184),
        BaseElement::new(3631871650260640512),
        BaseElement::new(18375880622847003104),
        BaseElement::new(18446743700047396865),
        BaseElement::new(7491684177549991903),
        BaseElement::new(16328081945490043393),
        BaseElement::new(11603770461046470701),
        BaseElement::new(18433515290973110273),
        BaseElement::new(5185371269165997376),
    ],
    [
        BaseElement::new(5185371269165997376),
        BaseElement::new(14872856315882446081),
        BaseElement::new(4542937756286289441),
        BaseElement::new(18446742626305572865),
        BaseElement::new(12103477672291081763),
        BaseElement::new(2060678330275253760),
        BaseElement::new(2379906257383880660),
        BaseElement::new(13229869363167232),
        BaseElement::new(12103947580051285184),
        BaseElement::new(3631871650260640512),
        BaseElement::new(18375880622847003104),
        BaseElement::new(18446743700047396865),
        BaseElement::new(7491684177549991903),
        BaseElement::new(16328081945490043393),
        BaseElement::new(11603770461046470701),
        BaseElement::new(18433515290973110273),
    ],
    [
        BaseElement::new(18433515290973110273),
        BaseElement::new(5185371269165997376),
        BaseElement::new(14872856315882446081),
        BaseElement::new(4542937756286289441),
        BaseElement::new(18446742626305572865),
        BaseElement::new(12103477672291081763),
        BaseElement::new(2060678330275253760),
        BaseElement::new(2379906257383880660),
        BaseElement::new(13229869363167232),
        BaseElement::new(12103947580051285184),
        BaseElement::new(3631871650260640512),
        BaseElement::new(18375880622847003104),
        BaseElement::new(18446743700047396865),
        BaseElement::new(7491684177549991903),
        BaseElement::new(16328081945490043393),
        BaseElement::new(11603770461046470701),
    ],
    [
        BaseElement::new(11603770461046470701),
        BaseElement::new(18433515290973110273),
        BaseElement::new(5185371269165997376),
        BaseElement::new(14872856315882446081),
        BaseElement::new(4542937756286289441),
        BaseElement::new(18446742626305572865),
        BaseElement::new(12103477672291081763),
        BaseElement::new(2060678330275253760),
        BaseElement::new(2379906257383880660),
        BaseElement::new(13229869363167232),
        BaseElement::new(12103947580051285184),
        BaseElement::new(3631871650260640512),
        BaseElement::new(18375880622847003104),
        BaseElement::new(18446743700047396865),
        BaseElement::new(7491684177549991903),
        BaseElement::new(16328081945490043393),
    ],
    [
        BaseElement::new(16328081945490043393),
        BaseElement::new(11603770461046470701),
        BaseElement::new(18433515290973110273),
        BaseElement::new(5185371269165997376),
        BaseElement::new(14872856315882446081),
        BaseElement::new(4542937756286289441),
        BaseElement::new(18446742626305572865),
        BaseElement::new(12103477672291081763),
        BaseElement::new(2060678330275253760),
        BaseElement::new(2379906257383880660),
        BaseElement::new(13229869363167232),
        BaseElement::new(12103947580051285184),
        BaseElement::new(3631871650260640512),
        BaseElement::new(18375880622847003104),
        BaseElement::new(18446743700047396865),
        BaseElement::new(7491684177549991903),
    ],
    [
        BaseElement::new(7491684177549991903),
        BaseElement::new(16328081945490043393),
        BaseElement::new(11603770461046470701),
        BaseElement::new(18433515290973110273),
        BaseElement::new(5185371269165997376),
        BaseElement::new(14872856315882446081),
        BaseElement::new(4542937756286289441),
        BaseElement::new(18446742626305572865),
        BaseElement::new(12103477672291081763),
        BaseElement::new(2060678330275253760),
        BaseElement::new(2379906257383880660),
        BaseElement::new(13229869363167232),
        BaseElement::new(12103947580051285184),
        BaseElement::new(3631871650260640512),
        BaseElement::new(18375880622847003104),
        BaseElement::new(18446743700047396865),
    ],
    [
        BaseElement::new(18446743700047396865),
        BaseElement::new(7491684177549991903),
        BaseElement::new(16328081945490043393),
        BaseElement::new(11603770461046470701),
        BaseElement::new(18433515290973110273),
        BaseElement::new(5185371269165997376),
        BaseElement::new(14872856315882446081),
        BaseElement::new(4542937756286289441),
        BaseElement::new(18446742626305572865),
        BaseElement::new(12103477672291081763),
        BaseElement::new(2060678330275253760),
        BaseElement::new(2379906257383880660),
        BaseElement::new(13229869363167232),
        BaseElement::new(12103947580051285184),
        BaseElement::new(3631871650260640512),
        BaseElement::new(18375880622847003104),
    ],
    [
        BaseElement::new(18375880622847003104),
        BaseElement::new(18446743700047396865),
        BaseElement::new(7491684177549991903),
        BaseElement::new(16328081945490043393),
        BaseElement::new(11603770461046470701),
        BaseElement::new(18433515290973110273),
        BaseElement::new(5185371269165997376),
        BaseElement::new(14872856315882446081),
        BaseElement::new(4542937756286289441),
        BaseElement::new(18446742626305572865),
        BaseElement::new(12103477672291081763),
        BaseElement::new(2060678330275253760),
        BaseElement::new(2379906257383880660),
        BaseElement::new(13229869363167232),
        BaseElement::new(12103947580051285184),
        BaseElement::new(3631871650260640512),
    ],
    [
        BaseElement::new(3631871650260640512),
        BaseElement::new(18375880622847003104),
        BaseElement::new(18446743700047396865),
        BaseElement::new(7491684177549991903),
        BaseElement::new(16328081945490043393),
        BaseElement::new(11603770461046470701),
        BaseElement::new(18433515290973110273),
        BaseElement::new(5185371269165997376),
        BaseElement::new(14872856315882446081),
        BaseElement::new(4542937756286289441),
        BaseElement::new(18446742626305572865),
        BaseElement::new(12103477672291081763),
        BaseElement::new(2060678330275253760),
        BaseElement::new(2379906257383880660),
        BaseElement::new(13229869363167232),
        BaseElement::new(12103947580051285184),
    ],
    [
        BaseElement::new(12103947580051285184),
        BaseElement::new(3631871650260640512),
        BaseElement::new(18375880622847003104),
        BaseElement::new(18446743700047396865),
        BaseElement::new(7491684177549991903),
        BaseElement::new(16328081945490043393),
        BaseElement::new(11603770461046470701),
        BaseElement::new(18433515290973110273),
        BaseElement::new(5185371269165997376),
        BaseElement::new(14872856315882446081),
        BaseElement::new(4542937756286289441),
        BaseElement::new(18446742626305572865),
        BaseElement::new(12103477672291081763),
        BaseElement::new(2060678330275253760),
        BaseElement::new(2379906257383880660),
        BaseElement::new(13229869363167232),
    ],
    [
        BaseElement::new(13229869363167232),
        BaseElement::new(12103947580051285184),
        BaseElement::new(3631871650260640512),
        BaseElement::new(18375880622847003104),
        BaseElement::new(18446743700047396865),
        BaseElement::new(7491684177549991903),
        BaseElement::new(16328081945490043393),
        BaseElement::new(11603770461046470701),
        BaseElement::new(18433515290973110273),
        BaseElement::new(5185371269165997376),
        BaseElement::new(14872856315882446081),
        BaseElement::new(4542937756286289441),
        BaseElement::new(18446742626305572865),
        BaseElement::new(12103477672291081763),
        BaseElement::new(2060678330275253760),
        BaseElement::new(2379906257383880660),
    ],
    [
        BaseElement::new(2379906257383880660),
        BaseElement::new(13229869363167232),
        BaseElement::new(12103947580051285184),
        BaseElement::new(3631871650260640512),
        BaseElement::new(18375880622847003104),
        BaseElement::new(18446743700047396865),
        BaseElement::new(7491684177549991903),
        BaseElement::new(16328081945490043393),
        BaseElement::new(11603770461046470701),
        BaseElement::new(18433515290973110273),
        BaseElement::new(5185371269165997376),
        BaseElement::new(14872856315882446081),
        BaseElement::new(4542937756286289441),
        BaseElement::new(18446742626305572865),
        BaseElement::new(12103477672291081763),
        BaseElement::new(2060678330275253760),
    ],
    [
        BaseElement::new(2060678330275253760),
        BaseElement::new(2379906257383880660),
        BaseElement::new(13229869363167232),
        BaseElement::new(12103947580051285184),
        BaseElement::new(3631871650260640512),
        BaseElement::new(18375880622847003104),
        BaseElement::new(18446743700047396865),
        BaseElement::new(7491684177549991903),
        BaseElement::new(16328081945490043393),
        BaseElement::new(11603770461046470701),
        BaseElement::new(18433515290973110273),
        BaseElement::new(5185371269165997376),
        BaseElement::new(14872856315882446081),
        BaseElement::new(4542937756286289441),
        BaseElement::new(18446742626305572865),
        BaseElement::new(12103477672291081763),
    ],
    [
        BaseElement::new(12103477672291081763),
        BaseElement::new(2060678330275253760),
        BaseElement::new(2379906257383880660),
        BaseElement::new(13229869363167232),
        BaseElement::new(12103947580051285184),
        BaseElement::new(3631871650260640512),
        BaseElement::new(18375880622847003104),
        BaseElement::new(18446743700047396865),
        BaseElement::new(7491684177549991903),
        BaseElement::new(16328081945490043393),
        BaseElement::new(11603770461046470701),
        BaseElement::new(18433515290973110273),
        BaseElement::new(5185371269165997376),
        BaseElement::new(14872856315882446081),
        BaseElement::new(4542937756286289441),
        BaseElement::new(18446742626305572865),
    ],
];

// ROUND CONSTANTS
// ================================================================================================

/// Tip5 round constants;
const ARK: [[BaseElement; STATE_WIDTH]; NUM_ROUNDS] = [
    [
        BaseElement::new(13630775303355457758),
        BaseElement::new(16896927574093233874),
        BaseElement::new(10379449653650130495),
        BaseElement::new(1965408364413093495),
        BaseElement::new(15232538947090185111),
        BaseElement::new(15892634398091747074),
        BaseElement::new(3989134140024871768),
        BaseElement::new(2851411912127730865),
        BaseElement::new(8709136439293758776),
        BaseElement::new(3694858669662939734),
        BaseElement::new(12692440244315327141),
        BaseElement::new(10722316166358076749),
        BaseElement::new(12745429320441639448),
        BaseElement::new(17932424223723990421),
        BaseElement::new(7558102534867937463),
        BaseElement::new(15551047435855531404),
    ],
    [
        BaseElement::new(17532528648579384106),
        BaseElement::new(5216785850422679555),
        BaseElement::new(15418071332095031847),
        BaseElement::new(11921929762955146258),
        BaseElement::new(9738718993677019874),
        BaseElement::new(3464580399432997147),
        BaseElement::new(13408434769117164050),
        BaseElement::new(264428218649616431),
        BaseElement::new(4436247869008081381),
        BaseElement::new(4063129435850804221),
        BaseElement::new(2865073155741120117),
        BaseElement::new(5749834437609765994),
        BaseElement::new(6804196764189408435),
        BaseElement::new(17060469201292988508),
        BaseElement::new(9475383556737206708),
        BaseElement::new(12876344085611465020),
    ],
    [
        BaseElement::new(13835756199368269249),
        BaseElement::new(1648753455944344172),
        BaseElement::new(9836124473569258483),
        BaseElement::new(12867641597107932229),
        BaseElement::new(11254152636692960595),
        BaseElement::new(16550832737139861108),
        BaseElement::new(11861573970480733262),
        BaseElement::new(1256660473588673495),
        BaseElement::new(13879506000676455136),
        BaseElement::new(10564103842682358721),
        BaseElement::new(16142842524796397521),
        BaseElement::new(3287098591948630584),
        BaseElement::new(685911471061284805),
        BaseElement::new(5285298776918878023),
        BaseElement::new(18310953571768047354),
        BaseElement::new(3142266350630002035),
    ],
    [
        BaseElement::new(549990724933663297),
        BaseElement::new(4901984846118077401),
        BaseElement::new(11458643033696775769),
        BaseElement::new(8706785264119212710),
        BaseElement::new(12521758138015724072),
        BaseElement::new(11877914062416978196),
        BaseElement::new(11333318251134523752),
        BaseElement::new(3933899631278608623),
        BaseElement::new(16635128972021157924),
        BaseElement::new(10291337173108950450),
        BaseElement::new(4142107155024199350),
        BaseElement::new(16973934533787743537),
        BaseElement::new(11068111539125175221),
        BaseElement::new(17546769694830203606),
        BaseElement::new(5315217744825068993),
        BaseElement::new(4609594252909613081),
    ],
    [
        BaseElement::new(3350107164315270407),
        BaseElement::new(17715942834299349177),
        BaseElement::new(9600609149219873996),
        BaseElement::new(12894357635820003949),
        BaseElement::new(4597649658040514631),
        BaseElement::new(7735563950920491847),
        BaseElement::new(1663379455870887181),
        BaseElement::new(13889298103638829706),
        BaseElement::new(7375530351220884434),
        BaseElement::new(3502022433285269151),
        BaseElement::new(9231805330431056952),
        BaseElement::new(9252272755288523725),
        BaseElement::new(10014268662326746219),
        BaseElement::new(15565031632950843234),
        BaseElement::new(1209725273521819323),
        BaseElement::new(6024642864597845108),
    ],
];
