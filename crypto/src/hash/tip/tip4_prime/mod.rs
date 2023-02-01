// Copyright (c) Facebook, Inc. and its affiliates.
//
// This source code is licensed under the MIT license found in the
// LICENSE file in the root directory of this source tree.

use super::super::mds::mds_f64_12x12::mds_multiply;
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

/// Sponge state is set to 12 field elements or 96 bytes; 8 elements are reserved for rate and
/// the remaining 4 elements are reserved for capacity.
const STATE_WIDTH: usize = 12;

/// The rate portion of the state is located in elements 0 through 7.
const RATE_RANGE: Range<usize> = 0..8;
const RATE_WIDTH: usize = RATE_RANGE.end - RATE_RANGE.start;

const INPUT1_RANGE: Range<usize> = 0..4;
const INPUT2_RANGE: Range<usize> = 4..8;

/// The capacity portion of the state is located in elements 8 through 11.
const CAPACITY_RANGE: Range<usize> = 8..12;

/// The output of the hash function is a digest which consists of 4 field elements or 32 bytes.
///
/// The digest is returned from state elements 0 through 3 (the first four elements of the
/// rate portion).
const DIGEST_RANGE: Range<usize> = 0..4;
const DIGEST_SIZE: usize = DIGEST_RANGE.end - DIGEST_RANGE.start;

/// The number of rounds is set to 5 to target 128-bit security level
const NUM_ROUNDS: usize = 5;

// HASHER IMPLEMENTATION
// ================================================================================================

/// Implementation of [Hasher] trait for Tip4' hash function with 256-bit output.
#[allow(non_camel_case_types)]
pub struct Tip4p_256();

impl Hasher for Tip4p_256 {
    type Digest = ElementDigest;

    const COLLISION_RESISTANCE: u32 = 128;

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
            // state; if the rate is filled up, apply the Tip4' permutation and start absorbing
            // again from zero index.
            state[RATE_RANGE.start + i] += BaseElement::new(u64::from_le_bytes(buf));
            i += 1;
            if i % RATE_WIDTH == 0 {
                Self::apply_permutation(&mut state);
                i = 0;
            }
        }

        // if we absorbed some elements but didn't apply a permutation to them (would happen when
        // the number of elements is not a multiple of RATE_WIDTH), apply the Tip4' permutation.
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

impl ElementHasher for Tip4p_256 {
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
        // up; then apply the Tip4' permutation and start absorbing again; repeat until all
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
        // the number of elements is not a multiple of RATE_WIDTH), apply the Tip4' permutation.
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

impl Tip4p_256 {
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

    // TIP4' PERMUTATION
    // --------------------------------------------------------------------------------------------

    /// Applies Tip4' permutation to the provided state.
    pub fn apply_permutation(state: &mut [BaseElement; STATE_WIDTH]) {
        for i in 0..NUM_ROUNDS {
            Self::apply_round(state, i);
        }
    }

    /// Tip4' round function.
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
/// Tip4' Lookup table for the S-Box S
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
/// Tip4' MDS matrix
const MDS: [[BaseElement; STATE_WIDTH]; STATE_WIDTH] = [
    [
        BaseElement::new(7),
        BaseElement::new(23),
        BaseElement::new(8),
        BaseElement::new(26),
        BaseElement::new(13),
        BaseElement::new(10),
        BaseElement::new(9),
        BaseElement::new(7),
        BaseElement::new(6),
        BaseElement::new(22),
        BaseElement::new(21),
        BaseElement::new(8),
    ],
    [
        BaseElement::new(8),
        BaseElement::new(7),
        BaseElement::new(23),
        BaseElement::new(8),
        BaseElement::new(26),
        BaseElement::new(13),
        BaseElement::new(10),
        BaseElement::new(9),
        BaseElement::new(7),
        BaseElement::new(6),
        BaseElement::new(22),
        BaseElement::new(21),
    ],
    [
        BaseElement::new(21),
        BaseElement::new(8),
        BaseElement::new(7),
        BaseElement::new(23),
        BaseElement::new(8),
        BaseElement::new(26),
        BaseElement::new(13),
        BaseElement::new(10),
        BaseElement::new(9),
        BaseElement::new(7),
        BaseElement::new(6),
        BaseElement::new(22),
    ],
    [
        BaseElement::new(22),
        BaseElement::new(21),
        BaseElement::new(8),
        BaseElement::new(7),
        BaseElement::new(23),
        BaseElement::new(8),
        BaseElement::new(26),
        BaseElement::new(13),
        BaseElement::new(10),
        BaseElement::new(9),
        BaseElement::new(7),
        BaseElement::new(6),
    ],
    [
        BaseElement::new(6),
        BaseElement::new(22),
        BaseElement::new(21),
        BaseElement::new(8),
        BaseElement::new(7),
        BaseElement::new(23),
        BaseElement::new(8),
        BaseElement::new(26),
        BaseElement::new(13),
        BaseElement::new(10),
        BaseElement::new(9),
        BaseElement::new(7),
    ],
    [
        BaseElement::new(7),
        BaseElement::new(6),
        BaseElement::new(22),
        BaseElement::new(21),
        BaseElement::new(8),
        BaseElement::new(7),
        BaseElement::new(23),
        BaseElement::new(8),
        BaseElement::new(26),
        BaseElement::new(13),
        BaseElement::new(10),
        BaseElement::new(9),
    ],
    [
        BaseElement::new(9),
        BaseElement::new(7),
        BaseElement::new(6),
        BaseElement::new(22),
        BaseElement::new(21),
        BaseElement::new(8),
        BaseElement::new(7),
        BaseElement::new(23),
        BaseElement::new(8),
        BaseElement::new(26),
        BaseElement::new(13),
        BaseElement::new(10),
    ],
    [
        BaseElement::new(10),
        BaseElement::new(9),
        BaseElement::new(7),
        BaseElement::new(6),
        BaseElement::new(22),
        BaseElement::new(21),
        BaseElement::new(8),
        BaseElement::new(7),
        BaseElement::new(23),
        BaseElement::new(8),
        BaseElement::new(26),
        BaseElement::new(13),
    ],
    [
        BaseElement::new(13),
        BaseElement::new(10),
        BaseElement::new(9),
        BaseElement::new(7),
        BaseElement::new(6),
        BaseElement::new(22),
        BaseElement::new(21),
        BaseElement::new(8),
        BaseElement::new(7),
        BaseElement::new(23),
        BaseElement::new(8),
        BaseElement::new(26),
    ],
    [
        BaseElement::new(26),
        BaseElement::new(13),
        BaseElement::new(10),
        BaseElement::new(9),
        BaseElement::new(7),
        BaseElement::new(6),
        BaseElement::new(22),
        BaseElement::new(21),
        BaseElement::new(8),
        BaseElement::new(7),
        BaseElement::new(23),
        BaseElement::new(8),
    ],
    [
        BaseElement::new(8),
        BaseElement::new(26),
        BaseElement::new(13),
        BaseElement::new(10),
        BaseElement::new(9),
        BaseElement::new(7),
        BaseElement::new(6),
        BaseElement::new(22),
        BaseElement::new(21),
        BaseElement::new(8),
        BaseElement::new(7),
        BaseElement::new(23),
    ],
    [
        BaseElement::new(23),
        BaseElement::new(8),
        BaseElement::new(26),
        BaseElement::new(13),
        BaseElement::new(10),
        BaseElement::new(9),
        BaseElement::new(7),
        BaseElement::new(6),
        BaseElement::new(22),
        BaseElement::new(21),
        BaseElement::new(8),
        BaseElement::new(7),
    ],
];

// ROUND CONSTANTS
// ================================================================================================

/// Tip4' round constants;
const ARK: [[BaseElement; STATE_WIDTH]; NUM_ROUNDS] = [
    [
        BaseElement::new(6389859829374159539),
        BaseElement::new(17229380705756925792),
        BaseElement::new(1685504065768360085),
        BaseElement::new(11495601389602964927),
        BaseElement::new(14590499332882967787),
        BaseElement::new(1025546974836270924),
        BaseElement::new(6534388361641189904),
        BaseElement::new(83528329195322609),
        BaseElement::new(16825958830920599305),
        BaseElement::new(14367535913531541706),
        BaseElement::new(9373882113169352939),
        BaseElement::new(15206069388815085038),
    ],
    [
        BaseElement::new(11737921204895515252),
        BaseElement::new(16675938809814109201),
        BaseElement::new(3081949179470873665),
        BaseElement::new(14408516998600831324),
        BaseElement::new(1734666691041050162),
        BaseElement::new(6965669371710145373),
        BaseElement::new(2627545038496711157),
        BaseElement::new(4131567233355311863),
        BaseElement::new(6663606854236214944),
        BaseElement::new(4356348818072343221),
        BaseElement::new(16977499875487476107),
        BaseElement::new(7133528005210469376),
    ],
    [
        BaseElement::new(6576072100995654613),
        BaseElement::new(10365773285517510429),
        BaseElement::new(5099169911046663151),
        BaseElement::new(8204666535893592250),
        BaseElement::new(8197747920434370952),
        BaseElement::new(14055144323705103704),
        BaseElement::new(3823557958183613205),
        BaseElement::new(13605446815893856691),
        BaseElement::new(5735413688512934979),
        BaseElement::new(6851805580612830374),
        BaseElement::new(15826549711290838538),
        BaseElement::new(12685129780251896487),
    ],
    [
        BaseElement::new(14630567176945838998),
        BaseElement::new(10914381549746387617),
        BaseElement::new(7421170053404255007),
        BaseElement::new(729729048984908966),
        BaseElement::new(14034826357100489727),
        BaseElement::new(1838990003271062491),
        BaseElement::new(1898220471987185721),
        BaseElement::new(16545038689555322876),
        BaseElement::new(8094685071437843378),
        BaseElement::new(12850344507116018299),
        BaseElement::new(12312176261866659209),
        BaseElement::new(15951038048410299964),
    ],
    [
        BaseElement::new(10272996019760781151),
        BaseElement::new(15604072713712758006),
        BaseElement::new(16041565754647052886),
        BaseElement::new(14050558224171946805),
        BaseElement::new(7797167814871903221),
        BaseElement::new(14703987605672997265),
        BaseElement::new(489387044347918799),
        BaseElement::new(16525276081598735217),
        BaseElement::new(18290527903441405997),
        BaseElement::new(14091853142054427371),
        BaseElement::new(10013792448977695701),
        BaseElement::new(5570882487753758060),
    ],
];
