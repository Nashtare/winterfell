// Copyright (c) Facebook, Inc. and its affiliates.
//
// This source code is licensed under the MIT license found in the
// LICENSE file in the root directory of this source tree.

use super::{BaseElement, ElementHasher, FieldElement, Tip4p_256, MDS, STATE_WIDTH};
use proptest::prelude::*;

#[test]
fn apply_permutation() {
    let mut state: [BaseElement; STATE_WIDTH] = [
        BaseElement::new(0),
        BaseElement::new(1),
        BaseElement::new(2),
        BaseElement::new(3),
        BaseElement::new(4),
        BaseElement::new(5),
        BaseElement::new(6),
        BaseElement::new(7),
        BaseElement::new(8),
        BaseElement::new(9),
        BaseElement::new(10),
        BaseElement::new(11),
    ];

    Tip4p_256::apply_permutation(&mut state);

    let expected = vec![
        BaseElement::new(8086224146445274039),
        BaseElement::new(12620228105612859910),
        BaseElement::new(4429745645163147655),
        BaseElement::new(12827206290147492018),
        BaseElement::new(7103575185686863209),
        BaseElement::new(5938996934280238338),
        BaseElement::new(7458235737397060283),
        BaseElement::new(127950926479970750),
        BaseElement::new(433935175963827303),
        BaseElement::new(11405496933068372192),
        BaseElement::new(4026696970861104429),
        BaseElement::new(6779880475047698803),
    ];

    assert_eq!(expected, state);
}

#[test]
fn hash() {
    let state: [BaseElement; STATE_WIDTH] = [
        BaseElement::new(0),
        BaseElement::new(1),
        BaseElement::new(2),
        BaseElement::new(3),
        BaseElement::new(4),
        BaseElement::new(5),
        BaseElement::new(6),
        BaseElement::new(7),
        BaseElement::new(8),
        BaseElement::new(9),
        BaseElement::new(10),
        BaseElement::new(11),
    ];

    let result = Tip4p_256::hash_elements(&state);

    let expected = vec![
        BaseElement::new(10976450459897344865),
        BaseElement::new(91268347875488615),
        BaseElement::new(11593875934785294243),
        BaseElement::new(444105276802301255),
    ];

    assert_eq!(expected, result.as_elements());
}

#[inline(always)]
fn apply_mds_naive(state: &mut [BaseElement; STATE_WIDTH]) {
    let mut result = [BaseElement::ZERO; STATE_WIDTH];
    result.iter_mut().zip(MDS).for_each(|(r, mds_row)| {
        state.iter().zip(mds_row).for_each(|(&s, m)| {
            *r += m * s;
        });
    });
    *state = result;
}

proptest! {
    #[test]
    fn mds_freq_proptest(a in any::<[u64;STATE_WIDTH]>()) {

        let mut v1 = [BaseElement::ZERO;STATE_WIDTH];
        let mut v2;

        for i in 0..STATE_WIDTH {
            v1[i] = BaseElement::new(a[i]);
        }
        v2 = v1.clone();

        apply_mds_naive(&mut v1);
        Tip4p_256::apply_mds(&mut v2);

        prop_assert_eq!(v1, v2);
    }
}
