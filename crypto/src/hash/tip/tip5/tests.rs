// Copyright (c) Facebook, Inc. and its affiliates.
//
// This source code is licensed under the MIT license found in the
// LICENSE file in the root directory of this source tree.

use super::{BaseElement, ElementHasher, FieldElement, Tip5_320, MDS, RATE_WIDTH, STATE_WIDTH};
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
        BaseElement::new(12),
        BaseElement::new(13),
        BaseElement::new(14),
        BaseElement::new(15),
    ];

    Tip5_320::apply_permutation(&mut state);

    let expected = vec![
        BaseElement::new(10175906631542820923),
        BaseElement::new(3855840548738021448),
        BaseElement::new(2025673606217227566),
        BaseElement::new(12844804840549141103),
        BaseElement::new(14675156848853604811),
        BaseElement::new(14826089615852996663),
        BaseElement::new(10347352866268213224),
        BaseElement::new(11152115716500330301),
        BaseElement::new(5111971323170855292),
        BaseElement::new(8786047983351640631),
        BaseElement::new(13278464272689444323),
        BaseElement::new(6663957460006323152),
        BaseElement::new(2554291996928436336),
        BaseElement::new(1966561718193927968),
        BaseElement::new(14725516211967935089),
        BaseElement::new(12331578240098214586),
    ];
    assert_eq!(expected, state);
}

#[test]
fn hash() {
    let state: [BaseElement; RATE_WIDTH] = [
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
    ];

    let result = Tip5_320::hash_elements(&state);

    let expected = vec![
        BaseElement::new(11606474556183478127),
        BaseElement::new(7414287774634619863),
        BaseElement::new(9612525545596271753),
        BaseElement::new(4029871457117740896),
        BaseElement::new(13014278221680054371),
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
        Tip5_320::apply_mds(&mut v2);

        prop_assert_eq!(v1, v2);
    }
}
