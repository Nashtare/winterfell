// Copyright (c) Facebook, Inc. and its affiliates.
//
// This source code is licensed under the MIT license found in the
// LICENSE file in the root directory of this source tree.

use super::{BaseElement, ElementHasher, FieldElement, Tip4_256, MDS, STATE_WIDTH};
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

    Tip4_256::apply_permutation(&mut state);

    let expected = vec![
        BaseElement::new(7296980799497271299),
        BaseElement::new(6932375651683702735),
        BaseElement::new(3822213789028045974),
        BaseElement::new(505530293833729654),
        BaseElement::new(17605892691164638210),
        BaseElement::new(1470910630400469657),
        BaseElement::new(4639730692758713667),
        BaseElement::new(14190377574028754122),
        BaseElement::new(17996505783824879381),
        BaseElement::new(961936667785551774),
        BaseElement::new(12287826717832661200),
        BaseElement::new(11356005505906191738),
        BaseElement::new(10120826246370087028),
        BaseElement::new(1060843996108630575),
        BaseElement::new(6359159595152998556),
        BaseElement::new(15525397670994864633),
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
        BaseElement::new(12),
        BaseElement::new(13),
        BaseElement::new(14),
        BaseElement::new(15),
    ];

    let result = Tip4_256::hash_elements(&state);

    let expected = vec![
        BaseElement::new(2049960550193020876),
        BaseElement::new(17337156043624831090),
        BaseElement::new(16635828423870957893),
        BaseElement::new(12786227377510387913),
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
        Tip4_256::apply_mds(&mut v2);

        prop_assert_eq!(v1, v2);
    }
}
