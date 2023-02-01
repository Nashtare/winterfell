// Copyright (c) Facebook, Inc. and its affiliates.
//
// This source &code is licensed under the MIT license found in the
// LICENSE file in the root directory of this source tree.

// FFT-BASED MDS MULTIPLICATION HELPER FUNCTIONS
// ================================================================================================

use math::{
    fields::f64::{self, BaseElement},
    FieldElement,
};

pub(crate) fn mds_multiply(state: &mut [BaseElement; 16]) {
    const SHIFTS: [u8; 16] = [4, 1, 4, 3, 3, 7, 0, 5, 1, 5, 0, 2, 6, 2, 4, 1];
    let mut array: [u128; 16] = [0; 16];
    ntt_noswap(state);

    for i in 0..16 {
        array[i] = (state[i].inner() as u128) << SHIFTS[i];
    }
    for i in 0..16 {
        state[i] = BaseElement::from_mont(f64::mont_red_cst(array[i]));
    }

    intt_noswap(state);
}

#[inline(always)]
fn ntt_noswap(x: &mut [BaseElement]) {
    const POWERS_OF_OMEGA_BITREVERSED: [BaseElement; 8] = [
        BaseElement::ONE,
        BaseElement::new(281474976710656),
        BaseElement::new(18446744069397807105),
        BaseElement::new(18446742969902956801),
        BaseElement::new(17293822564807737345),
        BaseElement::new(4096),
        BaseElement::new(4503599626321920),
        BaseElement::new(18446744000695107585),
    ];

    // outer loop iteration 1
    for j in 0..8 {
        let u = x[j];
        let v = x[j + 8] * BaseElement::ONE;
        x[j] = u + v;
        x[j + 8] = u - v;
    }

    // outer loop iteration 2
    for (i, zeta) in POWERS_OF_OMEGA_BITREVERSED.iter().enumerate().take(2) {
        let s = i * 8;
        for j in s..(s + 4) {
            let u = x[j];
            let v = x[j + 4] * *zeta;
            x[j] = u + v;
            x[j + 4] = u - v;
        }
    }

    // outer loop iteration 3
    for (i, zeta) in POWERS_OF_OMEGA_BITREVERSED.iter().enumerate().take(4) {
        let s = i * 4;
        for j in s..(s + 2) {
            let u = x[j];
            let v = x[j + 2] * *zeta;
            x[j] = u + v;
            x[j + 2] = u - v;
        }
    }

    // outer loop iteration 4
    for (i, zeta) in POWERS_OF_OMEGA_BITREVERSED.iter().enumerate().take(8) {
        let s = i * 2;
        let u = x[s];
        let v = x[s + 1] * *zeta;
        x[s] = u + v;
        x[s + 1] = u - v;
    }
}

#[inline(always)]
fn intt_noswap(x: &mut [BaseElement]) {
    const POWERS_OF_OMEGA_INVERSE: [BaseElement; 8] = [
        BaseElement::ONE,
        BaseElement::new(68719476736),
        BaseElement::new(1099511627520),
        BaseElement::new(18446744069414580225),
        BaseElement::new(18446462594437873665),
        BaseElement::new(18442240469788262401),
        BaseElement::new(16777216),
        BaseElement::new(1152921504606846976),
    ];

    // outer loop iteration 1
    {
        // while k < 16 as usize
        // inner loop iteration 1
        {
            let u = x[1];
            let v = x[0];
            x[1] = v - u;
            x[0] = v + u;
        }

        // inner loop iteration 2
        {
            let u = x[2 + 1];
            let v = x[2];
            x[2 + 1] = v - u;
            x[2] = v + u;
        }

        // inner loop iteration 3
        {
            let u = x[4 + 1];
            let v = x[4];
            x[4 + 1] = v - u;
            x[4] = v + u;
        }

        // inner loop iteration 4
        {
            let u = x[6 + 1];
            let v = x[6];
            x[6 + 1] = v - u;
            x[6] = v + u;
        }

        // inner loop iteration 5
        {
            let u = x[8 + 1];
            let v = x[8];
            x[8 + 1] = v - u;
            x[8] = v + u;
        }

        // inner loop iteration 6
        {
            let u = x[10 + 1];
            let v = x[10];
            x[10 + 1] = v - u;
            x[10] = v + u;
        }

        // inner loop iteration 7
        {
            let u = x[12 + 1];
            let v = x[12];
            x[12 + 1] = v - u;
            x[12] = v + u;
        }

        // inner loop iteration 7
        {
            let u = x[14 + 1];
            let v = x[14];
            x[14 + 1] = v - u;
            x[14] = v + u;
        }
    }

    // outer loop iteration 2
    {
        // while k < 16 as usize
        // inner loop iteration 1
        {
            for j in 0..2 {
                let zeta = POWERS_OF_OMEGA_INVERSE[4 * j];
                {
                    let u = x[j + 2] * zeta;
                    let v = x[j];
                    x[j + 2] = v - u;
                    x[j] = v + u;
                }
                // inner loop iteration 2
                {
                    let u = x[4 + j + 2] * zeta;
                    let v = x[4 + j];
                    x[4 + j + 2] = v - u;
                    x[4 + j] = v + u;
                }
                // inner loop iteration 3
                {
                    let u = x[8 + j + 2] * zeta;
                    let v = x[8 + j];
                    x[8 + j + 2] = v - u;
                    x[8 + j] = v + u;
                }
                // inner loop iteration 4
                {
                    let u = x[12 + j + 2] * zeta;
                    let v = x[12 + j];
                    x[12 + j + 2] = v - u;
                    x[12 + j] = v + u;
                }
            }
        }
    }

    // outer loop iteration 3
    {
        // while k < 16 as usize
        {
            for j in 0..4 {
                let zeta = POWERS_OF_OMEGA_INVERSE[2 * j];
                // inner loop iteration 1
                {
                    let u = x[j + 4] * zeta;
                    let v = x[j];
                    x[j + 4] = v - u;
                    x[j] = v + u;
                }
                // inner loop iteration 2
                {
                    let u = x[8 + j + 4] * zeta;
                    let v = x[8 + j];
                    x[8 + j + 4] = v - u;
                    x[8 + j] = v + u;
                }
            }
        }
    }

    // outer loop iteration 4
    {
        for j in 0..8 {
            let zeta = POWERS_OF_OMEGA_INVERSE[j];
            let u = x[j + 8] * zeta;
            let v = x[j];
            x[j + 8] = v - u;
            x[j] = v + u;
        }
    }
}
