// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! P2 algorithms for quantiles and histograms.

#![deny(
    missing_docs,
    missing_debug_implementations,
    missing_copy_implementations,
    trivial_casts,
    trivial_numeric_casts,
    unsafe_code,
    unstable_features,
    unused_import_braces,
    unused_qualifications
)]

extern crate simple_error;

pub use self::histogram::Histogram;
pub use self::quantile::Quantile;

mod histogram;
mod quantile;

fn parabolic(i: usize, d: f64, q: &[f64], n: &[f64]) -> f64 {
    return q[i]
        + d / (n[i + 1] - n[i - 1])
            * ((n[i] - n[i - 1] + d) * (q[i + 1] - q[i]) / (n[i + 1] - n[i])
                + (n[i + 1] - n[i] - d) * (q[i] - q[i - 1]) / (n[i] - n[i - 1]));
}

fn linear(i: usize, d: f64, q: &[f64], n: &[f64]) -> f64 {
    return q[i] + d * (q[(i as f64 + d) as usize] - q[i]) / (n[(i as f64 + d) as usize] - n[i]);
}

#[cfg(test)]
mod test_data {
    pub const OBS: [f64; 20] = [
        0.02, 0.15, 0.74, 3.39, 0.83, 22.37, 10.15, 15.43, 38.62, 15.92, 34.60, 10.28, 1.47, 0.40,
        0.05, 11.39, 0.27, 0.42, 0.09, 11.37,
    ];
    pub const MIN_RESULTS: [f64; 5] = [0.02, 0.15, 0.74, 0.83, 3.39];
    pub const FULL_RESULTS: [f64; 5] = [0.02, 0.493895, 4.44063, 17.2039, 38.62];
    pub const COUNT_RESULTS: [usize; 5] = [1, 6, 10, 16, 20];
}
