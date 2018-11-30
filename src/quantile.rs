// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! The p2 algorithm for quantiles using five markers.

use std::cmp::Ordering;
use std::f64;
use std::fmt;

const QUANTILE_MARKERS: usize = 5;

/// P2 Quantile Data Structure
#[derive(Clone, Copy)]
pub struct Quantile {
    q: [f64; QUANTILE_MARKERS],
    n: [f64; QUANTILE_MARKERS], // this is an integer but to avoid a lot of casting it is made a float
    n1: [f64; QUANTILE_MARKERS],
    p: f32,
    cnt: u8,
}

impl fmt::Debug for Quantile {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({} p_quantile)", self.p)
    }
}

impl Quantile {
    /// Resets the quantile back to its initial state.
    pub fn clear(&mut self) -> &mut Self {
        self.q[0] = 0.0;
        self.q[1] = 0.0;
        self.q[2] = 0.0;
        self.q[3] = 0.0;
        self.q[4] = 0.0;

        self.n[0] = 1.0;
        self.n[1] = 2.0;
        self.n[2] = 3.0;
        self.n[3] = 4.0;
        self.n[4] = 5.0;

        self.n1[0] = 1.0;
        self.n1[1] = 1.0 + 2.0 * self.p as f64;
        self.n1[2] = 1.0 + 4.0 * self.p as f64;
        self.n1[3] = 3.0 + 2.0 * self.p as f64;
        self.n1[4] = 5.0;

        self.cnt = QUANTILE_MARKERS as u8;

        self
    }

    /// Constructor taking the p_quantile to calculate (e.g. 0.5 == median)
    pub fn new(p: f32) -> simple_error::SimpleResult<Quantile> {
        if p <= 0.0 || p >= 1.0 {
            return Err(simple_error::SimpleError::new(
                "p_quantile out of range 0 < p < 1",
            ));
        }
        let mut q = Quantile {
            q: [0.0; QUANTILE_MARKERS],
            n: [0.0; QUANTILE_MARKERS],
            n1: [0.0; QUANTILE_MARKERS],
            p: p,
            cnt: 0,
        };
        q.clear();
        Ok(q)
    }

    /// Adds a value to a histogram, NAN is ignored.
    pub fn add(&mut self, x: f64) -> f64 {
        if x.is_nan() {
            if self.cnt == 0 {
                return self.q[2];
            }
            return x;
        }
        if self.cnt > 0 {
            self.cnt -= 1;
            self.q[self.cnt as usize] = x;
            if self.cnt == 0 {
                self.q
                    .sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
                return self.q[2];
            }
            return f64::NAN;
        }

        let mut k = 0;
        if x < self.q[0] {
            self.q[0] = x;
            k = 1;
        } else if self.q[0] <= x && x < self.q[1] {
            k = 1;
        } else if self.q[1] <= x && x < self.q[2] {
            k = 2;
        } else if self.q[2] <= x && x < self.q[3] {
            k = 3;
        } else if self.q[3] <= x && x <= self.q[4] {
            k = 4;
        } else if self.q[4] < x {
            self.q[4] = x;
            k = 4;
        }

        for i in k..QUANTILE_MARKERS {
            self.n[i] += 1.0;
        }

        self.n1[1] += self.p as f64 / 2.0;
        self.n1[2] += self.p as f64;
        self.n1[3] += (1.0 + self.p as f64) / 2.0;
        self.n1[4] += 1.0;

        for i in 1..QUANTILE_MARKERS - 1 {
            let mut d = self.n1[i] - self.n[i];
            if (d >= 1.0 && self.n[i + 1] - self.n[i] > 1.0)
                || (d <= -1.0 && self.n[i - 1] - self.n[i] < -1.0)
            {
                if d > 0.0 {
                    d = 1.0;
                } else {
                    d = -1.0;
                }
                let q1 = ::parabolic(i, d, &self.q, &self.n);
                if self.q[i - 1] < q1 && q1 < self.q[i + 1] {
                    self.q[i] = q1;
                } else {
                    self.q[i] = ::linear(i, d, &self.q, &self.n);
                }
                self.n[i] += d;
            }
        }
        return self.q[2];
    }

    /// Returns the estimated quantile value for the specified marker.
    /// * 1 = min
    /// * 2 = p/2
    /// * 3 = p
    /// * 4 = (1+p)/2
    /// * 5 = max
    pub fn estimate(&self, marker: usize) -> f64 {
        if self.cnt != 0 {
            return f64::NAN;
        }
        return self.q[marker - 1];
    }

    /// Returns the number of observations that are less than or equal to the marker (see estimate).
    pub fn count(&self, marker: usize) -> usize {
        if self.cnt != 0 {
            return 0;
        }
        return self.n[marker - 1] as usize;
    }
}

#[cfg(test)]
mod tests {
    use super::Quantile;
    use std::f64;
    use test_data as td;

    #[test]
    fn test_quantile_min() {
        let mut q = Quantile::new(0.5).unwrap();
        assert!(q.add(f64::NAN).is_nan());
        let mut i = 1;
        while i < super::QUANTILE_MARKERS {
            assert!(q.estimate(i).is_nan());
            assert_eq!(q.count(i), 0);
            q.add(td::OBS[i - 1]);
            i += 1;
        }
        q.add(td::OBS[i - 1]);
        let result = std::panic::catch_unwind(|| q.estimate(6));
        assert!(result.is_err());
        let result = std::panic::catch_unwind(|| q.count(6));
        assert!(result.is_err());

        assert_eq!(q.add(f64::NAN), td::MIN_RESULTS[2]);

        i = 0;
        for x in &td::MIN_RESULTS {
            let rpq = q.estimate(i + 1);
            assert!(
                (rpq - x).abs() < 0.00001,
                format!("marker: {} received:{} expected:{}", i + 1, rpq, x)
            );
            i += 1;
            assert_eq!(i, q.count(i));
        }
        q.clear();
        assert_eq!(0, q.count(5));
    }

    #[test]
    fn test_quantile() {
        let mut q = Quantile::new(0.5).unwrap();
        let mut median = 0.0;
        for x in &td::OBS {
            median = q.add(*x);
        }

        let mut i = 1;
        for x in &td::FULL_RESULTS {
            let rpq = q.estimate(i);
            assert!(
                (rpq - x).abs() < 0.00001,
                format!("marker: {} received:{} expected:{}", i, rpq, x)
            );
            if i == 3 {
                assert!(
                    (rpq - median).abs() < 0.00001,
                    format!("add() received:{} expected:{}", rpq, median)
                );
            }
            i += 1;
        }

        i = 1;
        for x in &td::COUNT_RESULTS {
            let cnt = q.count(i);
            assert_eq!(cnt, *x);
            i += 1;
        }
    }
}
