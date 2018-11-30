// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

//! The p2 algorithm for histograms using equiprobable cells.

use std::cmp::Ordering;
use std::f64;
use std::fmt;
use std::u16;
use std::vec::Vec;

/// P2 Histogram Data Structure
#[derive(Clone)]
pub struct Histogram {
    q: Vec<f64>,
    n: Vec<f64>, // this is an integer but to avoid a lot of casting it is made a float
    b: u16,
    cnt: u16,
}

impl fmt::Debug for Histogram {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "({} buckets)", self.b)
    }
}

impl Histogram {
    /// Resets the histogram back to its initial state.
    pub fn clear(&mut self) -> &mut Self {
        self.cnt = self.b + 1;
        for i in 0..self.cnt {
            self.q[i as usize] = 0.0;
        }

        for i in 0..self.cnt {
            self.n[i as usize] = (i + 1) as f64;
        }

        self
    }

    /// Constructor taking the number of buckets to allocate in the histogram ( 3 < buckets < u16::MAX - 1).
    pub fn new(buckets: u16) -> simple_error::SimpleResult<Histogram> {
        if buckets < 4 || buckets == u16::MAX {
            return Err(simple_error::SimpleError::new(
                "buckets out of range 3 < buckets < u16 MAX - 1",
            ));
        }
        let mut h = Histogram {
            q: vec![0.0; (buckets + 1) as usize],
            n: vec![0.0; (buckets + 1) as usize],
            b: buckets,
            cnt: 0,
        };
        h.clear();
        Ok(h)
    }

    /// Adds a value to a histogram, NAN is ignored.
    pub fn add(&mut self, x: f64) {
        if x.is_nan() {
            return;
        }
        if self.cnt > 0 {
            self.cnt -= 1;
            self.q[self.cnt as usize] = x;
            if self.cnt == 0 {
                self.q
                    .sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
            }
            return;
        }

        let mut k = 0;
        if x < self.q[0] {
            self.q[0] = x;
            k = 1;
        } else {
            for i in 0..self.b - 1 {
                if self.q[i as usize] <= x && x < self.q[(i + 1) as usize] {
                    k = i + 1;
                    break;
                }
            }
        }

        if k == 0 {
            if self.q[(self.b - 1) as usize] <= x && x <= self.q[self.b as usize] {
                k = self.b;
            } else if self.q[self.b as usize] < x {
                self.q[self.b as usize] = x;
                k = self.b;
            }
        }

        for i in k..self.b + 1 {
            self.n[i as usize] += 1.0;
        }

        for i in 1..self.b as usize {
            let n1 = 1.0 + i as f64 * (self.n[self.b as usize] - 1.0) / self.b as f64;
            let mut d = n1 - self.n[i];
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
    }

    /// Returns the estimated quantile value for the specified marker.
    /// * 1 = min
    /// * bucket + 1 = max
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
    use super::Histogram;
    use test_data as td;

    #[test]
    fn test_histogram_min() {
        let mut q = Histogram::new(4).unwrap();
        let mut i = 1;
        while i < 5 {
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
    fn test_histogram() {
        let mut q = Histogram::new(4).unwrap();
        for x in &td::OBS {
            q.add(*x);
        }

        let mut i = 1;
        for x in &td::FULL_RESULTS {
            let rpq = q.estimate(i);
            assert!(
                (rpq - x).abs() < 0.00001,
                format!("marker: {} received:{} expected:{}", i, rpq, x)
            );
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
