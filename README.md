# p2

[![Build Status](https://travis-ci.org/trink/p2.svg?branch=master)](https://travis-ci.org/trink/p2)
[![Coverage Status](https://coveralls.io/repos/github/trink/p2/badge.svg?branch=master)](https://coveralls.io/github/trink/p2?branch=master)

Piecewise Parabolic Prediction (P2)

The [p2](http://www.cs.wustl.edu/~jain/papers/ftp/psqr.pdf) algorithm for
dynamic calculation of quantiles and histograms without storing observation.

## Installation

```
cargo install --git https://github.com/trink/parquetfmt

```
## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
p2 = "0.1"
```

[Full Documentation](https://trink.github.io/p2).


## License

MPL. See [LICENSE](LICENSE).  
Copyright (c) 2018 Mike Trinkala <trink@acm.org>

