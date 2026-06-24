# Objective Function for `make_blocks()`

Internal helper used by
[`make_blocks()`](https://8ginette8.github.io/gbif.range/reference/make_blocks.md)
to distribute clusters across folds as evenly as possible.

## Usage

``` r
optme(x, nms, grps, tot)
```

## Arguments

- x:

  Candidate allocation of the remaining clusters to folds.

- nms:

  Named vector of cluster sizes.

- grps:

  Current fold assignments built before optimization.

- tot:

  Total number of observations across all clusters.
