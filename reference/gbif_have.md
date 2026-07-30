# Is a Network-Dependent Object Usable?

Small predicate used to guard help examples that depend on GBIF or on a
downloaded ecoregion layer. Following CRAN policy,
[`get_gbif()`](https://8ginette8.github.io/gbif.range/reference/get_gbif.md)
and
[`read_ecoreg()`](https://8ginette8.github.io/gbif.range/reference/read_ecoreg.md)
fail gracefully rather than raising an error when the remote source is
unavailable: they return `NULL` or an empty table. The examples
therefore have to check before using the result, or an upstream outage
would turn into a check error.

## Usage

``` r
gbif_have(...)
```

## Arguments

- ...:

  One or more objects to test. Typically the value returned by
  [`get_gbif()`](https://8ginette8.github.io/gbif.range/reference/get_gbif.md)
  or
  [`read_ecoreg()`](https://8ginette8.github.io/gbif.range/reference/read_ecoreg.md).

## Value

`TRUE` if every object is non-`NULL` and, where it has rows, has at
least one; `FALSE` otherwise.

## Details

Accepts several objects at once and returns `TRUE` only if all of them
are present and non-empty.
