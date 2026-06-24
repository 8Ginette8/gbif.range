# Draw a Custom Legend Bar

Internal plotting helper adapted from Philipp Brun.

## Usage

``` r
cscl(
  colors,
  crds,
  horiz = FALSE,
  zrng = c(0, 100),
  at = 10 * 0:10,
  labs = NA,
  tickle = 0.2,
  title = 1,
  lablag = 1,
  titlag = 2,
  box = TRUE,
  breaks,
  cx = 0.8,
  tria = "n"
)
```

## Arguments

- colors:

  Vector of fill colors.

- crds:

  Numeric vector of length four giving the plotting coordinates of the
  legend box.

- horiz:

  Logical. Should the legend be drawn horizontally?

- zrng:

  Numeric vector of length two giving the value range represented by the
  legend.

- at:

  Numeric values at which tick marks should be drawn.

- labs:

  Optional labels for `at`. If `NA`, `at` is used.

- tickle:

  Numeric tick-mark length scaling factor.

- title:

  Legend title.

- lablag:

  Numeric offset multiplier for tick labels.

- titlag:

  Numeric offset multiplier for the title.

- box:

  Logical. Should a border be drawn around the legend?

- breaks:

  Optional custom break points. Must have length `length(colors) + 1` if
  supplied.

- cx:

  Numeric text expansion factor.

- tria:

  Character flag controlling triangular tips at the legend ends.
