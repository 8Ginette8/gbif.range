# Draw a Figure Label at a Standard Position

Draw a Figure Label at a Standard Position

## Usage

``` r
fig_label(
  text,
  region = "figure",
  pos = "topleft",
  cex = NULL,
  margin = 0.02,
  ...
)
```

## Arguments

- text:

  Label text.

- region:

  Character string specifying whether the label should be placed
  relative to the `"figure"`, `"plot"`, or `"device"`.

- pos:

  Character string giving the label position.

- cex:

  Optional text expansion factor.

- margin:

  Numeric margin used to offset the label from the selected boundary.

- ...:

  Additional arguments passed to
  [`graphics::text()`](https://rdrr.io/r/graphics/text.html).

## Value

Invisibly, a numeric vector `c(x, y)` giving the x and y boundary
coordinates of the selected `region` used to position the label. Called
primarily for its side effect of drawing the label on the current plot
device.
