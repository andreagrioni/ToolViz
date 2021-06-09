ToolViz
================
Grioni, Andrea
2021-06-09

``` r
library(ToolViz)
#> Warning: replacing previous import 'dplyr::collapse' by 'glue::collapse' when
#> loading 'ToolViz'
data("test_data")
```

# Rain Cloud Plot

## Visualize time-dependent protein concentration

The function is adapted from Cédric Scherer’s blog
[here](https://www.cedricscherer.com/2021/06/06/visualizing-distributions-with-raincloud-plots-and-how-to-create-them-with-ggplot2/)

Function arguments:

-   data: tibble class object as sample X features. Features can be gene
    names, probes ids, or gene sets.
-   target: column name corresponding to the feature to be visualize
    (e.g. Probe ID, Gene Name),
-   annotation: if features are probes ids, you can provide a tibble
    class object with probes ids and corresponding EntrezGeneSymbol.
-   x: the column name of the time event (e.g. Visit, Days)..
-   group\_ids: how to separate samples in the dataset
    (e.g. Treatment\_code, ARM code)..
-   spaghetti: add lines to connect each observation over time.
-   colspaghetti: color lines according to a feature (e.g. SEX).

The annotation must have a column named `SeqId` with ProbeIds as values
and a column named `EntrezGeneSymbol` with the corresponding gene name.

``` r
plot_raincloud(
  test_data,
  #annotation=anno_v4,
  target='seq.10008.43',
  x='Sample...Visit.Description..clin.event.',
  group_ids = 'ARM',
  spaghetti="Sample...Screening.ID",
  colspaghetti='PlateId'
)
#> Coordinate system already present. Adding new coordinate system, which will replace the existing one.
```

![](/tmp/RtmpUCrPex/preview-c12c69c03a4a.dir/ToolViz_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
#ggsave("PATH/TO/FILE.png", dpi=300, width = 15, height=15)
```

# Distribution Plot

## Visualize Protein Concentration Changes over time as Violin Plot

Function arguments:

-   data: tibble class object as sample X features. Features can be gene
    names, probes ids, or gene sets.
-   target: column name corresponding to the feature to be visualize
    (e.g. Probe ID, Gene Name),
-   annotation: if features are probes ids, you can provide a tibble
    class object with probes ids and corresponding EntrezGeneSymbol.
-   x: the column name of the time event (e.g. Visit, Days)..
-   group\_ids: how to separate samples in the dataset
    (e.g. Treatment\_code, ARM code)..
-   spaghetti: add lines to connect each observation over time.
-   colspaghetti: color lines according to a feature (e.g. SEX).

The annotation must have a column named `SeqId` with ProbeIds as values
and a column named `EntrezGeneSymbol` with the corresponding gene name.

``` r
ToolViz::plot_dist(
  test_data,
  #annotation=anno_v4,
  target='seq.10008.43',
  x='Sample...Visit.Description..clin.event.',
  #rows_by="ARM",
  cols_by='ARM',
  kind="grid",
  x_rot=45,
  scales_grid='free_y',
  spaghetti="Sample...Screening.ID",
  colspaghetti='SEX',
  stats=TRUE,
  x_lab="time")
#> Warning: `fun.y` is deprecated. Use `fun` instead.
```

![](/tmp/RtmpUCrPex/preview-c12c69c03a4a.dir/ToolViz_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
#ggsave("PATH/TO/FILE.png", dpi=300, width = 15, height=15)
```

# draw\_heatmap

## Visually represent the plate position of each sample in the dataset

Function useful to debug experiment where some of the samples were
corrupted. The heatmap helps to identify if corrupted samples originated
from a cluster/region of the plate.

The function looks for the column named `PlatePosition` and separate it
into `column_location` and `row_location`.

Function arguments:

-   table = table obtained from somascan adat file

``` r
draw_heatmap(
  table=test_data,
  color_by="ARM",
  text_by="SEX",
  title="heatmap 96-well plate position",
  text_size=4)
```

![](/tmp/RtmpUCrPex/preview-c12c69c03a4a.dir/ToolViz_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
#ggsave("PATH/TO/FILE.png", dpi=300, width = 15, height=15)
```
