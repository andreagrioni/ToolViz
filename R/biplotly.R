#' PCA biplot interactive with plotly
#'
#' @param table matrix of aptamer values and metadata
#' @param x x value
#' @param y y value
#' @param color color by feature
#' @param pca_obj PCAtools::pca object
#' @return plotly scatterplot of components
#' @examples
#' biplotly(t_wdf$df, x="PC1", y="PC2", color="ARM", pca_obj = t_wdf$pca_obj)
#' @import plotly glue
#' @export
biplotly <- function(table, x="PC1", y="PC2", color="ARM", pca_obj=NULL, title=NULL) {
  x_lab = if (!is.null(pca_obj)) {glue("{x} % variance: {round(pca_obj$variance[[x]], 2)}") } else {glue("{x}")}
  y_lab = if (!is.null(pca_obj)) {glue("{y} % variance: {round(pca_obj$variance[[y]], 2)}") } else {glue("{y}")}
  title_ = if (is.null(title)) { glue("PCs Comp {x} vs {y}") } else { title }
  # draw plot
  plot_ly(
    data=table,
    x = ~.data[[x]],
    y = ~.data[[y]],
    type = 'scatter',
    mode = 'markers',
    color = ~.data[[color]]
  ) %>%
    plotly::add_trace(
      data=table,
      hoverinfo = 'text',
      showlegend = F,
      text=~paste(
        '</br> PlateId: ', PlateId,
        '</br> ScannerID: ', ScannerID,
        '</br> SampleId: ', SampleId,
        '</br> ARM: ', ARM,
        '</br> AGE: ', AGE,
        '</br> SEX: ', SEX
      )) %>%
    plotly::layout(
      title = glue(title_),
      xaxis = list(title = glue(x_lab)),
      yaxis = list(title = glue(y_lab))
    )
}
