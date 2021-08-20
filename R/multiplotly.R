#' multiplot PCs interactive with plotly
#'
#' @param table matrix of aptamer values and metadata
#' @param color color by feature
#' @param pca_obj PCAtools::pca object
#' @param marker_size marker size
#' @param show_legend Def. TRUE
#' @return multiple plotly scatterplot of components
#' @examples
#' multiplotly(t_wdf$df, color="ARM", pca_obj = t_wdf$pca_obj)
#' @import plotly glue
#' @export
multiplotly <- function(table, color="ARM", pca_obj=NULL, marker_size=4, show_legend=TRUE) {

  pc1_lab = if (!is.null(pca_obj)) {glue("PC1\n{round(pca_obj$variance[['PC1']], 2)}%") } else {"PC1"}
  pc2_lab = if (!is.null(pca_obj)) {glue("PC2\n{round(pca_obj$variance[['PC2']], 2)}%") } else {"PC2"}
  pc3_lab = if (!is.null(pca_obj)) {glue("PC3\n{round(pca_obj$variance[['PC3']], 2)}%") } else {"PC3"}
  pc4_lab = if (!is.null(pca_obj)) {glue("PC4\n{round(pca_obj$variance[['PC4']], 2)}%") } else {"PC4"}
  pc5_lab = if (!is.null(pca_obj)) {glue("PC5\n{round(pca_obj$variance[['PC5']], 2)}%") } else {"PC5"}

  fig <- table %>%
    plot_ly() %>%
    add_trace(
      type = 'splom',
      dimensions = list(
        list(label=glue("{pc1_lab}"), values=~PC1),
        list(label=glue("{pc2_lab}"), values=~PC2),
        list(label=glue("{pc3_lab}"), values=~PC3),
        list(label=glue("{pc4_lab}"), values=~PC4),
        list(label=glue("{pc5_lab}"), values=~PC5)
      ),
      color=~.data[[color]],
      marker = list(
        size = marker_size
      ),
      data=table,
      hoverinfo = 'text',
      showlegend = show_legend,
      text=~paste(
        '</br> PlateId: ', PlateId,
        '</br> ScannerID: ', ScannerID,
        '</br> SampleId: ', SampleId,
        '</br> ARM: ', ARM,
        '</br> AGE: ', AGE,
        '</br> SEX: ', SEX
      )
    ) %>%
    style(
      diagonal = list(visible = F),
      showlowerhalf = F
    )
  return(fig)
}
