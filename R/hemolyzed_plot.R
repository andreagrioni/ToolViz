#' plot hemolyzed samples vs normal.
#'
#' @param accepted table of accepted samples
#' @param rm_visual table of filtered by visual
#' @param rm_hba table of filtered by aptamer value
#' @param hemo_id aptamer colname
#' @param interactive draw an interactive plotly
#' @param title set main title
#' @return boxplot of aptamer distribution by class
#' @importFrom ggplot2 ggplot geom_boxplot
#' @import dplyr
#' @import plotly
#' @export
hemo_plot <- function (
  accepted,
  rm_visual,
  rm_hba,
  hemo_id="seq.4915.64",
  interactive=TRUE,
  title
  )
  {
    target_cols <- c("SampleId", hemo_id)
    accepted_ <- accepted %>%
      select(all_of(target_cols)) %>%
      mutate(flag = "accepted") %>%
      mutate(status="empty note")
    rm_visual_ <- rm_visual %>%
      select(all_of(target_cols)) %>%
      mutate(flag = "removed") %>%
      mutate(status="visual check")
    rm_hba_ <- rm_hba %>%
      select(all_of(target_cols)) %>%
      mutate(flag = "removed") %>%
      mutate(status="hemo check")

    tmp_data <- rbind(accepted_, rm_visual_, rm_hba_)

    if (interactive) {

      fig <- plot_ly(
        tmp_data,
        x = ~flag,
        y = ~log2(.data[[hemo_id]]),
        type = "box",
        marker = list( size = 10),
        boxpoints = "all",
        jitter = 0.4,
        pointpos = 0,
        hoverinfo = "text",
        text = ~paste('</br> SampleId: ', SampleId,
                      '</br> Flag: ', flag,
                      '</br> Status:', status,
                      '</br> RFU:', .data[[hemo_id]])) %>%
        layout(title=title)
    } else {
      fig <- tmp_data %>%
        ggplot(aes(x = flag, y = log2(.data[[hemo_id]]), hue = flag)
        ) + geom_boxplot() + ggtitle(title)
    }
    return(fig)
  }
