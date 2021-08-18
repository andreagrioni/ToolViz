#' plot hemolyzed samples vs normal.
#'
#' @param accepted table of accepted samples
#' @param rm_visual table of filtered by visual
#' @param rm_hba table of filtered by aptamer value
#' @param hemo_id aptamer colname
#' @return boxplot of aptamer distribution by class
#' @examples
#' @importFrom ggplot2 ggplot geom_boxplot
#' @import dplyr
#' @export
hemo_plot <- function(
  accepted,
  rm_visual,
  rm_hba,
  hemo_id="seq.4915.64"
  ) {

  target_cols <- c("SampleId", hemo_id)
  accepted_ <- accepted %>%
    select(all_of(target_cols)) %>%
    mutate(flag="accepted")
  rm_visual_ <- rm_visual %>%
    select(all_of(target_cols)) %>%
    mutate(flag="visual")
  rm_hba_ <- rm_hba %>%
    select(all_of(target_cols)) %>%
    mutate(flag="hemo")

  tmp_data <- rbind(accepted_, rm_visual_, rm_hba_)

  plot <- tmp_data %>%
    ggplot(aes(x=flag, y=.data[[hemo_id]], hue=flag)) +
    geom_boxplot()
  return(plot)
}
