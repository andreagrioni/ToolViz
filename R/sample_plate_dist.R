#' Visualize 96-wells plate and sample positions.
#' 
#' @param table digested_adat table.
#' @param color_by color samples by column.
#' @param text_by annotate samples by column.
#' @param title ggplot title.
#' @param nrow TODO.
#' @param ncol TODO.
#' @param text_size text fond size.
#' @param debug verbose (bool).
#' @return ggplot class object
#' @examples
#' draw_heatmap(digested_adat$complete, "PlateId")
#' @export
#' @importFrom forcats as_factor
draw_heatmap <- function(
  table, 
  color_by,
  text_by=NULL,
  title="heatmap",
  nrow=NULL,
  ncol=NULL,
  text_size=4
)
{
  
  working_tbl <- table %>% dplyr::rename( color_col = .data[[color_by]])
  # # fix colors
  n_colors <- length(working_tbl %>% pull(color_col) %>% unique())
  colors <- setNames(object = scales::hue_pal()(n_colors), nm = working_tbl %>% 
                       pull(color_col) %>% 
                       unique())
  targets <- c(
    "PlateId", "row_location", "column_location", "PlatePosition", "color_col", text_by)
  
  working_tbl %<>% separate(
    PlatePosition,
    into=c("column_location", "row_location"),
    sep="(?<=[A-Za-z])(?=[0-9])", remove=FALSE) %>%
    mutate(row_location=as.numeric(row_location)) %>%
    arrange(desc(column_location), row_location) %>% 
    mutate(row_location=as.character(row_location))
  
  column_location = rev(LETTERS[seq( from = 1, to = 8 )])
  row_location = seq(c(1:12)) %>% as.character()
  dummy_grid <- expand.grid(
    column_location=column_location, row_location=row_location
  )
  long_df <- tibble()
  
  for (plate in working_tbl %>% pull(PlateId) %>% unique()) {
    original <- working_tbl %>% filter(PlateId == plate)
    tmp <- dummy_grid %>%
      left_join(original, by=c("column_location", "row_location")) %>%
      select( all_of(targets))
    long_df %<>% bind_rows(tmp)
  }
  
  multiplot_ <- long_df %>%
    filter(!is.na(PlateId)) %>%
    ggplot(
      aes(
        x=as_factor(row_location),
        y=column_location,
        fill=color_col
      )
    ) + 
    geom_tile() +
    #scale_fill_manual(values = colors) +
    ggtitle(title) + 
    labs(
      x="row position",
      y="column position"
    )
  
  multiplot_ <- tryCatch(
    {
      scales::train_discrete(table$AGE)
      multiplot_ + scale_fill_manual(values = colors)
    },
    error=function(cond) {
      multiplot_
    })
  
  if (!is.null(text_by)) {
    multiplot_ <- multiplot_ + 
      geom_text(aes(label=.data[[text_by]]), size=text_size)
  }
  
  multiplot_ <- multiplot_ + facet_wrap(~ PlateId)
  
  return(multiplot_)
}