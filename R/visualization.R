#' reshape somascan table from wide to long
#'
#' @param data tibble class object as sample X features.
#' @param x x-axis value (def. aptamers).
#' @param y y-axis value (def. rfu).
#' @param target EntrezGeneSymbol or Aptamer SeqId.
#' @param annotation Somascan Annotation Table.
#' @param group_ids fill color of basic ggplot (aes) (def. NULL).
#' @param cols_by facet_by col by feature (def. groups).
#' @param rows_by facet_by row by feature (def. aptamers).
#' @param spaghetti group samples by feature (allow spaghetti plot).
#' @param colspaghetti color spaghetti line by feature (def. NULL).
#' @return object tidy table long format
#' @examples
#' wide_to_long_soma(data=table,
#' annotation=anno, target='CXCL10',
#' x='time_idx', rows_by="aptamers",
#' cols_by='groups', spaghetti='SUBJID')
#' @export
#' @importFrom tidyr pivot_longer
#' @importFrom purrr is_empty
#' @import glue dplyr
wide_to_long_soma <- function(
  data,
  x=NULL,
  y=NULL,
  target=NULL,
  annotation=NULL,
  group_ids=NULL,
  rows_by=NULL,
  cols_by=NULL,
  sample_id='SampleId',
  spaghetti=NULL,
  colspaghetti=NULL
  ) {

  meta = c(
    x, y,
    group_ids,
    rows_by,
    cols_by,
    sample_id,
    spaghetti,
    colspaghetti)

  if (!is.null(target)) {
    if (target %in% colnames(data) ) {
      target_ids <- target
      target_names <- NULL
    } else {
      target_names <- target
      target_ids <- NULL
    }
  } else {
    print("please provide aptamer SeqId or EntrezGeneSymbol")
  }
  # find SeqId from EntrezGeneSymbol
  if (!is_empty(target_names)) {
    if (!is_empty(annotation)) {
      target_ids <- annotation %>%
        filter(str_detect(
          EntrezGeneSymbol,
          paste(target_names, collapse = '|'))) %>%
        pull(SeqId)
    } else {
      print(
        glue("target {target_names} not in data columns. Please provide annotation table to match {target_names} with aptamer SeqId"))
      return()
    }
  }
  # retrieve targets only
  subtable <- data %>%
    select(
      all_of(target_ids), all_of(meta)) %>%
    pivot_longer(
      cols=all_of(target_ids),
      names_to='aptamers',
      values_to='rfu')
  # add annotation
  if (!is_empty(annotation)) {
    subtable <- subtable %>%
      inner_join(annotation,
                 by=c("aptamers" = "SeqId"))
    subtable <- subtable %>%
      unite(
        aptamers, aptamers, EntrezGeneSymbol)
  }
  return(subtable)
}
#' generate distribution plot aptamers
#'
#' @param data tibble class object as samplexfeatures.
#' @param x x-axis value (def. aptamers).
#' @param y y-axis value (def. rfu).
#' @param target EntrezGeneSymbol or Aptamer SeqId.
#' @param annotation Somascan Annotation Table.
#' @param group_ids TODO.
#' @param cols_by facet_by col by feature (def. groups).
#' @param rows_by facet_by row by feature (def. aptamers).
#' @param kind if split plot with facet_by (def. grid)
#' @param title main plot title.
#' @param x_lab x axis label.
#' @param y_lab y axis label.
#' @param x_rot rotation x axis label (def. 90).
#' @param debug verbose level (bool).
#' @param scales_grid Are scales shared across all facets (the default, "fixed"), or do they vary across rows ("free_x"), columns ("free_y"), or both rows and columns ("free")?
#' @param spaghetti group samples by feature (allow spaghetti plot).
#' @param colspaghetti color spaghetti line by feature (def. NULL).
#' @param stats add stat_summary (median) (bool, def. FALSE).
#' @return gglot object
#' @examples
#' plot_dist(data=table,
#' annotation=anno, target='CXCL10',
#' x='time_idx', rows_by="aptamers",
#' cols_by='groups', kind="grid", x_rot=90,
#' scales_grid = 'free_y', spaghetti='SUBJID',
#' stats=TRUE)
#' @export
#' @import glue ggplot2 dplyr
plot_dist <- function(
  data,
  x=NULL,
  y=NULL,
  target=NULL,
  annotation=NULL,
  group_ids=NULL,
  rows_by=NULL,
  cols_by=NULL,
  sample_id=NULL,
  kind='grid',
  title='',
  x_lab=x,
  y_lab=y,
  x_rot=90,
  debug=FALSE,
  scales_grid='fixed',
  spaghetti=NULL,
  colspaghetti=NULL,
  stats=FALSE
  ) {

  subtable <- wide_to_long_soma(
    data=data,
    x=x,
    y=y,
    target=target,
    annotation=annotation,
    group_ids=group_ids,
    rows_by=rows_by,
    cols_by=cols_by,
    sample_id=sample_id,
    spaghetti=spaghetti,
    colspaghetti=colspaghetti
  )

  if (is.null(y)) {
    y = 'rfu'
  }

  if (is.null(rows_by)) {
    rows_by = 'aptamers'
  }

  gg_ <- subtable %>%
    ggplot(aes(x=.data[[x]], y=.data[[y]]))

  if (!is_empty(group_ids)) {
    gg_ <- gg_ + geom_violin(
      aes(color=.data[[group_ids]]))
  } else {
    gg_ <- gg_ + geom_violin()
  }
  # add group lines to make spaghetti
  if (!is.null(spaghetti)) {
    if (!is.null(colspaghetti)) {
      gg_ <- gg_ + geom_line(
        aes(
          group=.data[[spaghetti]],
          color=.data[[colspaghetti]]
        ),
        size=0.5,
        alpha=0.5)
    } else {
      gg_ <- gg_ + geom_line(
        aes(group=.data[[spaghetti]]),
        size=0.5,
        alpha=0.5,
        color='grey')
    }
  }
  # make grid with facet_grid
  if (kind=='grid') {
    gg_ <- gg_ +
      facet_grid(
        scales=scales_grid,
        cols=vars(.data[[cols_by]]),
        rows=vars(.data[[rows_by]]))
  }
  if (stats) {
    gg_ <- gg_ + stat_summary(
      fun.y=median, geom="point", size=1, color="red")
  }
  # add labels
  #see https://ggplot2.tidyverse.org/reference/labs.html
  gg_ = gg_ + labs(
    title=title,
    x=x_lab,
    y=y_lab)
  # add theme
  gg_ = gg_ + theme(
    axis.text.x = element_text(angle = x_rot)) + theme_minimal()
  return(gg_)
  }
#' generate raincloud plot aptamers
#' inspired by https://www.cedricscherer.com/2021/06/06/visualizing-distributions-with-raincloud-plots-and-how-to-create-them-with-ggplot2/
#'
#' @param data tibble class object as samplexfeatures.
#' @param x x-axis value (def. aptamers).
#' @param y y-axis value (def. rfu).
#' @param target EntrezGeneSymbol or Aptamer SeqId.
#' @param annotation Somascan Annotation Table.
#' @param group_ids TODO.
#' @param title main plot title.
#' @param x_lab x axis label.
#' @param y_lab y axis label.
#' @param x_rot rotation x axis label (def. 90).
#' @param spaghetti group samples by feature (allow spaghetti plot).
#' @param colspaghetti color spaghetti line by feature (def. NULL).
#' @param stats add stat_summary (median) (bool, def. FALSE).
#' @return gglot object
#' @examples
#' raincloud(data=table,
#' annotation=anno, target='CXCL10',
#' x='time_idx',  group_ids="ARM", x_rot=90,
#' spaghetti='SUBJID', stats=TRUE)
#' @export
#' @import glue ggplot2 dplyr
plot_raincloud <- function(
  data,
  x=NULL,
  y=NULL,
  target=NULL,
  annotation=NULL,
  group_ids=NULL,
  title='',
  x_lab=x,
  y_lab=y,
  x_rot=90,
  spaghetti=NULL,
  colspaghetti=NULL,
  stats=FALSE
  ) {
  # retrive data from input table
  subtable <- wide_to_long_soma(
    data=data,
    x=x,
    y=y,
    target=target,
    annotation=annotation,
    group_ids=group_ids,
    spaghetti=spaghetti,
    colspaghetti=colspaghetti
  )

  if (is.null(y)) {
    y = 'rfu'
  }

  # basic plot
  gg_ <- subtable %>%
    ggplot(aes(x=.data[[x]], y=.data[[y]], fill=.data[[group_ids]]))
  # let it rain
  gg_ <- gg_ + ggdist::stat_halfeye(
    adjust = .5,
    width = .6,
    .width = 0,
    justification = -.3,
    point_colour = NA) +
    geom_boxplot(
      width = .25,
      outlier.shape = NA
    ) +
    geom_point(
      size = 1.3,
      alpha = .3,
      position = position_jitter(
        seed = 1, width = .1
      )) +
        coord_cartesian(clip = "off") +
        coord_flip() +
        theme_minimal()
  # add group lines to make spaghetti
  if (!is.null(spaghetti)) {
    if (!is.null(colspaghetti)) {
      gg_ <- gg_ + geom_line(
        aes(
          group=.data[[spaghetti]],
          color=.data[[colspaghetti]]
        ),
        size=0.2,
        alpha=0.8)
    } else {
      gg_ <- gg_ + geom_line(
        aes(group=.data[[spaghetti]]),
        size=0.2,
        alpha=0.8,
        color='grey')
    }
  }
  # add labels
  #see https://ggplot2.tidyverse.org/reference/labs.html
  gg_ = gg_ + labs(
    title=title,
    x=x_lab,
    y=y_lab)
  # add theme
  gg_ = gg_ + theme(
    axis.text.x=element_text(angle = x_rot))
  return(gg_)
}
#' Visualize 96-wells plate and sample positions.
#'
#' @param table digested_adat table.
#' @param color_by color samples by column.
#' @param text_by annotate samples by column.
#' @param title ggplot title.
#' @param text_size text fond size.
#' @return ggplot class object
#' @examples
#' draw_heatmap(digested_adat$complete, "PlateId")
#' @export
#' @importFrom forcats as_factor
#' @import tidyr dplyr
draw_heatmap <- function(
  table,
  color_by,
  text_by=NULL,
  title="heatmap",
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
  
  working_tbl %>% separate(
    PlatePosition,
    into=c("column_location", "row_location"),
    sep="(?<=[A-Za-z])(?=[0-9])", remove=FALSE) %>%
    mutate(row_location=as.numeric(row_location)) %>%
    arrange(desc(column_location), row_location) %>%
    mutate(row_location=as.character(row_location)) -> working_tbl
  
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
    long_df %>% bind_rows(tmp) -> long_df
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
#' create distribution plot for aptamers
#' 
#' @param table a table with expression values.
#' @param title main plot title.
#' @param colors_fix named vector with color groups.
#' @param color_by color by column value.
#' @return distribution plot
#' @examples
#' plot_distribution(digested_adat$complete, "a dist plot", color_by='ARMCD')
#' @export
#' @importFrom dplyr select starts_with
#' @importFrom tidyr pivot_longer 
#' @importFrom glue glue
#' @import ggplot2
plot_distribution <- function (
  table, title=NULL, colors_fix=NULL, color_by=NULL) {
  table %>% 
    select(starts_with('seq.'), `color_by` ) %>%
    pivot_longer(!`color_by`) %>% 
    ggplot() + 
    geom_density(aes_string(x="value", color=color_by )) + 
    ggtitle(glue({title})) + 
    xlab("intensity") + ylab("density") + 
    scale_color_manual(values=colors_fix)
}
#' compute correlation between aptamers in dataframe A and dataframe B
#' 
#' @param df_a dataframe A
#' @param df_b dataframe B.
#' @param corr_type correlation method (see ?correlate).
#' @param comp_name column name for correlation values.
#' @param sample_aptamers autocorrelate for N aptamers (1000)
#' @return correlation dataframe
#' @examples
#' autocorrelation(table_1, table_2)
#' @export
#' @importFrom corrr correlate
#' @importFrom dplyr pull
#' @importFrom tibble as_tibble add_column
autocorrelation <- function(
  df_a, df_b, corr_type="spearman", comp_name='correlation', sample_aptamers=1000) {

  get_corr <- function(aptamer, df_a, df_b, corr_type) {
    correlate(
      df_a %>% pull(aptamer),
      df_b %>% pull(aptamer), 
      method=corr_type, 
      quiet=TRUE
    ) %>% 
      pull(x)
  }
  aptamers_name = df_a %>% 
    select(starts_with('seq')) %>% 
    colnames() %>% 
    sample(size=sample_aptamers)
  correlation <- sapply(
    aptamers_name, get_corr, df_a=df_a, df_b=df_b, corr_type=corr_type )
  corr_data <- as_tibble(correlation) %>% 
   add_column('SeqId'=names(correlation)) %>%
   add_column('aptamers'=comp_name)
  return(corr_data)
}
#' Scatterplot of two input tables.
#' 
#' @param x first input table with aptamer seq columns
#' @param y second input table with aptamer seq columns
#' @param sample_size fit line on n sample (default=5)
#' @param title main plot title
#' @param aptamers_name vector of aptamers col names (extract from table if NULL)
#' @return ggplot object.
#' @examples
#' plot_distscatter(digested_adat$complete_log2, digested_adat$postproces)
#' @export
#' @import dplyr
#' @import ggplot2
#' @importFrom stats lm
plot_distscatter <- function(
  x=table_a, y=table_b, sample_size=5, title="scatterplot", aptamers_name=NULL) {
  if (is.null(aptamers_name)) {
    x %>% 
      select(starts_with('seq.')) %>% 
      colnames -> aptamers_name
  }
  x %>% 
    select(SampleId, all_of(aptamers_name)) %>% 
    slice_head(n=sample_size ) %>%
    pivot_longer(
      cols=!SampleId, 
      names_to = 'aptamers', 
      values_to = 'value_x'
    ) -> raw
  y %>% 
    select(SampleId, all_of(aptamers_name)) %>% 
    slice_head(n=sample_size) %>%
    pivot_longer(
      cols=!SampleId, 
      names_to='aptamers', 
      values_to='value_y'
    ) -> postprocess
  raw %>% 
    bind_cols(
      postprocess %>% 
        select(value_y)
    ) -> methodData
  fit <- stats::lm(
    methodData$value_x ~ methodData$value_y
  )
  r2_value = format(
    summary(fit)$adj.r.squared, digits=4
  )
  methodData %>% ggplot() + 
    geom_point(aes(x=value_x, y=value_y)) + 
    ggtitle(glue("{title}\tr2={r2_value}"))
}
#' Lineplot of calibrators/norm and alignments.
#' 
#' @param table digested_adat table
#' @param target input is calibrators or raw/normalized RFU table
#' @param log2 do log2 transformation of RFU values
#' @param title main plot title
#' @param x_lab label for x axis
#' @param sample_f subset table to build dist_plot norm def. 0.2 (float between 0 and 1)
#' @param log write to log file (def.TRUE)
#' @return ggplot object.
#' @examples
#' plot_distline(digested_adat$calibrators, target="calibrators", log2=TRUE, title="calibrators", x_lab="RFU")
#' @export
#' @import dplyr
#' @import ggplot2
#' @importFrom forcats fct_reorder
plot_distline <- function(
  table, 
  target='calibrators', 
  log2=FALSE, 
  title=target, 
  x_lab=NULL,
  aptamer_tag='seq.',
  sample_f = 0.2,
  log=TRUE
  ) { 
  # reassign column names
  if (target=='calibrators') {
    table %>% rename(Samples = Barcode) -> data
  } else if (target=='norm') {
    table %>% 
    rename(Samples = SampleId) %>% 
    sample_frac(sample_f) -> data 
  } else {
    if (log) {
      log_info("unknown value for target:\t{target}")
    }
    return(NULL)
  }
  # pivot longer
  data %>% 
    select(
      starts_with(aptamer_tag), PlateId, Samples
    ) %>% 
    pivot_longer(!c("PlateId", "Samples")
    ) -> data_long
  if (log2) {
    data_long  %>%
      mutate(RFU=log2(value)) -> data_long
    if (is.null(x_lab)) {
      x_lab = glue("log2(RFU)")
    }
    else {
      x_lab = glue("log2({x_lab})")
    }
  } else {
    data_long  %>%
      mutate(RFU=value) -> data_long
    if (is.null(x_lab)) {
      x_lab = glue("RFU")
    }
  }
  data_long %>% 
    mutate(
      Samples=fct_reorder(
        Samples, desc(PlateId))) -> data_long
  # generate plot
  data_long %>%
    ggplot(aes(y=Samples, x=RFU)) + 
    geom_line(aes(color=PlateId))  +
    theme_minimal() +
    labs(
      x = glue("{x_lab} [red dot: median]"),
      y = "Samples",
      title = glue("Alignment of {title}")
    ) + theme(
      axis.ticks.y=element_blank(),
      axis.text.y=element_blank(),
      legend.position="right"
    ) + stat_summary(
      fun=median,
      geom="point",
      shape=20,
      size=2,
      color="red",
      fill="red"
    ) -> plot
  return(plot)
}