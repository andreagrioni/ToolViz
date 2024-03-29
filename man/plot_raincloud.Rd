% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plot_raincloud}
\alias{plot_raincloud}
\title{generate raincloud plot aptamers
inspired by https://www.cedricscherer.com/2021/06/06/visualizing-distributions-with-raincloud-plots-and-how-to-create-them-with-ggplot2/}
\usage{
plot_raincloud(
  data,
  x = NULL,
  y = NULL,
  target = NULL,
  annotation = NULL,
  group_ids = NULL,
  title = "",
  x_lab = x,
  y_lab = y,
  x_rot = 90,
  spaghetti = NULL,
  colspaghetti = NULL,
  stats = FALSE
)
}
\arguments{
\item{data}{tibble class object as samplexfeatures.}

\item{x}{x-axis value (def. aptamers).}

\item{y}{y-axis value (def. rfu).}

\item{target}{EntrezGeneSymbol or Aptamer SeqId.}

\item{annotation}{Somascan Annotation Table.}

\item{group_ids}{TODO.}

\item{title}{main plot title.}

\item{x_lab}{x axis label.}

\item{y_lab}{y axis label.}

\item{x_rot}{rotation x axis label (def. 90).}

\item{spaghetti}{group samples by feature (allow spaghetti plot).}

\item{colspaghetti}{color spaghetti line by feature (def. NULL).}

\item{stats}{add stat_summary (median) (bool, def. FALSE).}
}
\value{
gglot object
}
\description{
generate raincloud plot aptamers
inspired by https://www.cedricscherer.com/2021/06/06/visualizing-distributions-with-raincloud-plots-and-how-to-create-them-with-ggplot2/
}
\examples{
raincloud(data=table,
annotation=anno, target='CXCL10',
x='time_idx',  group_ids="ARM", x_rot=90,
spaghetti='SUBJID', stats=TRUE)
}
