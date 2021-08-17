#' generate a quantile-quantile plot from limma fitted model.
#' 
#' @param model a limma class model.
#' @param coeff vector of coeffient names or indexes.
#' @return ggplot class object.
#' @examples
#' get_qqplot(fit2, c("PBOvsTRT"))
#' @export
get_qqplot <- function(model, coeff) {
  topTable(
    model,
    number=Inf,
    sort.by="none",
    coef=coeff
  ) %>%
    as_tibble -> wdata
  qq_ <- make_qqplot(
    -log10(wdata$P.Value),
    c(0.5, .75),
    glue("QQ Plot {coeff}"))
  return(qq_)
}

make_qqplot <- function(pvals, quant, title){  
  len = length(pvals)
  res=qqplot(-log10((1:len)/(1+len)),pvals,plot.it=F)
  plot(
    res$x,
    res$y, 
    main=title, 
    xlab="Theoretical", 
    ylab="Actual", 
    col=ifelse(res$y>as.numeric(quantile(res$y, quant[1])), 
               ifelse(res$y>as.numeric(quantile(res$y, quant[2])), 
                      "red", "blue"), "black")
  )
  abline(0, 1)
}
#' generate a bar plot of p-values.
#' 
#' @param model a limma class model.
#' @param coeff vector of coeffient names or indexes.
#' @return plot class object.
#' @examples
#' pval_distribution(fit2, c("PBOvsTRT"))
#' @export
pval_distribution <- function (model, coeff) {
  deg_genes <- topTable(
    model,
    number=Inf,
    sort.by="none",
    coef=coeff) # you need to avoid sorting so results will be in the same order as input.
  histogram <- hist(deg_genes$P.Value, main=glue("coeff {coeff}"))
  return(histogram)
}