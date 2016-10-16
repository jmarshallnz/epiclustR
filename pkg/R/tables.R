#' Table of outbreaks from an epiclustR model
#' 
#' @param mod The model to table
#' @param data The data used for the model
#' @param period The time period to use. Defaults to NULL which is all time
#' @param threshold The threshold for an outbreak. Only outbreaks that exceed this threshold will be displayed. Defaults to 0.2
#' @export
table_outbreaks <- function(mod, data, period = NULL, threshold = 0.2) {
  # pull out the outbreaks
  X <- ssapply(mod, extract_variable, 'X')
  mX = apply(X, 1:2, mean)
  mX = data.frame(mX, check.names = FALSE)
  mX$ReportWeek = as.Date(row.names(data$cases))
  regions = tidyr::gather_(mX, key_col='Region', value_col='P', gather_cols=setdiff(names(mX), 'ReportWeek'), convert=TRUE)

  # match these up with cases
  case_outbreaks <- data$case_list %>% dplyr::left_join(data$spat_list, by="Spatial") %>%
    dplyr::left_join(regions)
  case_outbreaks <- case_outbreaks[,c('CaseID', 'ReportWeek', 'Region', 'P')]

  # filter out the date period
  if (!is.null(period)) {
    # TODO: Period should be date-ized
    case_outbreaks <- case_outbreaks[case_outbreaks$ReportWeek %in% period]
  }

  # filter out based on the threshold
  case_outbreaks <- case_outbreaks[case_outbreaks$P >= threshold,]

  # pull out the cases that correspond to these outbreaks
  tab <- case_outbreaks[order(case_outbreaks$ReportWeek),]
  dupes <- duplicated(tab[,c('Region', 'ReportWeek', 'P')])
  tab[dupes, c('Region', 'ReportWeek', 'P')] <- ''

  tab[, c('CaseID', 'ReportWeek', 'Region', 'P')]
}
