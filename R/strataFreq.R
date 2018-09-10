strataFreq <- function(g) {
  g@data %>% 
    dplyr::group_by(stratum) %>% 
    dplyr::summarize(num.ind = n_distinct(id)) %>% 
    dplyr::ungroup()
}