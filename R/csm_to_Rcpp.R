csm_to_Rcpp <- function(x){
  if(is.function(x)){
    x_args <- args(x) |>
      names() |>
      head(-1)
    x_body <-
      body(x) |>
      as.list() |>
      tail(-1) |>
      lapply(csm_to_Rcpp)
  # }else if("{" %in% class(x)){
  #   as.list(x) |>
  #     lapply(csm_to_Rcpp)
  }else if(is.expression(x)){

  }else if(is.name(x)){
    return(x)
  }
}
