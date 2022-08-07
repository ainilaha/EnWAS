

#' Environment‐Wide Association Study (EnWAS)
#'
#' @param base_model a string of base model
#' @param exposure_vars the phenotype list
#' @param data_set data set
#' @param inv_norm
#'
#' @return the model list and EnWAS result
#' @export
#'
#' @examples linear_model <- 'BMXWAIST ~ RIDAGEYR*RIAGENDR + BMXBMI'
#' linear_res <- enwas(linear_model, exposure_vars, data)
enwas <-
  function(base_model,
           exposure_vars,
           data_set,
           inv_norm = FALSE) {
    num_var <- length(exposure_vars)

    model_list <- vector(mode = "list", length = num_var)
    num_cols <- sapply(data_set[, exposure_vars], is.numeric)
    num_cols <- names(num_cols[num_cols == TRUE])
    if (inv_norm) {
      data_set[, num_cols] <- lapply(data_set[, num_cols], invNorm)
    }




    association_list <- data.frame()
    factor_terms <- c() # to hold the factors terms
    num_var <- length(exposure_vars)
    for (i in 1:num_var) {
      exposure <- exposure_vars[i]

      if (is.factor(data_set[, exposure])) {
        factor_terms <-
          c(factor_terms, paste0(exposure, levels(data_set[, exposure])))
      }

      model <- build_formula(base_model, exposure)

      mod <- lm(model, data_set)
      model_list[[i]] <- mod
      mod_df <- broom::tidy(mod)
      association_list <-
        rbind(association_list, mod_df) # step 3b of algorithm
    }


    xwas_list <-
      association_list[association_list$term %in% c(exposure_vars, factor_terms), ]


    # sd_x_list <-  sapply(data_set[,num_cols],sd)
    sd_x_list <-  sapply(data_set, function(x)
      sd(as.numeric(x)))
    for (var in exposure_vars) {
      xwas_list[grepl(var, xwas_list$term), c('estimate', 'std.error')] <-
        xwas_list[grepl(var, xwas_list$term), c('estimate', 'std.error')] * sd_x_list[var]
    }



    xwas_list$fdr <-
      p.adjust(xwas_list$p.value, method = 'BY')



    xwas_list <- xwas_list %>%
      mutate(lower = estimate - 1.96 * std.error)  %>%
      mutate(upper = estimate + 1.96 * std.error)

    return (list(model_list = model_list, enwas_res = xwas_list))


  }

#' Build formula for EnWAS
#'
#' @param base_model the formula of the model with the covariates
#' @param exposure_var exposure variables (aka. phenotype lists or predictors)
#' @param inv Boolean flag for use inverse normal transformation
#'
#' @return formula for EnWAS
#' @export
#'
#' @examples build_formula("BMXWAIST ~ RIDAGEYR*RIAGENDR + BMXBMI", phenotype_lists)

build_formula <- function(base_model,exposure_var,inv=FALSE) {
  var <- "+ %s"
  if(inv) var <- " + invNorm(%s)"

  as.formula(
    sprintf(paste0(base_model, var),
            exposure_var
    )
  )
}


