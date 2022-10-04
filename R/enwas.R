

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

    # model_list <- vector(mode = "list", length = num_var)
    num_cols <- sapply(data_set[, exposure_vars], is.numeric)
    num_cols <- names(num_cols[num_cols == TRUE])
    if (inv_norm) {
      data_set[, num_cols] <- lapply(data_set[, num_cols], invNorm)
    }

    qc_cols <- c('terms',"null.deviance","df.null","logLik",
                 "AIC","BIC","deviance","df.residual", "nobs")
    qc_mtx <- matrix(0, nrow = num_var, ncol = length(qc_cols))
    colnames(qc_mtx) <- qc_cols
    qc_mtx[,1] <- exposure_vars

    lrt_cols <- c('terms',"Resid. Dev","Df","Deviance","Ratio","Pr(>Chi)",
                  "omega","p_omega","LRTstat","p_LRT")
    lrt_mtx <- matrix(0, nrow = num_var, ncol = length(lrt_cols))
    colnames(lrt_mtx) <- lrt_cols
    lrt_mtx[,1] <- exposure_vars

    base_m <- glm(formula = as.formula(base_model),data_set,family = gaussian())
    association_list <- data.frame()
    # factor_terms <- c() # to hold the factors terms
    for (i in 1:num_var) {
      exposure <- exposure_vars[i]

      # if (is.factor(data_set[, exposure])) {
      #   factor_terms <-
      #     c(factor_terms, paste0(exposure, levels(data_set[, exposure])))
      # }
      # model_list[[i]] <- mod

      model <- build_formula(base_model, exposure)
      mod <- glm(model, data_set,family = gaussian())
      mod_df <- broom::tidy(mod)
      association_list <- rbind(association_list, mod_df)

      qc_mtx[i,2:9] <- round(unlist(broom::glance(mod)),3)
      #-----------------ANOVA LRT---------------------------
      ano_res <- anova(base_m,mod,test="LRT")
      lrt_mtx[i,2:4] <- round(unlist(ano_res[2,c("Resid. Dev","Df","Deviance")]),3)
      lrt_mtx[i,5] <- paste0(round(ano_res$"Deviance"[2]/ano_res$"Resid. Dev"[2]*100,3),"%")
      p_value <- ano_res$"Pr(>Chi"[2]
      lrt_mtx[i,6] <- if(is.na(p_value) | p_value<1e-4) "<1e-4" else round(p_value,4)

      #-----------------vuongtest---------------------------
      vong <- nonnest2::vuongtest(base_m,mod,nested = TRUE)
      lrt_mtx[i,7:9] <- round(unlist(vong[c("omega","p_omega","LRTstat")]),3)
      p_value <- vong$p_LRT$A
      lrt_mtx[i,10] <- if(is.na(p_value) | p_value<1e-3) "<1e-4" else round(p_value,4)


    }


    # xwas_list <-
    #   association_list[association_list$term %in% c(exposure_vars, factor_terms), ]

    xwas_list <-
      association_list[association_list$term %in% exposure_vars, ]


    if(inv_norm==FALSE){
      sd_x_list <-  sapply(data_set, function(x)
        sd(as.numeric(x)))
      for (var in exposure_vars) {
        xwas_list[grepl(var, xwas_list$term), c('estimate', 'std.error')] <-
          xwas_list[grepl(var, xwas_list$term), c('estimate', 'std.error')] * sd_x_list[var]
      }

    }



    xwas_list$fdr <-
      p.adjust(xwas_list$p.value, method = 'BY')

    xwas_list$lower <- xwas_list$estimate - 1.96 * xwas_list$std.error
    xwas_list$upper <- xwas_list$estimate + 1.96 * xwas_list$std.error

    qc_mtx <- as.data.frame(qc_mtx)
    qc_mtx <- qc_mtx[qc_mtx$terms %in%xwas_list$term,]

    lrt_mtx <- as.data.frame(lrt_mtx)
    lrt_mtx <- lrt_mtx[lrt_mtx$terms %in%xwas_list$term,]

    return (list(qc_mtx = qc_mtx,lrt_mtx = lrt_mtx, enwas_res = xwas_list))


  }



#' Environment‐Wide Association Study (EnWAS) with testing
#'
#' @param base_model a string of base model
#' @param exposure_vars the phenotype list
#' @param train_set taring set
#' @param test_set test set
#' @param lab_col label column
#' @param inv_norm
#'
#' @return the model list and EnWAS result
#' @export
#'
#' @examples linear_model <- 'BMXWAIST ~ RIDAGEYR*RIAGENDR + BMXBMI'
#' linear_res <- enwas_cv(linear_model, exposure_vars, train_set,test_set)
enwas_cv <-
  function(base_model,
           exposure_vars,
           train_set,
           test_set,
           lab_col="diastolic",
           inv_norm = FALSE) {

    num_var <- length(exposure_vars)

    model_list <- vector(mode = "list", length = num_var)
    names(model_list) <- exposure_vars
    num_cols <- sapply(train_set[, exposure_vars], is.numeric)
    num_cols <- names(num_cols[num_cols == TRUE])
    if (inv_norm) {
      train_set[, num_cols] <- lapply(train_set[, num_cols], invNorm)
      test_set[, num_cols] <- lapply(test_set[, num_cols], invNorm)
    }

    mse <- rep(NA,num_var)

    for (i in 1:num_var) {
      model_str <- build_formula(base_model, exposure_vars[i])
      model <- lm(model_str, train_set)
      pred_vals <- predict(model, test_set)
      mse[i] <- mean((pred_vals - test_set[,lab_col]) ^ 2)
    }

    mse

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

#' Inverse Normal Transformation
#'
#' @param x
#'
#' @return transformed data
#' @export
#'
#' @examples invNorm(nhanes$BMXWAIST)
invNorm <- function(x) {
  qnorm((rank(x) - 3/8)/(length(x) +1 - 6/8))
  }



