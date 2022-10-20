#' Cross Validation for Base Models
#'
#' @param model_list a list of strings for the models with names
#' @param label outcomes
#' @param group_col columns for group data in CV
#' @param df data set
#'
#' @return MSE matrix
#' @export
#'
#' @examples
#' lm_str <- "diastolic ~ RIDAGEYR*RIAGENDR + BMXBMI + RIDRETH1"
#' ns_str <- "diastolic ~ ns(RIDAGEYR, knots = seq(30, 80, by = 10), Boundary.knots=c(20,90)) * RIAGENDR + ns(BMXBMI,knots = c(seq(15, 45, by = 5),seq(45,65,by=10)),Boundary.knots=c(10,85)) + RIDRETH1"
#' model_list <- c(lm_str,ns_str)
#' names(model_list) <- c("linear","spline")
#' mse_mtx <- cv_base_m(model_list,label="diastolic",group_col="years",df=nhanes)
cv_base_m <- function(model_list,label,group_col,df=nhanes) {

  groups <- levels(df[,group_col])

  len_group <- length(groups)

  mse_mtx <- matrix(0, nrow = len_group, ncol = length(model_list))
  rownames(mse_mtx) <- groups
  colnames(mse_mtx) <- names(model_list)

  for (i in 1:len_group) {
    for(j in 1:length(model_list)){
      cv_train <- df[df[,group_col] != groups[i], ]
      cv_test <- df[df[,group_col] == groups[i], ]
      ns_model <- lm(as.formula(model_list[j]), cv_train)
      pred_vals <- predict(ns_model, cv_test)
      mse_mtx[i,j] <- mean((pred_vals - cv_test[,label]) ^ 2)

    }
  }


  mse_mtx
}

#' Produce Table of ANOVA LRT by iterative over the terms
#'
#' @param form_str string of the formula
#' @param data data set
#'
#' @return results as data frame
#' @export
#'
#' @examples lm_str <- 'diastolic ~ RIDAGEYR*RIAGENDR + BMXBMI + RIDRETH1'
#' anova_lrt(lm_str)
anova_lrt <- function(form_str, data=nhanes){
  full_model <- glm(formula = as.formula(form_str),data,family = gaussian())

  form_str <- gsub(" ", "", form_str, fixed = TRUE)
  terms <- gsub("[+~*]", " ", form_str)
  terms <- c(unlist(strsplit(terms," ")),"Full Model")
  cols <- c("Terms",'RSS','Df',"Sum of Sq","Ratio","Pr(>Chi)")
  n_terms <- length(terms)-1
  lrt_mtx <- matrix(0, nrow = n_terms, ncol = length(cols))
  colnames(lrt_mtx) <- cols
  lrt_mtx[,1]<-ifelse(nchar(terms[-1]) > 13, paste0(substring(terms[-1], 1, 10), "..."), terms[-1])

  for (i in 1:n_terms){

    remove_term <- if (i+1==2) paste0(terms[i+1],"[+~*]") else paste0("[+~*]",terms[i+1])
    remove_term <- gsub("\\(","\\\\(",remove_term)
    remove_term <- gsub("\\)","\\\\)",remove_term)
    sub_str <- gsub(remove_term, "", form_str)
    sub_model <- glm(formula = as.formula(sub_str),data,family = gaussian())

    ano_res <- anova(sub_model,full_model,test="LRT")
    lrt_mtx[i,2] <- round(ano_res$"Resid. Dev"[1],3)
    lrt_mtx[i,3] <- round(ano_res$Df[2],3)
    lrt_mtx[i,4] <- round(ano_res$"Deviance"[2],3)
    lrt_mtx[i,5] <- paste0(round(ano_res$"Deviance"[2]/ano_res$"Resid. Dev"[2]*100,3),"%")
    p_value <- ano_res$"Pr(>Chi)"[2]
    lrt_mtx[i,6] <- if(is.na(p_value) | p_value<1e-3) "<1e-3" else round(p_value,3)

  }

  as.data.frame(lrt_mtx)

}

