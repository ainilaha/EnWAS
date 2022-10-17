# library(DBI)
# nhanes_db <- DBI::dbConnect(RSQLite::SQLite(), "C:\\projects\\data_nhanes\\nhanes_new.sqlite")
# dbListTables(nhanes_db)
#
#
# base_df <- dbGetQuery(nhanes_db, "SELECT demo.SEQN,
#                                     (BPXDI1+BPXDI2)/2 AS DIASTOLIC,RIAGENDR,RIDAGEYR,RIDRETH1,BMXBMI
#                                   FROM
#                                     demo
#                                   INNER JOIN body_measures ON demo.SEQN=body_measures.SEQN
#                                   INNER JOIN blood_pressure ON demo.SEQN=blood_pressure.SEQN
#                                   WHERE
#                                       RIDAGEYR>20
#                                     AND BMXBMI is not NULL
#                                     AND BPXDI1 IS NOT NULL and BPXDI1 <> 0
#                                     AND BPXDI2 IS NOT NULL and BPXDI2 <> 0")


# exposure_cols <- paste(exposure_vars, collapse = ", ")

#' Title
#'
#' @param exprs
#' @param tb_name
#'
#' @return
#' @export
#'
#' @examples
load_exprs<-function(exprs,tb_name,nhanes_db=nhanes_db){
  exprs_cols <- paste(exprs, collapse = ", ")
  sql = paste0("SELECT demo.SEQN,",exprs_cols,"
               FROM ",
               tb_name,
               " INNER JOIN demo ON ",tb_name,".SEQN=demo.SEQN
               WHERE
                    RIDAGEYR>20
               ")
  expr_data <- dbGetQuery(nhanes_db,sql)
  expr_data <- phesant(expr_data)
  expr_data
}


# diet_expr <- load_exprs(exposure_cols,"DietaryInterviewTotalNutrientIntakesFirstDay")
#
# data <- merge(base_df,diet_expr$data,by="SEQN")



# lm_str <- 'DIASTOLIC ~ RIDAGEYR*RIAGENDR + BMXBMI + RIDRETH1'
# lm_base <- glm(formula = as.formula(lm_str), family = gaussian,data,na.action=na.omit)



