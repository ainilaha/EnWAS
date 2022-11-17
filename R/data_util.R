
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




