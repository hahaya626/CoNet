#----------CPNT----------#
#' Network regression analysis for survival outcomes in transcriptome-wide association studies
#'
#' @param data1 a data.frame in which to represent  the relationship among nodes:if one node links another, the relationship between them is 1, otherwise is 0
#' @param data2 a data.frame in which includes all variables needed in formula
#' @param num is the number of nodes
#' @param formula a formula to be passed to function coxph(). For more details, please refer to package survival.
#'
#' @return a data.frame
#' @export
#'
#' @examples
#' library(survival)
#' data("data1")
#' data("data2")
#' num <- dim(data1[1])
#' formula <-Surv(ctime,censor)~V1+V2+V3+V4+V5+V6+V7+V1_2+V2_3+V3_4+V3_5+V3_6+V3_7+PC1+PC2+PC3+PC4+PC5
#' result <- CPNT(data1, data2, num, formula)
CPNT <- function(data1, data2, num, formula){
  #----------EDGE----------#
  b<-unlist(data1)
  sigma<-matrix(b,nrow = num)
  sigma0<-sigma
  sigma0[lower.tri(sigma0,diag=T)]<- 0
  EDGE <- which(sigma0!=0,arr.ind=T)
  #----------PM----------#
  MDIPreData<-apply(EDGE,1,MDIF,data=DataSet_)
  colname<-apply(EDGE,1,RENAME,varname="V")
  colnames(MDIPreData)<-colname
  MDIFitData<-cbind(DataSet_,MDIPreData)
  MDIFitData<-data.frame(MDIFitData)
  #----------regression----------#
  MDIfit<-coxph(formula,data=MDIFitData)
  MDItemp<-summary(MDIfit)
  MDItem<-MDItemp$coefficients
  return(MDItem)
}
