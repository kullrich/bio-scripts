library(biogeo)
library(geodist)
library(readr)
library(dplyr)

get_pairwise_geodist<-function(file){
  options(scipen=22)
  input_df<-readr::read_delim(file, col_types = cols(.default = "c"))
  input_df$LAT_dec <- unlist(lapply(strsplit(input_df$LAT,","),function(x) 
    biogeo::dms2dd(
    as.numeric(x[1]),
    as.numeric(x[2]),
    as.numeric(x[3]),
    x[4])))
  input_df$LONG_dec <- unlist(lapply(strsplit(input_df$LONG,","),function(x) 
    biogeo::dms2dd(
    as.numeric(x[1]),
    as.numeric(x[2]),
    as.numeric(x[3]),
    x[4])))
  gdist_mat<-geodist::geodist(dplyr::select(input_df,LONG_dec,LAT_dec))
  colnames(gdist_mat)<-input_df$ID
  rownames(gdist_mat)<-input_df$ID
  return(gdist_mat)
}
