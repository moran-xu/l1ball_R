require("ggplot2")
require("reshape2")

plotBoxplot<- function(fit,p_show = 50, outlier =FALSE){
  
  df <- melt(fit$theta[,1:p_show])
  colnames(df)<- c("name","index","value")
  df$index<- as.factor(df$index)
  # Basic boxplot
  df['y']<- rep(y[1:p_show], each=nrow(df)/p_show)
  
  p<-   ggplot(df, aes(x=index, y=value))
  if(!outlier){
    p+ geom_boxplot(outlier.shape = NA)
  }else{
    p  + geom_boxplot()
  }
  #+geom_point(aes(x=index,y=y),col='red')
}