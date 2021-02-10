#######Utility Functions Written by Enrico Crema#########

#Compute Great Arc Distances
distMat<-function(data)
  {
require(argosfilter)
res<-matrix(0,nrow=dim(data)[1],ncol=dim(data)[1])
for (i in 1:dim(data)[1])
  {
    for (j in 1:dim(data)[1])
      {
        if(i!=j){
        res[i,j]<-distance(lat1=data$y[i],lat2=data$y[j],lon1=data$x[i],lon2=data$x[j])
      }
      }
  }
res[is.nan(res)]=0
return(as.dist(res))
  }

BinaryCulture<-function(cultureList) #written by Crema 
    {
res<-matrix(NA,ncol=length(cultureList),nrow=length(cultureList))

for (i in 1:length(cultureList))
    {
        for (j in 1:length(cultureList))
            {
                if (i>j)
                    {
                        res[i,j]=as.numeric(cultureList[i]!=cultureList[j])
                    }
            }
    }
return(res)
}

#utility function by Darwin https://stackoverflow.com/a/28845828
write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}