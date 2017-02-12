#' SimBCarray.R Simulate Array Data
#'
#' The functions simulates array data under a block diagonal
#' covariance matrix and a simple linear model for the sample
#' specific means.
#'
#' @details add some details here
#'
#' @param SimBCdes A SimBCdes class object which contains:\cr
#'           \itemize{
#'           \item{N:} number of samples\cr
#'           \item{Pblock:} number of feature blocks\cr
#'           \item{BlockCovDF:} Pblock length DF\cr
#'            \itemize{
#'             \item{p -} number of features in block\cr
#'             \item{rho -} correlation\cr
#'             \item{scale -} scaling factor for covariance matrix\cr
#'             }
#'           \item{DesignMat:} N x Q  design matrix [Non-Full Rank for convenience]\cr
#'             N samples and Q covariate design columns\cr
#'           \item{BlockBetaMat:} Q x Pblock  matrix of model coefficients\cr
#'               valid per block.\cr
#'           \item{DesignDF:} design dataframe with Q rows and columns:\cr
#'             \itemize{
#'                \item{name -} variable name\cr
#'                \item{type -} factor of numeric\cr
#'                }
#'           \item{SampleNames:} N length character vector for the sample names\cr
#'           \item{FeatNames:} P length character vector for the feature names\cr
#'           }
#'
#' @return A list containing:\cr
#'         \itemize{
#'           \item{X:} P x N data matrix
#'           \item{PhenoDF:} N x Qunq data frame, where Qunq is the number\cr
#'              of unique variable names in DesignDF
#'               }
#'
#' @export
SimBCarray=function(SimBCdes){

    SimBlock=function(i){
        Sigma=array(rep(SimBCdes$BlockCovDF$rho[i],SimBCdes$BlockCovDF$p[i]^2),
                    dim=c(SimBCdes$BlockCovDF$p[i],SimBCdes$BlockCovDF$p[i]))
        diag(Sigma)=1
        Sigma=SimBCdes$BlockCovDF$scale[i]*Sigma
        mu=SimBCdes$DesignMat%*%SimBCdes$BlockBetaMat[,i]

        mvfun=function(j) mvrnorm(1,mu=rep(mu[j],SimBCdes$BlockCovDF$p[i]),Sigma=Sigma)

        return(mapply(mvfun, 1:SimBCdes$N))  }


    tdat=mapply(SimBlock,1:SimBCdes$Pblock)

    if(SimBCdes$Pblock==1) X=tdat[[1]]

    if(SimBCdes$Pblock>1) X=rbind(tdat[[1]],tdat[[2]])
    if(SimBCdes$Pblock>2){
        for(j in 3:(SimBCdes$Pblock)){
            X=rbind(X,(tdat[[j]]))
        }}


    colnames(X)=SimBCdes$SampleNames
    rownames(X)=SimBCdes$FeatNames

    # now, populate PhenoDF
    unqF=as.character(unique(SimBCdes$DesignDF$name))
    nunqF=length(unqF)
    Ftype=rep("",nunqF)
    PhenoDF=array("",dim=c(SimBCdes$N,nunqF))
    for(i in 1:nunqF){
        Fdx=which(SimBCdes$DesignDF$name==unqF[i])
        if(SimBCdes$DesignDF$type[Fdx[1]]=="factor"){
                 PhenoDF[,i]=(colnames(SimBCdes$DesignMat)[Fdx])[
                    apply(SimBCdes$DesignMat[,Fdx], 1, function(x) which(x==1))
            ]
                 Ftype[i]="factor"
        } # end "factor"
        if(SimBCdes$DesignDF$type[Fdx[1]]=="numeric"){
         PhenoDF[,i]=as.character(SimBCdes$DesignMat[,Fdx])
         Ftype[i]="numeric"
        } # end "numeric"
    }# end loop over unque covariates

    PhenoDF=as.data.frame(PhenoDF)
    for(i in which(Ftype=="numeric")) PhenoDF[,i]=as.numeric(PhenoDF[,i])

    colnames(PhenoDF)=unqF
    rownames(PhenoDF)=colnames(X)

    BCarray=list(X=X,PhenoDF=PhenoDF)
    class(BCarray)="BCarray"
    return(BCarray)
}

