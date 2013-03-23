
######################################################bring in the data we need
setwd('/media/jenya/PATRIOT/project')
coef = read.table("coef.csv", header = TRUE)
cov.vals = read.table("cov.vals.csv", header = TRUE)

library(nlme)
library(contrast)
library(parallel)
##################### Extractor Functions from Jose Pinheiro ##############
#
### These two functions (varRan and varWithin) extract the random effects for a
### given level (varRan) or the level-1 error (varWithin)
varRan <- function(object, level = 1)
  ### NOTE: level = 1 is the largest grouping level (participant) and
  ### level = 2 is the next largest (day)
{
  sigE <- object$sig^2
  sigE*pdMatrix(object$modelStruct$reStruct)[[level]]
}
varWithin <- function(object, wch)
  ### NOTE: wch specifies "which" individual's level-1 error term you wish to
  ### extract. For the models we considered, these will be similar across
  ### individuals but can vary if the level-1 error term is extended
{
  getMat <- function(wch, grps, ugrps, cS, corM, vF, std) {
    wgrps <- grps[grps == ugrps[wch]]
    nG <- length(wgrps)
    if (!is.null(cS)) {
      val <- corM[[wch]]
    } else {
      val <- diag(nG)
    }
    if (!is.null(vF)) {
      std <- diag(std[[wch]])
      val <- std %*% (val %*% std)
    } else {
      val <- std*std*val
    }
    val
  }
  grps <- getGroups(object)
  if (is.null(grps)) { # gls/gnls case
    grps <- rep(1, length(resid(object)))
    ugrps <- "1"
  } else {
    ugrps <- unique(as.character(grps))
  }
  if (!is.null(vF <- object$modelStruct$varStruct)) {
    ## weights from variance function, if present
    std <- split(object$sigma/varWeights(vF), grps)[ugrps]
  } else {
    std <- object$sigma
  }
  if (!is.null(cS <- object$modelStruct$corStruct)) {
    corM <- corMatrix(cS)[ugrps]
  } else {
    corM <- NULL
  }
  if (!missing(wch)) { # particular group
    return(getMat(wch, grps, ugrps, cS, corM, vF, std))
  } else {
    if (length(ugrps) == 1) { # single group
      return(getMat(1, grps, ugrps, cS, corM, vF, std))
    }
    val <- vector("list", length(ugrps))
    names(val) <- ugrps
    for(i in 1:length(ugrps)) {
      val[[i]] <- getMat(i, grps, ugrps, cS, corM, vF, std)
    }
    return(val)
  }
}





################################arguments to simulateData function if the model was present

simulateData(N = 30,
            beta = c(fixef(model)),
            sigma1 = sqrt(varWithin(model, wch=1)[1,1]),
            sigma2 = sqrt(varRan(model, level=2)),
            tau = varRan(model, level=1))







simulateData = function(N,  beta,sigma1, sigma2, tau){
  
  ###subject-level random-effects
  bj = rnorm(N, mean = 0, sd=sqrt(tau))
  b0j = rep(bj, each = 24*8) # vector of length N*24*8
  
  # vector of length N*8
  ### day-level random intercept
  epsilon2 = rnorm(N*8, mean = 0, sd = sigma2)# vector of length N*8
  epsilon2 = rep(epsilon2, each = 24) # vector of length N*24*8
  
  ###level-1 error
  epsilon1 = rnorm(N*8*24, mean = 0, sd = sigma1)# vector of length N*24*8
  
  
  ######################################################################create treat "one"
  
  t=data.frame(value = beta[1,2]*cov.vals$cov1 + 
                 beta[2,2]*cov.vals$cov2 + 
                 rep(beta[3:26,2],times =8))           
  
  one = data.frame(value = rep(t$value,times = N) + b0j  + epsilon2 + epsilon1)
  one$treat = factor("one")
  one$id = factor(rep(1:N, each = 24*8))
  #####################################################################create treat "two"
  
  t.2 = data.frame(value  = t$value + 
                     rep(beta[28,2], times = 24*8) + 
                     rep(beta[53:76,2],times = 8))
  
  two = data.frame(value = rep(t.2$value, times = N) + b0j + epsilon2 + epsilon1)
  two$treat = factor("two")
  two$id = factor(rep(1:N, each = 24*8))
  #####################################################################create treat "three"
  
  t.3 = data.frame(value = t$value + 
                     rep(beta[27,2], times = 24*8) + 
                     rep(beta[29:52,2],times = 8))
  
  three = data.frame(value = rep(t.3$value, times = N) + b0j + epsilon2 + epsilon1)
  three$treat = factor("three")
  three$id = factor(rep(1:N, each = 24*8))
  simulatedData = data.frame(rbind(one, two, three))
  simulatedData$hour = factor(c("01","02","03","04","05","06","07","08","09","10",
                                "11","12","13","14","15","16","17","18","19","20",
                                "21","22","23","24"))
  simulatedData$Day = factor(rep(1:8, each = 24))
  
  
  
  simulatedData$cov1 = cov.vals$cov1
  simulatedData$cov2 = cov.vals$cov2
  simulatedData$value[simulatedData$value < 0] = NA
  
  return(simulatedData)
}



######################################################computations in parallel

RNGkind("L'Ecuyer-CMRG")


cl = makeCluster(detectCores(), type = "FORK")

clusterSetRNGStream(cl, 1985)


########################check what's inhereted on each node
clusterEvalQ(cl, ls())


system.time( clusterApplyLB(cl, 1:100,
                         
                         function(N=25, beta = coef, 
                                   sigma1= 0.924259217846522, 
                                   sigma2= 0.2636957, 
                                   tau = 0.3348744)
{ 
                           
                           data = simulateData(N,
                                               beta,
                                               sigma1,
                                               sigma2,
                                               tau)
  model = try(
    lme(value ~ cov1 + cov2 + hour*treat-1,
        data = data,
        na.action=na.omit,
        control= lmeControl(opt="optim"),
        random = ~1|id/Day, keep.data=TRUE))
  
  if(class(model)=="try-error" ||
       summary(model)$apVar[[1]] == "Non-positive definite approximate variance-covariance"){
    contrast.one.two = NA
  } else {
    
    contrast.one.two = contrast(model,
                                list(hour = c("01","02","03","04","05","06","07","08","09","10",
                                              "11","12","13","14","15","16","17","18","19","20",
                                              "21","22","23","24"),
                                     cov1 = 0, cov2 = 0,
                                     treat ="one"),
                                list(hour = c("01","02","03","04","05","06","07","08","09","10",
                                              "11","12","13","14","15","16","17","18","19","20",
                                              "21","22","23","24"),
                                     cov1 = 0, cov2 = 0,
                                     treat ="two"),
                                type="average",
                                X = TRUE )
    
  }})
)

stopCluster(cl)

mat = matrix(0, nrow = 200, ncol=2)

for (i in 1:8){
  
  mat[i,1] = as.numeric(results[[i]][1])
  mat[i,2]= as.numeric(results[[i]][7])*3

}







##############################allocate a matrix to save the output

mat = matrix(0, nrow = 100, ncol = 2)

#############################iterative function

system.time(
  for (i in 1:100){
    
    data = simulateData(N = 25,
                        beta = coef,
                        sigma1 = 0.924259217846522,
                        sigma2 = 0.2636957,
                        tau = 0.3348744)
    
    model = try(
      lme(value ~ cov1 + cov2 + hour*treat-1,
          data = data,
          na.action=na.omit,
          control= lmeControl(opt="optim"),
          random = ~1|id/Day, keep.data=TRUE))
    
    if(class(model)=="try-error" ||
         summary(model)$apVar[[1]] == "Non-positive definite approximate variance-covariance"){
      output.mat[i,1:2] == NA
    } else {
      
      contrast.one.two = contrast(model,
                                  list(hour = c("01","02","03","04","05","06","07","08","09","10",
                                                "11","12","13","14","15","16","17","18","19","20",
                                                "21","22","23","24"),
                                       cov1 = 0, cov2 = 0,
                                       treat ="one"),
                                  list(hour = c("01","02","03","04","05","06","07","08","09","10",
                                                "11","12","13","14","15","16","17","18","19","20",
                                                "21","22","23","24"),
                                       cov1 = 0, cov2 = 0,
                                       treat ="three"),
                                  type="average",
                                  X = TRUE )
      
      mat[i, 1] = round(contrast.one.two[[1]],5)
      mat[i, 2] = round(contrast.one.two[[7]]*3,10) 
      #####pvalue multiplied by 3 to correct for the # of contrasts
      
      
    }
    ### print the iteration number, so we know where we at
    
    iteration.number = i
    print(paste("Iteration #", iteration.number))
  }
  
)



df = data.frame(mat);names(df)=c("coef","pvalue")


sum(length(df$pvalue[df$pvalue <= 0.05]))/200

k=df[which(df$pvalue %in% NA),]











