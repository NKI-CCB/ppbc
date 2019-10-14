
## install.packages(c("DEoptim","snow","doParallel","doSNOW","foreach","doMC","sva"))


## minimisation function of the differential evolution algorithm ##
## there are 5 different metrics that can be used, pearson correlation, spearman correlation, euclidean distance, manhattan distance and rmse ##

min_function <- function(b,metric,cur_target,ref) {
  reg_data <- as.vector(as.matrix(ref) %*% b)
  if (metric %in% c("pearson","spearman")){
    return(1-cor(reg_data,cur_target,method=metric))
  }else if(metric == "rmse"){
    return(sqrt(mean((reg_data - cur_target)**2)))
  }else if(metric %in% c("euclidean","manhattan")){
    return(dist(rbind(reg_data,cur_target),method=metric)[1])
  }
  
}

smart.round <- function(x, digits = 0) {
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y / up
}


core_alg_decon <- function(i,target,ref,metric=metric){
  cur_target <- target[,i]
  n <- ncol(ref)
  mapfun <- function(x){smart.round(x,4)}
  
  ## TCD works with differential evolution to find the optimal combination of cell types ##
  ## DEoptim.control sets the variables for the differential evolution ##
  DEopt <- DEoptim.control(strategy=2,itermax=10000,NP=20,reltol=1e-8,steptol=1000,trace=F)
  ## Differntial evolution to minimize distance between the reference profile combination and the sample ##
  de_fit <- suppressWarnings(DEoptim(fn=min_function,lower=rep(0,n),upper=rep(10,n),metric=metric,cur_target=cur_target,ref=ref,control=DEopt,fnMap=mapfun))
  fit <- list()
  fit$value <- de_fit$optim$bestval
  fit$par <- de_fit$optim$bestmem
  fit$iter <- de_fit$optim$iter
  if (fit$iter < 10000){
    fit$convergence <- "yes"
    fit$message <- "change smaller than reltol (1e-08)"
  }else{
    fit$convergence <- "no"
    fit$message <- "maximum number of iterations reached"
  }
  #	fractions <- fit$par/sum(fit$par)
  fractions <- fit$par
  
  ## calculate some diagnostics on the similarity between reconstructed and real tumor ##
  cur_result <- as.matrix(ref) %*% fractions
  pearson_test <- cor.test(cur_target,cur_result,method="pearson")
  spearman_test <- suppressWarnings(cor.test(cur_target,cur_result,method="spearman"))
  rmse <- sqrt(mean((cur_target - cur_result)**2))
  linfit <- lm(cur_result ~ cur_target)
  linfit_force_origin <- lm(cur_result ~ 0 + cur_target)
  fit_results <- c(colnames(target)[i],fit$value,fit$convergence,fit$message,fit$iter,pearson_test$estimate,pearson_test$p.value,spearman_test$estimate,spearman_test$p.value,linfit$coefficients[2],linfit$coefficients[1],summary(linfit)$r.squared,linfit_force_origin$coefficients[1],summary(linfit_force_origin)$r.squared,rmse)
  return(c(fit_results,fractions))
}


## deconvolution wrapper function ##

deconvolve <- function(target=NULL,ref=NULL,sig_genes=NULL,rescale=F,metric="manhattan",do_parallel=F,no_cores=4,combat=F){
  if(!(metric %in% c("pearson","spearman","rmse","euclidean","manhattan"))){
    stop("metric has to be one of the following: pearson, spearman, rmse, euclidean, manhattan")
  }else if(is.null(target) | is.null(ref)){
    stop("target data and ref data have to be provided")
  }
  
  nsamples <- ncol(target)
  genes_in_both <- intersect(rownames(target),rownames(ref))
  if (is.null(sig_genes)){
    target <- target[genes_in_both,]
    ref <- ref[genes_in_both,]
    
  }else{
    if(!is.character(sig_genes) | !is.vector(sig_genes)){
      stop("sig_genes has to be a character vector")
    }else{
      sig_genes_in_both <- intersect(genes_in_both,sig_genes)
      if ( length(sig_genes_in_both) < 100){
        warning(paste("Only ",length(sig_genes_in_both),"signature genes overlapping between target data and ref data"))
      }
      target <- target[sig_genes_in_both,]
      ref <- ref[sig_genes_in_both,]
    }
  }
  if (rescale){
    ref <- apply(ref,2,function(x) (x/sum(x))*1000000)
    target <- apply(target,2,function(x) (x/sum(x))*1000000)
  }
  ## if needed combat can be used to decrease batch effects between the reference profiles and the sample sets ##
  ## the adjusted reference progiles are saved ##
  if (combat){
    target <- target[apply(target,1,function(x) sum(x>0)>0),]
    genes_in_both <- intersect(rownames(target),rownames(ref))
    target <- target[genes_in_both,]
    ref <- ref[genes_in_both,]
    ref <- apply(ref,2,function(x) (x/sum(x))*1000000)
    target <- apply(target,2,function(x) (x/sum(x))*1000000)	
    require(sva)
    t1 <- cbind(ref,target)
    t2 <- ComBat(t1,c(rep("ref",ncol(ref)),rep("target",ncol(target))),par.prior=T)
    #t2 <- apply(t2,2,function(x) (x/sum(x))*1000000)
    ref <- t2[,c(1:ncol(ref))]
    if (rescale){
      ref <- apply(ref,2,function(x) (x/sum(x))*1000000)
    }
    write.table(ref,paste("ref_combat_",metric,".csv",sep=""),sep="\t")
  }	
  n <- ncol(ref)
  
  if (do_parallel){
    library(foreach)
    library(doMC)
    library(snow)
    library(doSNOW)
    cl <- makeCluster(no_cores)
    registerDoSNOW(cl)
    
    fit_res <- t(foreach(i=c(1:nsamples), .combine= data.frame, .export=c("core_alg_decon","min_function","smart.round"), .packages=c("DEoptim","snow")) %dopar% {core_alg_decon(i,target=target,ref=ref,metric=metric)})
    colnames(fit_res) <- c("sample_name","distance","convergence_code","message","iterations","pearson_R",
                           "pearson_pval","spearman_R","spearman_pval","lin_regres_slope","lin_regres_intercept","lin_regres_r_squared",
                           "lin_regres_origin_slope","lin_regres_origin_r_squared","rmse",colnames(ref))
    stopCluster(cl)
    
  }else{
    
    fit_res <- as.data.frame(matrix(nrow=nsamples,ncol=(15+ncol(ref))))
    colnames(fit_res) <- c("sample_name","distance","convergence_code","message","iterations","pearson_R",
                           "pearson_pval","spearman_R","spearman_pval","lin_regres_slope","lin_regres_intercept","lin_regres_r_squared",
                           "lin_regres_origin_slope","lin_regres_origin_r_squared","rmse",colnames(ref))
    for ( i in 1:nsamples){
      
      fit_res[i,] <- core_alg_decon(i,target=target,ref=ref,metric=metric)
      
    }
  }
  reg_matrix <- fit_res[,c(16:ncol(fit_res))]
  reg_matrix <- apply(reg_matrix,2,as.numeric)
  reg_matrix <- t(apply(reg_matrix,1,function(x) x/sum(x)) )
  rownames(reg_matrix) <- fit_res[,1]
  fit_results <- as.data.frame(fit_res[,c(1:15)])
  
  return(list(cell_fractions = reg_matrix, fit_results=fit_results))
}



