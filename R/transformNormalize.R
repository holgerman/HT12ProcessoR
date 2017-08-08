transformNormalize <- function(total_nobkgd_eset_goodind, methodtransform = "quantile", dolog2 = T) {
    
    
    if (methodtransform == "quantile") {
        
        time2 <- Sys.time()
        total_nobkgd_eset_ql <- affyPLM::normalize.ExpressionSet.quantiles(total_nobkgd_eset_goodind, transfn = "none")
        exprs(total_nobkgd_eset_ql)[1:11, 1:6]
        if (dolog2) 
            exprs(total_nobkgd_eset_ql) = log2(exprs(total_nobkgd_eset_ql))
        exprs(total_nobkgd_eset_ql)[1:11, 1:6]
        Sys.time() - time2  #foro 19
        
        total_nobkgd_eset_ql
        
    } else if (methodtransform == "rsn") 
      {
        
        time2 <- Sys.time()
        total_nobkgd_eset_ql <- total_nobkgd_eset_goodind
        exprs(total_nobkgd_eset_ql)[1:11, 1:6]
        if (dolog2) 
            exprs(total_nobkgd_eset_ql) = log2(exprs(total_nobkgd_eset_ql))
        exprs(total_nobkgd_eset_ql)[1:11, 1:6]
        
        total_nobkgd_eset_ql = lumi::lumiN(total_nobkgd_eset_ql, method = "rsn", ifPlot = TRUE)
        exprs(total_nobkgd_eset_ql)[1:11, 1:6]
        
        total_nobkgd_eset_ql
        Sys.time() - time2  #arnor 30s
        
    } else stop("Parameter 'normalisation_method' must be either 'quantile' or 'rsn'.")
    
    
   
    
    total_nobkgd_eset_ql
}
