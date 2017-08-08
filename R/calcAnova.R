calcAnova <- function(eset = total_nobkgd_eset_ql_combat, toadjust.hyb, doplot = T, adjust4subgroup = F, mycovarspalte = "subgroup") {
    # toadjust.hyb='bigbatch' mycovarspalte = 'subgroup'
    
    
    esetname = deparse(substitute(eset))
    message("Doing calcAnova() on ", esetname)
    time_anova <- Sys.time()
    
    if (adjust4subgroup == F) {
        anovazeile = parse(text = paste("modelsum_after <- esApply(eset, 1, function(x) anova(lm(x~", paste(toadjust.hyb, collapse = " + "), ")))"))
        print(anovazeile)
        eval(anovazeile)
        time_anova = Sys.time() - time_anova
        print(time_anova)  #
        
        
        
        # modelparameter extrahieren
        
        all_pval = sapply(modelsum_after, function(x) {
            rn = row.names(x)
            y = x[, "Pr(>F)"]
            names(y) = rn
            return(y)
        })
        all_pval = t(all_pval)
        all_pval = all_pval[, toadjust.hyb]
        print(head(all_pval, 5))
        
        if (doplot == T) {
            
            
            ## qqplots
            ggd.qqplot(all_pval, main = toadjust.hyb)
            mtext(paste(esetname, "\nRot: x=y"), cex = 0.6)
            
        }
        return(all_pval)
    }
    
    if (adjust4subgroup == T) {
        
        message("... adjusting ANOVA for subgroup ...")
        toadjust.hyb_mod = paste(c(mycovarspalte, toadjust.hyb), collapse = " + ")
        anovazeile = parse(text = paste("modelsum_after <- esApply(eset, 1, function(x) anova(lm(x~", paste(toadjust.hyb_mod, collapse = " + "), ")))"))
        # anovazeile = parse(text=paste('modelsum_after <- esApply(eset[1:100,], 1, function(x) anova(lm(x~', paste(toadjust.hyb_mod, collapse=' + '), ')))'))
        print(anovazeile)
        eval(anovazeile)
        time_anova = Sys.time() - time_anova
        print(time_anova)  #
        
        
        
        # modelparameter extrahieren
        
        all_pval = sapply(modelsum_after, function(x) {
            rn = row.names(x)
            y = x[, "Pr(>F)"]
            names(y) = rn
            return(y)
        })
        all_pval = t(all_pval)
        all_pval = all_pval[, toadjust.hyb]
        print(head(all_pval, 5))
        
        if (doplot == T) {
            
            
            ## qqplots
            ggd.qqplot(all_pval, main = toadjust.hyb)
            mtext(paste(esetname, "\nRot: x=y"), cex = 0.6)
            
        }
        return(all_pval)
    }
    
    stopifnot(adjust4subgroup %in% c(T, F))
    
}
