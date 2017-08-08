get_badExpressed_Inds <- function(sample_overview_l5_subset, numgen_cutoff, subgroupname) {
    genesfound <- sample_overview_l5_subset$Detected.Genes..0.01.
    summary(genesfound)
    
    low_cutoff <- median(genesfound) - numgen_cutoff * IQR(genesfound)
    low_cutoff
    high_cutoff <- median(genesfound) + numgen_cutoff * IQR(genesfound)
    high_cutoff
    
    bad_genenum <- sample_overview_l5_subset[sample_overview_l5_subset$Detected.Genes..0.01. < low_cutoff | sample_overview_l5_subset$Detected.Genes..0.01. > high_cutoff, 
        "new_ID"]
    
    
    
    # gefunden
    par(mfrow = c(2, 1))
    boxplot(genesfound, main = "Detected Genes according to Illumina (Parameter 'Detected Genes  0.01')", cex.main = 0.9, ylim = c(min(c(low_cutoff, genesfound)), max(c(high_cutoff, 
        genesfound))))
    abline(h = c(low_cutoff, high_cutoff, median(genesfound)), col = c("red", "red", "blue"))
    mtext(paste(length(bad_genenum), " individuals outside Median + - ", numgen_cutoff, "*IQR (", round(low_cutoff), " - ", round(high_cutoff), ")"))
    hist(genesfound, xlab = "Detected Genes according to illumina p level 0.01", breaks = min(floor(dim(sample_overview_l5_subset)[1]/1.3), 200), main = unique(subgroupname), 
        cex.main = 0.9, xlim = c(min(c(low_cutoff, genesfound)), max(c(high_cutoff, genesfound))))
    abline(v = c(low_cutoff, high_cutoff, median(genesfound)), col = c("red", "red", "blue"))
    
    res = c()
    res$bad_genenum = bad_genenum
    res$high_cutoff = high_cutoff
    res$low_cutoff = low_cutoff
    res$median = median(genesfound)
    res$genesfound = genesfound
    res
}
