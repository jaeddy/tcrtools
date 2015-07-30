library(dplyr)
library(stringr)
library(tidyr)
library(reshape2)
library(ggplot2)

### mixcr analysis

rm(list = ls()) # clear workspace

# Analysis functions ------------------------------------------------------

# Function to combine individual * mixcrClns.txt files using unix: 1) add header
# from one file to newfile. 2) grep contents of all files and 3) append them to
# newfile
combine_mixcr_clns <- function(mixcrDir) {
    # Select one file to grab the column headers
    mixcrTmpFile <- list.files(mixcrDir, full.names = TRUE) %>% 
        str_extract(".*Clns.txt") %>% 
        na.omit() %>% 
        .[1]
    
    # Build and use Unix commands
    mixcrCombinedFile <- file.path(mixcrDir, "compiled_mixcr_output.txt")
    
    headerCmd <- sprintf("head -1 %s > %s", mixcrTmpFile, mixcrCombinedFile)
    system(headerCmd)
    
    compileCmd <- sprintf("grep '*' %s/*Clns.txt >> %s",
                          mixcrDir, mixcrCombinedFile)
    system(compileCmd)
    
    return(mixcrCombinedFile)
}

# Function to read in MiXCR data then extract & format variables
format_mixcr_clns <- function(mixcrCombinedFile) {
    # Read in MiXCR clones
    mixcrClns <- read.delim(mixcrCombinedFile)
    
    # Extract & format key variables
    mixcrClns <- mixcrClns %>% 
        transmute(libID = as.character(str_match(Clone.count, "lib[0-9]+")),
                  clnCount = as.numeric(str_match(Clone.count, "(?<=:)[0-9]+")),
                  vGene = str_extract(All.V.hits,
                                      "TR[A-Z]+[0-9]*(\\-[0-9])*(DV[0-9]+)*"),
                  jGene = str_extract(All.J.hits,
                                      "TR[A-Z]+[0-9]*(\\-[0-9][A-Z]*)*"),
                  junction = as.character(AA..seq..CDR3))
    
    return(mixcrClns)
}

# Function to filter MiXCR clones
filter_mixcr_clns <- function(mixcrClns, minCount = 0, minLength = 6) {
    
    mixcrClns <- mixcrClns %>% 
        filter(str_detect(vGene, "^((?![C-G]).)*$"),
               str_detect(junction, "^C"),
               str_detect(junction, "^((?!(\\*|_)).)*$"),
               !duplicated(.[, c(1, 3:5)]),
               clnCount > minCount,
               str_length(junction) > minLength)
    
    return(mixcrClns)
}

# Function to plot distribution of MiXCR clone counts & junction lengths
plot_mixcr_summary <- function(mixcrClns) {
    mixcrSummary <- mixcrClns %>% 
        mutate(log2_clnCount = log2(clnCount),
               junctionLength = str_length(junction)) %>% 
        melt(measure.vars = c("log2_clnCount", "junctionLength")) %>% 
        ggplot(aes(x = value)) +
        geom_histogram() +
        facet_wrap(~ variable)
    print(mixcrSummary)
}

# Function to filter IMGT clones
filter_imgt_clns <- function(imgtClns, minLength = 6) {
    
    imgtClns <- imgtClns %>% 
        filter(str_detect(V.GENE, "^((?![C-G]).)*$"),
               str_detect(J.GENE, "^((?![C-G]).)*$"),
               str_detect(JUNCTION, "^C"),
               str_detect(JUNCTION, "^((?!(\\*|_)).)*$"),
               !duplicated(.[, c(1:4)]),
               str_length(JUNCTION) > 6)
    
    return(imgtClns)
}

# Compare TCR overlap at different clone count cutoffs
compare_tcrs <- function(imgtClns, mixcrClns, maxCutoff = 40) {
    imgtTcrAll <- imgtClns %>% 
        unite_("cloneStr", c("libID", "JUNCTION"), sep = "::")
    
    clnCountRange <- seq(0, maxCutoff, 1)
    tcrBoth <- rep(0, maxCutoff + 1)
    tcrImgt <- rep(0, maxCutoff + 1)
    tcrMixcr <- rep(0, maxCutoff + 1)
    
    for (cutoff in clnCountRange) {
        mixcrTcrAll <- mixcrClns %>% 
            filter(clnCount > cutoff) %>% 
            select(-clnCount) %>% 
            unite_("cloneStr", c("libID", "junction"), sep = "::")
        
        tcrBoth[cutoff+1] <- intersect(imgtTcrAll$cloneStr, 
                                       mixcrTcrAll$cloneStr) %>% length()
        tcrImgt[cutoff+1] <- setdiff(imgtTcrAll$cloneStr, 
                                     mixcrTcrAll$cloneStr) %>% length()
        tcrMixcr[cutoff+1] <- setdiff(mixcrTcrAll$cloneStr, 
                                      imgtTcrAll$cloneStr) %>% length()
    }
    
    compareClns <- data_frame(cloneCount = clnCountRange,
                              imgtAll = nrow(imgtTcrAll),
                              maxIntersect = max(tcrBoth),
                              both = tcrBoth,
                              imgtOnly = tcrImgt, mixcrOnly = tcrMixcr)
    
    return(compareClns)
}

# Plot overlap values over range of clone count cutoffs
plot_tcr_mixcr <- function(mixcrClns, maxCutoff = 40) {
    clnCountRange <- seq(0, maxCutoff, 1)
    tcrMixcr <- rep(0, maxCutoff + 1)
    
    for (cutoff in clnCountRange) {
        mixcrTcrNum <- mixcrClns %>% 
            filter(clnCount > cutoff) %>% 
            nrow()
        
        tcrMixcr[cutoff+1] <- mixcrTcrNum
    }
    
    mixcrClnCounts <- data_frame(cloneCount = clnCountRange,
                                 clonotypes = tcrMixcr)
    
    cutoffVal <- mixcrClnCounts$clonotypes[which(mixcrClnCounts$cloneCount == 25)]
    mixcrClnCounts %>% 
        ggplot(aes(x = cloneCount, y = clonotypes)) +
        geom_line(size = 1.5) +
        geom_vline(x = 25) +
        geom_hline(y = cutoffVal)
}

# Plot overlap values over range of clone count cutoffs
plot_tcr_comp <- function(compareClns) {
    compareClnsLong <- melt(compareClns, id.vars = "cloneCount", 
                            variable.name = "tcrSource", 
                            value.name = "clonotypes")
    
    compareClnsLong %>% 
        ggplot(aes(x = cloneCount, y = clonotypes)) +
        geom_line(aes(colour = tcrSource), size = 1.5)
}

# Wrapper function for individual steps above
run_comparison <- function(mixcrDir, imgtFile) {
    imgtClns <- read.delim(imgtFile) %>% 
        filter_imgt_clns() %>% 
        filter(!str_detect(libID, "lib2"))
    
    mixcrCombinedFile <- combine_mixcr_clns(mixcrDir)
    mixcrClns <- format_mixcr_clns(mixcrCombinedFile) %>% 
        filter_mixcr_clns()
    
    compareClns <- compare_tcrs(imgtClns, mixcrClns)
    
    print(plot_tcr_comp(compareClns))
    return(compareClns)
}





