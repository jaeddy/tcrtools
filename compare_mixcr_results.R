### mixcr analysis

rm(list = ls()) # clear workspace

# Analysis functions ------------------------------------------------------

# Function to combine individual * mixcrClns.txt files using unix: 1) add header
# from one file to newfile. 2) grep contents of all files and 3) append them to
# newfile
combine_mixcr_clns <- function(mixcrDir, outputDir, project) {
    # Select one file to grab the column headers
    mixcrTmpFile <- data_frame(file = list.files(mixcrDir, full.names = TRUE)) %>% 
        mutate(size = file.size(file)) %>% 
        filter(str_detect(file, "Clns.txt"),
               size > 0) %>% 
        slice(1) %>% 
        select(file) %>% 
        unlist()
    
    # Build and use Unix commands
    mixcrCombinedFile <- file.path(outputDir, 
                                   paste(project, "compiled_mixcr_output.txt", 
                                         sep = "_"))
    if (!file.exists(mixcrCombinedFile)) {
        headerCmd <- sprintf("head -1 %s > %s", mixcrTmpFile, mixcrCombinedFile)
        system(headerCmd)
        
        compileCmd <- sprintf("grep '*' %s/*Clns.txt >> %s",
                              mixcrDir, mixcrCombinedFile)
        system(compileCmd)
    }

    return(mixcrCombinedFile)
}

# Function to compile and format IMGT results using R
compile_imgt_clns <- function(summaryFile, outputDir, project) {
#     # Select summary file
#     summaryFile <- list.files(imgtDir, full.names = TRUE) %>% 
#         str_extract(".*Summary.*") %>% 
#         na.omit()
    
    # Build and use Unix commands
    imgtFile <- file.path(outputDir, 
                          paste(project, "compiled_imgt_output.txt", 
                                sep = "_"))
    
    if (!file.exists(imgtFile)) {
        imgtSummary <- read.delim(summaryFile, stringsAsFactors = FALSE)
        
        imgtCompiled <- imgtSummary %>% 
            filter(Functionality == "productive") %>% 
            select(Sequence.ID, V.GENE.and.allele, J.GENE.and.allele, 
                   AA.JUNCTION) %>% 
            transmute(libID = str_extract(Sequence.ID, "lib[0-9]+"),
                      V.GENE = str_extract(V.GENE.and.allele, 
                                           "TR.*?(?=(\\*))"),
                      J.GENE = str_extract(J.GENE.and.allele, 
                                           "TR.*?(?=(\\*))"),
                      JUNCTION = AA.JUNCTION)
        
        write.table(imgtCompiled, imgtFile, 
                    sep = "\t", quote = FALSE, row.names = FALSE)
    }
    
    return(imgtFile)
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
               str_detect(jGene, "^((?![C-G]).)*$"),
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
    # mixcrAll <- rep(0, maxCutoff + 1)
    tcrBoth <- rep(0, maxCutoff + 1)
    tcrImgt <- rep(0, maxCutoff + 1)
    tcrMixcr <- rep(0, maxCutoff + 1)
    
    for (cutoff in clnCountRange) {
        mixcrTcrAll <- mixcrClns %>% 
            filter(clnCount >= cutoff) %>% 
            select(-clnCount) %>% 
            unite_("cloneStr", c("libID", "junction"), sep = "::")
        
        # mixcrAll[cutoff+1] <- nrow(mixcrTcrAll)
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
plot_tcr_mixcr <- function(mixcrClns, maxCutoff = 40, cutoffLine = 25) {
    clnCountRange <- seq(0, maxCutoff, 1)
    tcrMixcr <- rep(0, maxCutoff + 1)
    
    for (cutoff in clnCountRange) {
        mixcrTcrNum <- mixcrClns %>% 
            filter(clnCount >= cutoff) %>% 
            nrow()
        
        tcrMixcr[cutoff+1] <- mixcrTcrNum
    }
    
    mixcrClnCounts <- data_frame(cloneCount = clnCountRange,
                                 clonotypes = tcrMixcr)
    
    cutoffVal <- mixcrClns %>%
        filter(clnCount >= cutoffLine) %>% 
        nrow()
    
    mixcrClnCounts %>% 
        ggplot(aes(x = cloneCount, y = clonotypes)) +
        geom_line(size = 1.5) +
        geom_vline(x = cutoffLine) +
        geom_hline(y = cutoffVal)
}

# Plot overlap values over range of clone count cutoffs
plot_tcr_comp <- function(compareClns, cutoffLine = 25) {
    compareClnsLong <- melt(compareClns, id.vars = "cloneCount", 
                            variable.name = "tcrSource", 
                            value.name = "clonotypes")
    
    compareClnsLong %>% 
        ggplot(aes(x = cloneCount, y = clonotypes)) +
        geom_line(aes(colour = tcrSource), size = 1.5) +
        geom_vline(x = cutoffLine)
}

# Function to combine TRAV and TRBV genes into TCRs
construct_tcrs <- function(clnsDat, cutoff = 10, method = c("mixcr", "imgt")) {
    if (method == "mixcr") {
        travDat <- clnsDat %>% 
            filter(clnCount > cutoff,
                   str_detect(vGene, "TRAV")) %>% 
            select(-jGene) %>% 
            rename(clnCountA = clnCount, TRAV = vGene, 
                   TRAVjunction = junction)
        
        trbvDat <- clnsDat %>% 
            filter(clnCount > cutoff,
                   str_detect(vGene, "TRBV")) %>% 
            select(-jGene) %>% 
            rename(clnCountB = clnCount, TRBV = vGene,
                   TRBVjunction = junction)
    } else {
        travDat <- clnsDat %>% 
            filter(str_detect(V.GENE, "TRAV")) %>% 
            select(-J.GENE) %>% 
            rename(TRAV = V.GENE, TRAVjunction = JUNCTION)
        
        trbvDat <- clnsDat %>% 
            filter(str_detect(V.GENE, "TRBV")) %>% 
            select(-J.GENE) %>% 
            rename(TRBV = V.GENE, TRBVjunction = JUNCTION)
    }
    
    tcrDat <- inner_join(travDat, trbvDat, by = c("libID" = "libID"))
    return(tcrDat)
}


# sankey ntwork -----------------------------------------------------------


build_sankey_network <- function(tcrsAll, chain = c("A", "B", "both")) {
    tcrsAll <- tcrsAll %>% 
        mutate(tcr = str_c(TRAV, TRBV, sep = ":"),
               value = ifelse(tcrSource == "IMGT", 0.49, 0.51))
    
    if (chain == "A") {
        tcrSankey <- tcrsAll %>% 
            select(source = libID, target = TRAV, value) %>% 
            bind_rows(tcrsAll %>% 
                          select(source = TRAV, target = TRAVjunction, value))
    } else if (chain == "B") {
        tcrSankey <- tcrsAll %>% 
            select(source = libID, target = TRBV, value) %>% 
            bind_rows(tcrsAll %>% 
                          select(source = TRBV, target = TRBVjunction, value))
    } else if (chain == "both") {
        tcrSankey <- tcrsAll %>% 
            select(source = libID, target = TRAV, value) %>% 
            bind_rows(tcrsAll %>% 
                          select(source = libID, target = TRBV, value)) %>% 
            bind_rows(tcrsAll %>% 
                          select(source = TRAV, target = TRAVjunction, value)) %>% 
            bind_rows(tcrsAll %>% 
                          select(source = TRBV, target = TRBVjunction, value))
    } else {
        tcrSankey <- tcrsAll %>% 
            mutate(libID = str_c(libID, "[TRAV]", sep = " ")) %>% 
            select(source = libID, target = TRAV, value) %>% 
            bind_rows(tcrsAll %>% 
                          select(source = TRAV, target = TRBV, value)) %>% 
            bind_rows(tcrsAll %>% 
                          mutate(libID = str_c(libID, "[TRBV]", sep = " ")) %>% 
                          select(source = TRBV, target = libID, value)) 
    }
    
}



# sankey fxn --------------------------------------------------------------


sankey_plot <- function(tcrSankey, sankeyHeight = 600) {
    sankeyPlot <- rCharts$new()
    sankeyPlot$setLib('~/code/github/resources/rCharts_d3_sankey/libraries/widgets/d3_sankey')
    sankeyPlot$setTemplate(script = "~/code/github/resources/rCharts_d3_sankey/libraries/widgets/d3_sankey/layouts/tcrChart.html")

    sankeyPlot$set(
        data = tcrSankey,
        nodeWidth = 15,
        nodePadding = 10,
        layout = 32,
        width = 800,
        height = sankeyHeight 
    )
    return(sankeyPlot)
}




