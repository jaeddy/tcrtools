

# Function to plot distribution of MiXCR clone counts & junction lengths
plot_mixcr_summary <- function(mixcr_jxns) {
    mixcr_summary <- mixcr_jxns %>% 
        mutate(log2_cln_count = log2(cln_count),
               jxn_length = str_length(junction)) %>% 
        melt(measure.vars = c("log2_cln_count", "jxn_length")) %>% 
        ggplot(aes(x = value)) +
        geom_histogram() +
        facet_wrap(~ variable, scales = "free_x") +
        theme_gray() +
        theme(axis.title.x = element_blank()) +
        ylab("num_junctions")
    return(mixcr_summary)
}

# Plot overlap values over range of clone count cutoffs
plot_mixcr_jxn_dist <- function(mixcr_jxns, max_cutoff = 40, cutoff_line = 25) {
    cln_count_range <- seq(0, max_cutoff, 1)
    mixcr_jxn_dist <- rep(0, max_cutoff + 1)
    
    for (cutoff in cln_count_range) {
        mixcr_jxn_num <- mixcr_jxns %>% 
            filter(cln_count >= cutoff) %>% 
            nrow()
        
        mixcr_jxn_dist[cutoff+1] <- mixcr_jxn_num
    }
    
    mixcr_cln_counts <- data_frame(clone_count = cln_count_range,
                                   num_junctions = mixcr_jxn_dist)
    
    cutoff_val <- mixcr_jxns %>%
        filter(cln_count >= cutoff_line) %>% 
        nrow()
    
    mixcr_cln_counts %>% 
        ggplot(aes(x = clone_count, y = num_junctions)) +
        geom_line(size = 1.5) +
        geom_vline(x = cutoff_line) +
        geom_hline(y = cutoff_val) +
        theme_gray()
}

# Compare junction overlap at different clone count cutoffs
compare_jxn_dists <- function(imgt_jxns, mixcr_jxns, max_cutoff = 40) {
    imgt_jxns_all <- imgt_jxns %>% 
        unite_("jxn_str", c("lib_id", "junction"), sep = "::")
    
    cln_count_range <- seq(0, max_cutoff, 1)
    # mixcr_all <- rep(0, max_cutoff + 1)
    jxns_both <- rep(0, max_cutoff + 1)
    jxns_imgt <- rep(0, max_cutoff + 1)
    jxns_mixcr <- rep(0, max_cutoff + 1)
    
    for (cutoff in cln_count_range) {
        mixcr_jxns_all <- mixcr_jxns %>% 
            filter(cln_count >= cutoff) %>% 
            select(-cln_count) %>% 
            unite_("jxn_str", c("lib_id", "junction"), sep = "::")
        
        # mixcr_all[cutoff+1] <- nrow(mixcr_tcr_all)
        jxns_both[cutoff+1] <- intersect(imgt_jxns_all$jxn_str, 
                                         mixcr_jxns_all$jxn_str) %>% length()
        jxns_imgt[cutoff+1] <- setdiff(imgt_jxns_all$jxn_str, 
                                       mixcr_jxns_all$jxn_str) %>% length()
        jxns_mixcr[cutoff+1] <- setdiff(mixcr_jxns_all$jxn_str, 
                                        imgt_jxns_all$jxn_str) %>% length()
    }
    
    compare_jxns <- data_frame(clone_count = cln_count_range,
                               imgt_all = nrow(imgt_jxns_all),
                               max_intersect = max(jxns_both),
                               both = jxns_both,
                               imgt_only = jxns_imgt, mixcr_only = jxns_mixcr)
    
    return(compare_jxns)
}


# Plot overlap values over range of clone count cutoffs
plot_jxn_comp <- function(compare_jxns, cutoff_line = 20) {
    compare_jxns_melt <- melt(compare_jxns, id.vars = "clone_count", 
                              variable.name = "jxn_source", 
                              value.name = "num_junctions")
    
    compare_jxns_melt %>% 
        ggplot(aes(x = clone_count, y = num_junctions)) +
        geom_line(aes(colour = jxn_source), size = 1.5) +
        geom_vline(x = cutoff_line) +
        theme_gray()
}