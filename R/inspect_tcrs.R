
# Function to combine TRAV and TRBV genes into TCRs
construct_tcrs <- function(jxns_df, any = FALSE) {

    trav_df <- jxns_df %>% 
        filter(str_detect(v_gene, "TRAV")) %>% 
        select(-j_gene) %>% 
        dplyr::rename(trav_gene = v_gene, 
               trav_jxn = junction)
    
    trbv_df <- jxns_df %>% 
        filter(str_detect(v_gene, "TRBV")) %>% 
        select(-j_gene) %>% 
        dplyr::rename(trbv_gene = v_gene,
               trbv_jxn = junction)

    if (!any) {
        tcr_df <- inner_join(trav_df, trbv_df, by = "lib_id")
    } else {
        tcr_df <- full_join(trav_df, trbv_df, by = "lib_id") %>% 
            mutate(trav_gene = ifelse(is.na(trav_gene), "no_trav_gene", trav_gene),
                   trbv_gene = ifelse(is.na(trbv_gene), "no_trbv_gene", trbv_gene),
                   trav_jxn = ifelse(is.na(trav_jxn), "no_trav_jxn", trav_jxn),
                   trbv_jxn = ifelse(is.na(trbv_jxn), "no_trbv_jxn", trbv_jxn))
    }
    return(tcr_df)
}

build_sankey_network <- function(tcrs_all, chain = c("A", "B", "both")) {
    sources <- tcrs_all$tcr_source
    tcrs_all <- tcrs_all %>% 
        mutate(value = ifelse(tcr_source == sources[1], 
                              as.numeric(str_extract(lib_id, "[0-9]+")) + 0.1, 
                              as.numeric(str_extract(lib_id, "[0-9]+")) + 0.2))
    
    if (chain == "A") {
        tcr_sankey <- tcrs_all %>% 
            select(source = lib_id, target = trav_gene, value) %>% 
            distinct() %>% 
            bind_rows(tcrs_all %>% 
                          select(source = trav_gene, target = trav_jxn, value) %>% 
                          distinct())
    } else if (chain == "B") {
        tcr_sankey <- tcrs_all %>% 
            select(source = lib_id, target = trbv_gene, value) %>% 
            distinct() %>% 
            bind_rows(tcrs_all %>% 
                          select(source = trbv_gene, target = trbv_jxn, value) %>% 
                          distinct())
    } else if (chain == "both") {
        tcr_sankey <- tcrs_all %>% 
            select(source = lib_id, target = trav_gene, value) %>% 
            distinct() %>% 
            bind_rows(tcrs_all %>% 
                          select(source = lib_id, target = trbv_gene, value) %>% 
                          distinct()) %>% 
            bind_rows(tcrs_all %>% 
                          select(source = trav_gene, target = trav_jxn, value) %>% 
                          distinct()) %>% 
            bind_rows(tcrs_all %>% 
                          select(source = trbv_gene, target = trbv_jxn, value) %>% 
                          distinct()) %>% 
            distinct()
    } 
}

build_sankey_plot <- function(tcr_sankey, sankey_height = 600) {
    sankey_plot <- rCharts$new()
    sankey_plot$setLib('/Users/jaeddy/code/github/resources/rCharts_d3_sankey/libraries/widgets/d3_sankey')
    sankey_plot$setTemplate(script = "/Users/jaeddy/code/github/resources/rCharts_d3_sankey/libraries/widgets/d3_sankey/layouts/tcrChart.html")
    
    sankey_plot$set(
        data = tcr_sankey,
        nodeWidth = 15,
        nodePadding = 10,
        layout = 32,
        width = 800,
        height = sankey_height 
    )
    return(sankey_plot)
}
