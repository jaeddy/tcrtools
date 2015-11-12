
# Function to combine TRAV and TRBV genes into TCRs
construct_tcrs <- function(jxns_df) {

    trav_df <- jxns_df %>% 
        filter(str_detect(v_gene, "TRAV")) %>% 
        select(-j_gene) %>% 
        rename(trav_gene = v_gene, 
               trav_jxn = junction)
    
    trbv_df <- jxns_df %>% 
        filter(str_detect(v_gene, "TRBV")) %>% 
        select(-j_gene) %>% 
        rename(trbv_gene = v_gene,
               trbv_jxn = junction)

    tcr_df <- inner_join(trav_df, trbv_df, by = c("lib_id" = "lib_id"))
    return(tcr_df)
}

build_sankey_network <- function(tcrs_all, chain = c("A", "B", "both")) {
    tcrs_all <- tcrs_all %>% 
        mutate(tcr = str_c(trav_gene, trbv_gene, sep = ":"),
               value = ifelse(tcr_source == "IMGT", 0.49, 0.51))
    
    if (chain == "A") {
        tcr_sankey <- tcrs_all %>% 
            select(source = lib_id, target = trav_gene, value) %>% 
            bind_rows(tcrs_all %>% 
                          select(source = trav_gene, target = trav_jxn, value))
    } else if (chain == "B") {
        tcr_sankey <- tcrs_all %>% 
            select(source = lib_id, target = trbv_gene, value) %>% 
            bind_rows(tcrs_all %>% 
                          select(source = trbv_gene, target = trbv_jxn, value))
    } else if (chain == "both") {
        tcr_sankey <- tcrs_all %>% 
            select(source = lib_id, target = trav_gene, value) %>% 
            bind_rows(tcrs_all %>% 
                          select(source = lib_id, target = trbv_gene, value)) %>% 
            bind_rows(tcrs_all %>% 
                          select(source = trav_gene, target = trav_jxn, value)) %>% 
            bind_rows(tcrs_all %>% 
                          select(source = trbv_gene, target = trbv_jxn, value))
    } else {
        tcr_sankey <- tcrs_all %>% 
            mutate(lib_id = str_c(lib_id, "[TRAV]", sep = " ")) %>% 
            select(source = lib_id, target = trav_gene, value) %>% 
            bind_rows(tcrs_all %>% 
                          select(source = trav_gene, target = trbv_gene, value)) %>% 
            bind_rows(tcrs_all %>% 
                          mutate(lib_id = str_c(lib_id, "[TRBV]", sep = " ")) %>% 
                          select(source = trbv_gene, target = lib_id, value)) 
    }
    
}

build_sankey_plot <- function(tcr_sankey, sankey_height = 600) {
    sankey_plot <- rCharts$new()
    sankey_plot$setLib('~/code/github/resources/rCharts_d3_sankey/libraries/widgets/d3_sankey')
    sankey_plot$setTemplate(script = "~/code/github/resources/rCharts_d3_sankey/libraries/widgets/d3_sankey/layouts/tcrChart.html")
    
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
