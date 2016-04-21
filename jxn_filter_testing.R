# experimenting with junction filtering
 
srp060416_jxns_filtered <- mixcr$jxns

srp060416_jxns_all <- format_mixcr_jxns(mixcr$file)


srp060416_jxns_ab <- srp060416_jxns_all %>% 
    filter(str_detect(v_gene, "TR[A-B]V"))


srp060416_jxns_ab <- srp060416_jxns_ab %>% 
    filter(str_length(junction) >= 4)

srp060416_jxns_ab_match <- srp060416_jxns_ab %>% 
    filter(str_detect(v_gene, "TRA") & str_detect(j_gene, "TRA") |
               str_detect(v_gene, "TRB") & str_detect(j_gene, "TRB"))
                                                      

srp060416_jxns_ab_cys_phe_trp <- srp060416_jxns_ab_match %>% 
    filter(str_detect(junction, "^C"),
           str_detect(junction, "(F|W)$"))

srp060416_jxns_prod <- srp060416_jxns_ab_cys_phe_trp %>% 
    filter(str_detect(junction, "^((?!(\\*|_)).)*$"))

srp060416_jxns_top <- srp060416_jxns_prod %>% 
    select_top_jxns(allow_multi = TRUE)

prod_jxn_tally <- srp060416_jxns_prod %>% 
    group_by(lib_id) %>% 
    summarise(n_trav = sum(str_detect(v_gene, "TRA")),
              n_trbv = sum(str_detect(v_gene, "TRB")))
