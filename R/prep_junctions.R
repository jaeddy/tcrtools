library(stringr)
library(dplyr)

# Universal function to filter clonotypes, with options to specify minimum
# junction length (# amino acids) or clone count (for MiXCR), and whether to
# remove non-productive junctions or junctions without conserved leading C 
# residue or trailing F/W residue
filter_jxns <- function(jxn_df, min_length = 4, min_count = 1,
                        productive = TRUE, conserved = TRUE) {
    jxn_df <- jxn_df %>% 
        filter(str_detect(v_gene, "(?<=TR)[A-B]"), # TRA/B only
               str_detect(j_gene, "(?<=TR)[A-B]"), # TRA/B only
               str_extract(v_gene, "(?<=TR)[A-B]") == 
                   str_extract(j_gene, "(?<=TR)[A-B]")) # segments must match
               
    if(min_length > 0) {
        jxn_df <- jxn_df %>% 
            filter(str_length(junction) >= min_length)
    }

    
    if(productive | conserved) {
        jxn_df <- jxn_df %>% 
            filter(str_detect(junction, "^C"), # jxn AA must start with C
                   str_detect(junction, "(W|F)$")) # junction AA must end with W/F
    }
    
    if(productive) {
        jxn_df <- jxn_df %>% 
            filter(str_detect(junction, "^((?!(\\*|_|#)).)*$")) # productive only
    }
    
    if("cln_count" %in% names(jxn_df)) {
        jxn_df <- jxn_df %>% 
            filter(cln_count >= min_count)
    }
    
    return(jxn_df %>% 
               distinct(sample, v_gene, j_gene, junction, .keep_all = TRUE))
}


# There should be only one alpha and one beta junction for each lib; select
# the top hit for each lib, sorting first by clone count and second by
# alignment _score
select_top_jxns <- function(mixcr_jxns, allow_multi = FALSE, detail = FALSE) {

   mixcr_jxns <- mixcr_jxns %>% 
        mutate(a_b = ifelse(str_detect(v_gene, "TRA"), "TRA", "TRB")) %>% 
        group_by(lib_id, a_b) %>% 
        arrange(desc(cln_count), desc(v_gene_score)) 
   
   if (allow_multi) {
       mixcr_jxns <- mixcr_jxns %>% 
           slice(1:2) %>% 
           ungroup() %>% 
           group_by(lib_id, a_b) %>% 
           mutate(count_diff = c(0, diff(cln_count)) / lag(cln_count),
                  score_diff = c(0, diff(v_gene_score)) / lag(v_gene_score),
                  best = ifelse(is.na(count_diff), TRUE, FALSE),
                  second = ifelse(a_b == "TRA" & !best &
                                      (count_diff > -0.5 | score_diff > -0.1),
                                  TRUE, FALSE)) %>% 
           ungroup() %>% 
           filter(best | second)
   } else {
       mixcr_jxns <- mixcr_jxns %>% 
           slice(1) %>% 
           ungroup()
   }

   if (detail) {
       mixcr_jxns %>% 
           select(-a_b, -count_diff, -score_diff, -best, -second)
   } else {
       mixcr_jxns %>% 
           select(lib_id, v_gene, j_gene, junction)
   }

}

### DEPRECATED ###

# Function to filter IMGT junctions

# filter_imgt_jxns <- function(imgt_jxns, min_length = 6) {
#     
#     imgt_jxns <- imgt_jxns %>% 
#         filter(str_detect(v_gene, "^((?![C-G]).)*$"),
#                str_detect(j_gene, "^((?![C-G]).)*$"),
#                str_detect(junction, "^C"),
#                str_detect(junction, "^((?!(\\*|_)).)*$"),
#                !duplicated(.[, c(1:4)]),
#                str_length(junction) > min_length)
#     
#     return(imgt_jxns)
# }

# Function to filter MiXCR junctions

# filter_mixcr_jxns <- function(mixcr_jxns, min_count = 1, min_length = 4) {
#     
#     mixcr_jxns <- mixcr_jxns %>% 
#         filter(str_detect(v_gene, "TR[A-B]"), # TRA/B only
#                str_detect(j_gene, "TR[A-B]"), # TRA/B only
#                str_extract(v_gene, "(?<=TR)[A-B]") == 
#                    str_extract(j_gene, "(?<=TR)[A-B]"), # segments must match
#                str_detect(junction, "^C"), # jxn AA must start with C
#                str_detect(junction, "^((?!(\\*|_)).)*$"), # functional only
#                #                !duplicated(.[, c(1, 3:5)]), # remove dups
#                cln_count >= min_count,
#                str_length(junction) >= min_length)
#     
#     return(mixcr_jxns)
# }
