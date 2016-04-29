library(readr)
library(stringr)
library(dplyr)
library(parallel)

# clean column names of data frame
clean_headers <- function(df) {
    headers <- names(df)
    headers <- str_to_lower(headers)
    headers <- str_replace_all(headers, "( )+", "_")
    headers <- str_replace_all(headers, "[^[:alnum:]_]", "_")
    headers <- str_replace_all(headers, "_+", "_")
    headers <- str_replace_all(headers, "_$", "")
    
    names(df) <- headers
    return(df)
}

# read IMGT results Summary from file
read_imgt_summary <- function(file) {
    imgt_df <- read_tsv(file) %>% 
        clean_headers() %>% 
        filter(functionality != "No results")
    
    return(imgt_df)
}

# parse raw IMGT clonotype results from Summary file
parse_imgt_summary <- function(imgt_df) {
    imgt_df %>% 
        select(sequence_id, v_gene_and_allele, v_region_score, v_region_identity_nt,
               j_gene_and_allele, j_region_score, j_region_identity_nt,
               aa_junction, sequence) %>% 
        mutate(v_gene = str_extract(v_gene_and_allele, 
                                       "TR.*?(?=(\\*))"),
               v_region_score = as.integer(v_region_score),
               j_gene = str_extract(j_gene_and_allele, 
                                    "TR.*?(?=(\\*))"),
               j_region_score = as.integer(j_region_score),
               junction = aa_junction) %>% 
        select(one_of("sequence_id", 
                      "v_gene", "v_region_score", "v_region_identity_nt",
                      "j_gene", "j_region_score", "j_region_identity_nt",
                      "junction", "sequence"))
    
}

# read and parse IMGT results from list of archive (.txz) files
read_imgt <- function(file_list = NULL, folder = NULL, 
                      sample_regex = "(lib|SRR)[0-9]+") {
    if(is.null(file_list) & is.null(folder)) {
        stop("Input must be provided for either `file_list` or `folder` argument.")
    }
    
    if(!is.null(folder)) {
        file_list <- list.files(folder, full.names = TRUE) %>% 
            .[str_detect(tolower(.), ".txz")]
    }
    
    if(!is.list(file_list)) {
        file_list <- as.list(file_list)
    }
    
    jxn_df <- mclapply(file_list, function(x) {
        extract_cmd <- sprintf("tar -xf '%s' -O '1_Summary.txt'", x)
        pipe(extract_cmd) %>% 
            read_file() %>% 
            read_imgt_summary() %>% 
            parse_imgt_summary() %>% 
            mutate(sequence_id = str_extract(sequence_id, sample_regex)) %>% 
            rename(sample = sequence_id)
    }) %>% 
        bind_rows()
    
    return(jxn_df)
}

# read MiXCR results from file
read_mixcr_clones <- function(file) {
    mixcr_df <- read_delim(file, delim = "\t") %>% 
        clean_headers()
}

# parse raw MiXCR clonotype results
parse_mixcr_clones <- function(mixcr_df) {
    mixcr_df %>% 
        transmute(cln_count = clone_count,
                  v_gene = str_extract(all_v_hits,
                                       "TR[A-Z]+[0-9]*(\\-[0-9])*(DV[0-9]+)*"),
                  v_gene_score = str_extract(all_v_hits, "(?<=\\()[0-9]+") %>% 
                      as.integer(),
                  v_align = all_v_alignment,
                  j_gene = str_extract(all_j_hits, 
                                       "TR[A-Z]+[0-9]*(\\-[0-9][A-Z]*)*"),
                  j_align = all_j_alignment,
                  j_gene_score = str_extract(all_j_hits, "(?<=\\()[0-9]+") %>% 
                      as.integer(),
                  junction = as.character(aa_seq_cdr3)) %>% 
        rowwise() %>% 
        mutate(v_region_identity_nt = str_split(v_align, pattern = "\\|") %>% 
                   unlist() %>% 
                   .[length(.) - 2],
               v_region_score = str_split(v_align, pattern = "\\|") %>% 
                   unlist() %>% 
                   .[length(.)] %>% 
                   as.integer(),
               j_region_identity_nt = str_split(j_align, pattern = "\\|") %>% 
                   unlist() %>% 
                   .[length(.) - 2],
               j_region_score = str_split(j_align, pattern = "\\|") %>% 
                   unlist() %>% 
                   .[length(.)] %>% 
                   as.integer()) %>% 
        select(one_of(c( "cln_count", "v_gene", "v_gene_score",
                        "v_region_identity_nt", "v_region_score",
                        "j_gene", "j_gene_score",
                        "j_region_identity_nt",
                        "j_region_score",
                        "junction")))
}

# read and parse MiXCR results from list of files
read_mixcr <- function(file_list = NULL, folder = NULL, 
                       sample_regex = "(lib|SRR)[0-9]+") {
    if(is.null(file_list) & is.null(folder)) {
        stop("Input must be provided for either `file_list` or `folder` argument.")
    }
    
    if(!is.null(folder)) {
        file_list <- list.files(folder, full.names = TRUE) %>% 
            .[str_detect(tolower(.), "clns.txt")]
    }
    
    if(!is.list(file_list)) {
        file_list <- as.list(file_list)
    }
    
    jxn_df <- mclapply(file_list, function(x) {
        if(file.size(x)) {
            read_mixcr_clones(x) %>% 
                parse_mixcr_clones() %>% 
                mutate(sample = str_extract(basename(x), sample_regex))
        }
    }) %>% 
        bind_rows() %>% 
        select(one_of("sample", setdiff(names(.), "sample")))
    
    return(jxn_df)
}