### mixcr analysis

# Function to combine individual * mixcr_clns.txt files using unix: 1) add header
# from one file to newfile. 2) grep contents of all files and 3) append them to
# newfile
combine_mixcr_outputs <- function(mixcr_dir, output_dir, project) {
    # Select one file to grab the column headers
    mixcr_tmp_file <- data_frame(
        file = list.files(mixcr_dir, full.names = TRUE)) %>% 
        mutate(size = file.size(file)) %>% 
        filter(str_detect(file, "_clns.txt"),
               size > 0) %>% # make sure not to select empty files
        slice(1) %>% 
        select(file) %>% 
        unlist()
    
    # Build and use Unix commands
    mixcr_combined_file <- file.path(output_dir, 
                                     paste(project, "compiled_mixcr_output.txt", 
                                           sep = "_"))
    if (!file.exists(mixcr_combined_file)) {
        header_cmd <- sprintf("head -1 %s > %s", mixcr_tmp_file, mixcr_combined_file)
        system(header_cmd)
        
        compile_cmd <- sprintf("grep '*' %s/*_clns.txt >> %s",
                              mixcr_dir, mixcr_combined_file)
        system(compile_cmd)
    }
    
    return(mixcr_combined_file)
}

# Function to compile and format IMGT results using R
compile_imgt_output <- function(summary_file, output_dir, project) {
    # Build and use Unix commands
    imgt_file <- file.path(output_dir, 
                          paste(project, "compiled_imgt_output.txt", 
                                sep = "_"))
    
    if (!file.exists(imgt_file)) {
        imgt_summary <- read.delim(summary_file, stringsAsFactors = FALSE)
        
        imgt_compiled <- imgt_summary %>% 
            filter(Functionality == "productive") %>% 
            select(Sequence.ID, V.GENE.and.allele, J.GENE.and.allele, 
                   AA.JUNCTION) %>% 
            transmute(lib_id = str_extract(Sequence.ID, "lib[0-9]+"),
                      v_gene = str_extract(V.GENE.and.allele, 
                                           "TR.*?(?=(\\*))"),
                      j_gene = str_extract(J.GENE.and.allele, 
                                           "TR.*?(?=(\\*))"),
                      junction = AA.JUNCTION)
        
        write.table(imgt_compiled, imgt_file, 
                    sep = "\t", quote = FALSE, row.names = FALSE)
    }
    
    return(imgt_file)
}

# Function to read in MiXCR data then extract & format variables
format_mixcr_jxns <- function(mixcr_combined_file) {
    # Read in MiXCR clones
    mixcr_jxns <- read.delim(mixcr_combined_file)
    
    # Extract & format key variables
    mixcr_jxns <- mixcr_jxns %>% 
        transmute(lib_id = as.character(str_match(Clone.count, "lib[0-9]+")),
                  cln_count = as.numeric(str_match(Clone.count, "(?<=:)[0-9]+")),
                  v_gene = str_extract(All.V.hits,
                                      "TR[A-Z]+[0-9]*(\\-[0-9])*(DV[0-9]+)*"),
                  v_gene_score = str_extract(All.V.hits, "(?<=\\()[0-9]+"),
                  j_gene = str_extract(All.J.hits, 
                                      "TR[A-Z]+[0-9]*(\\-[0-9][A-Z]*)*"),
                  j_gene_score = str_extract(All.J.hits, "(?<=\\()[0-9]+"),
                  junction = as.character(AA..seq..CDR3))
    
    return(mixcr_jxns)
}

# Function to filter MiXCR junctions
filter_mixcr_jxns <- function(mixcr_jxns, min_count = 0, min_length = 6) {
    
    mixcr_jxns <- mixcr_jxns %>% 
        filter(str_detect(v_gene, "^((?![C-G]).)*$"),
               str_detect(j_gene, "^((?![C-G]).)*$"),
               str_detect(junction, "^C"),
               str_detect(junction, "^((?!(\\*|_)).)*$"),
               !duplicated(.[, c(1, 3:5)]),
               cln_count > min_count,
               str_length(junction) > min_length)
    
    return(mixcr_jxns)
}

# Function to filter IMGT junctions
filter_imgt_jxns <- function(imgt_jxns, min_length = 6) {
    
    imgt_jxns <- imgt_jxns %>% 
        filter(str_detect(v_gene, "^((?![C-G]).)*$"),
               str_detect(j_gene, "^((?![C-G]).)*$"),
               str_detect(junction, "^C"),
               str_detect(junction, "^((?!(\\*|_)).)*$"),
               !duplicated(.[, c(1:4)]),
               str_length(junction > min_length))
    
    return(imgt_jxns)
}

select_top_jxns <- function(mixcr_jxns) {
    # There should be only one alpha and one beta junction for each lib; select
    # the top hit for each lib, sorting first by clone count and second by
    # alignment _score
    mixcr_jxns %>% 
        filter(lib_id == "lib8472") %>% 
        mutate(a_b = ifelse(str_detect(v_gene, "TRA"), "TRA", "TRB")) %>% 
        group_by(lib_id, a_b) %>% 
        arrange(cln_count, v_gene_score) %>% 
        slice(1) %>% 
        ungroup() %>% 
        select(-cln_count, -v_gene_score, -j_gene_score, -a_b)
}

