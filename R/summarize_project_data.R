library(dplyr)
library(readr)
library(stringr)
library(readxl)
library(tidyr)
library(lubridate)
library(reshape2)

# Read in the pipeline log (spreadsheet, which I've been manually updating)
pipelineLogFile <- "~/Dropbox/BRI/rnaseqPipeline/PipelineLog.xlsx"
pipelineLog <- read_excel(pipelineLogFile)

# Remove reprocessing projects and projects with incomplete info
pipelineStandard <- pipelineLog %>% 
    filter(!str_detect(aligned_subfolder, "MultiFlowcell"),
           !is.na(results_folder))

# For each lib in the project, use a shell command to pull out the
# corresponding number of reads from the metrics file
count_reads <- Vectorize(function(metricsFile, libId) {
    getReadsCmd <- sprintf("grep %s %s | cut -d ',' -f 52", 
                           libId, metricsFile)
    numReads <- as.numeric(system(getReadsCmd, intern = TRUE)) * 1e-6
    return(numReads)
})

# For a given unaligned library folder, sum up the file size of all gzipped
# FASTQs and convert to GB
calc_file_size <- Vectorize(function(unalignedPath, libId) {
    libFolder <- list.dirs(unalignedPath) %>%
        .[str_detect(., libId)]
    lsCmd <- sprintf("ls -l %s | grep .gz | awk '{sum+=$5} END {print sum}'" , 
                     libFolder)
    fileSize <- as.numeric(system(lsCmd, intern = TRUE)) * 1e-9
    return(fileSize)
})

# Get list of projects
projects <- pipelineStandard$project

projectList <- list(proj = pipelineStandard$project, 
                    fc = pipelineStandard$flowcell)

# Do the work...
for (i in seq(1, length(projectList$proj))) {
    project_i <- projectList$proj[i]
    flowcell_i <- projectList$fc[i]
    print(sprintf("%s: %s", i, project_i))
    
    projectInfo <- pipelineStandard %>% 
        filter(project == project_i & flowcell == flowcell_i) %>% 
        mutate(path = str_c(aligned_data, aligned_subfolder, 
                            results_folder, sep = "/"),
               unaligned_path = str_c(unaligned_data, unaligned_subfolder,
                                      project_folder, sep = "/"))
    
    projectPath <- projectInfo$path %>% str_replace("mnt", "Volumes")
    projectLibs <- list.files(file.path(projectPath, "counts" )) %>% 
        str_extract("lib[0-9]+") %>% 
        na.omit()
    

    metricsFile <- list.files(file.path(projectPath, "metrics")) %>% 
        .[str_detect(., "combined")] %>% 
        file.path(projectPath, "metrics", .)
    
    unalignedPath <- projectInfo$unaligned_path %>% 
        str_replace("mnt", "Volumes")
    
    if (project_i == projects[1]) {
        pipelineThroughput <- data_frame(date = projectInfo$processed,
                                         flowcell = projectInfo$flowcell,
                                         project = project_i,
                                         prep = projectInfo$lib_prep,
                                         lib = projectLibs) %>% 
            mutate(preAlignReads = count_reads(metricsFile, lib),
                   fastqSize = calc_file_size(unalignedPath, lib))
    } else {
        ptTmp <- data_frame(date = projectInfo$processed,
                            flowcell = projectInfo$flowcell,
                            project = project_i,
                            prep = projectInfo$lib_prep,
                            lib = projectLibs) %>% 
            mutate(preAlignReads = count_reads(metricsFile, lib),
                   fastqSize = calc_file_size(unalignedPath, lib))
        pipelineThroughput <- pipelineThroughput %>% 
            bind_rows(ptTmp)
    }
}

# Summarize / create report
pipelineThroughput <- pipelineThroughput %>% 
    mutate(jobMonth = as.character(month(date, label = TRUE)))

truSeqSummary <- pipelineThroughput %>% 
    select(-date) %>% 
    filter(str_detect(prep, "TruSeq"),
           !is.na(preAlignReads)) %>% 
    melt(id.vars = c("preAlignReads", "fastqSize", "lib"), 
         measure.vars = c("jobMonth", "flowcell", "project"), 
         variable.name = "binType") %>% 
    group_by(binType, value) %>% 
    summarise(numReads = sum(preAlignReads),
              gbData = sum(fastqSize),
              numLibs = n()) %>% 
    group_by(binType) %>% 
    summarise_each(funs(min, max, median, sum),
                   numReads, gbData, numLibs) %>% 
    select(one_of(sort(names(.))))

truSeqLibs <- pipelineThroughput %>% 
    filter(str_detect(prep, "TruSeq")) %>% 
    mutate(numReads = preAlignReads,
           gbData = fastqSize,
           binType = "library") %>% 
    group_by(binType) %>% 
    summarise_each(funs(min, max, median, sum),
                   numReads, gbData) %>% 
    select(one_of(sort(names(.))))

truSeqSummary <- truSeqSummary %>% 
    full_join(truSeqLibs)

nexteraSummary <- pipelineThroughput %>% 
    select(-date) %>% 
    filter(str_detect(prep, "Nextera"),
           !is.na(preAlignReads)) %>% 
    melt(id.vars = c("preAlignReads", "fastqSize", "lib"), 
         measure.vars = c("jobMonth", "flowcell", "project"), 
         variable.name = "binType") %>% 
    group_by(binType, value) %>% 
    summarise(numReads = sum(preAlignReads),
              gbData = sum(fastqSize),
              numLibs = n()) %>% 
    group_by(binType) %>% 
    summarise_each(funs(min, max, median, sum),
                   numReads, gbData, numLibs) %>% 
    select(one_of(sort(names(.))))

nexteraLibs <- pipelineThroughput %>% 
    filter(str_detect(prep, "Nextera")) %>% 
    mutate(numReads = preAlignReads,
           gbData = fastqSize,
           binType = "library") %>% 
    group_by(binType) %>% 
    summarise_each(funs(min, max, median, sum),
                   numReads, gbData) %>% 
    select(one_of(sort(names(.))))

nexteraSummary <- nexteraSummary %>% 
    full_join(nexteraLibs)
