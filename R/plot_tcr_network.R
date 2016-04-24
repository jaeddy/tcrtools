library(reshape2)

imgt_mixcr$tcrs %>% melt(measure.vars = c("trav_jxn", "trbv_jxn"))

tcr <- list(name = "TCR", children = list(
    list(name = "TRAV", children = list()), 
    list(name = "TRBV", children = list())
    ))

# test --------------------------------------------------------------------


get_children <- function(network, node) {

    tree_list = list(name = node)
    
    children <- network %>% 
        filter(target == node) %>% 
        .$source %>% 
        unique()
    
    if (length(children)) {
        tree_list$children = list()
        for (i in seq(1, length(children))) {
            tree_list$children[[i]] = get_children(network, children[i])
        }
    }
    
    return(tree_list)
}

jxns <- c(imgt_mixcr$tcrs$trav_jxn, imgt_mixcr$tcrs$trbv_jxn) %>% 
    unique()
tcr <- list(name = "tcr",
            children = lapply(as.list(1:length(jxns)), function(x) { 
                return(
                    get_children(tcr_sankey, jxns[x]))}
            ))
tcr$children[1]
