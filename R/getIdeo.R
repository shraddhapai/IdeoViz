getIdeo <- function(ideoSource) {
    cytoTable <- NULL
    
    session <- browserSession();
    genome(session) <- ideoSource;
    cytoTable <- getTable(ucscTableQuery(session,"cytoBandIdeo"))
    
    return(cytoTable)
### (data.frame) ideogram table
}#, ex=function() {
#getIdeo("mm9")
#}
