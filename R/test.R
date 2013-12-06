# labkey.url.base <- "https://www.immunespace.org/"
# # labkey.url.path <- "/Yale/SDY162/"
# labkey.url.path <- "/Emory/SDY61/"
# e_m_a <- c("EXPM0001", "EXPM0002")
# # e_m_a <- c("EXPM0001")
# 
# eset <- read_exprs_mat(e_m_a)
# 
# # f2 <- function(lub = get('lub', .GlobalEnv)){
# lub <- "https://www.immunespace.org/"
# lub <- "https://posey.fhcrc.org/"
# 
# 
# f2 <- function(lub){
#   if(missing(lub)){
#     lub <- get('lub', .GlobalEnv)
#   }
#   return(lub)
# }
# 
# 
# 
# 
#     err <- tryCatch(expr={lub <- get('lub', .GlobalEnv)},
#                     error=function(ERR){
#                       message("'labkey.url.base' should either be present in the global environment or be passed as an argument")
#                       return(NA)
#                     }
#                     )
#     if(is.na(err)){
#       stop("")
#     }
#     
#   }
#   return(lub)
# }


# #test localhost
# l.u.b <- "http://localhost:8080/labkey"
# l.u.p <- "/home/WTFLK/"
# res <- labkey.selectRows(baseUrl = l.u.b, folderPath = l.u.p, schemaName = "lists",
#                   queryName = "gene_expression_reagents", colNameOpt = "rname")
# 
# #test importRows
# ti <- data.frame(GEA = paste("GEA", 1:5, sep="_"), description = "testImport")
# microbenchmark(import=
#                  Rlabkey:::labkey.importRows(baseUrl = l.u.b, folderPath = l.u.p,
#                   schemaName = "lists", queryName = "test_labkeyimportRows",
#                   toImport = ti),
#                insert=
#                  Rlabkey:::labkey.insertRows(baseUrl = l.u.b, folderPath = l.u.p,
#                   schemaName = "lists", queryName = "test_labkeyimportRows",
#                   toInsert = ti),
#               times=3)
# #Unit: milliseconds
# #expr        min         lq     median         uq        max neval
# #import   187.1982   209.5302   231.8622   243.7528   255.6434     3
# #insert 32720.7433 33143.1791 33565.6149 33785.6573 34005.6997     3
# 
# 
# nr <- 500
# ti <- data.frame(GEA = paste("GEA", 1:nr, sep="_"), description = "testImport")
# td <- data.frame(GEA = paste("GEA", 1:nr, sep="_"), stringsAsFactors=FALSE)
# res <- labkey.importRows(baseUrl = l.u.b, folderPath = l.u.p, schemaName = "lists",
#                          queryName = "test_deleteRow", toImport = ti)
# res <- labkey.deleteRows(baseUrl=l.u.b, folderPath=l.u.p,
#                     schemaName="lists", queryName="test_deleteRow",toDelete  = td)

# Delete 500 rows
#user  system elapsed 
#0.108   0.000  34.734 