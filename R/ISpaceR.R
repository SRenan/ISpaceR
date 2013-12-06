# labkey.url.base <- "https://www.immunespace.org/"
# labkey.url.path <- "/Yale/SDY162/"
# labkey.url.path <- "/Emory/SDY61/"
e_m_a <- c("EXPM0001", "EXPM0002")
# e_m_a <- c("EXPM0001")
##
#' Read expression matrix
#' 
#' Reads the gene_expression_matrix table and creates an ExpressionSet for the 
#' selected matrix or matrices.
#' 
#' @param e_m_a A \code{character}. The accession number of the matrix or matrices
#' to read.
#' 
#' @return A \code{list} of \code{ExpressionSet}, with one element for each 
#' given matrix.
#' 
#' @note This function uses \code{RCurl} to query ImmuneSpace.
#' Therefore, a valid .netrc file is needed when executed from a remote location.
#' 
#' @seealso \code{\link{limma_standard}}, \code{\link{write_gea}}
#' 
#' @export
#' @docType methods
#' @rdname read_exprs_mat
#' @aliases read_exprs_mat
#' 
#' @examples
#' labkey.url.base <- "https://www.immunespace.org/"
#' labkey.url.path <- "/Emory/SDY61/"
#' read_exprs_mat("EXPM0001")
#' 
#' @importFrom data.table fread
#' @importFrom RCurl getURL
#' @importFrom Rlabkey makeFilter labkey.selectRows
#' @importClassesFrom Biobase AnnotatedDataFrame ExpressionSet
##
read_exprs_mat <- function(e_m_a){
  # Get the main table
  em_filter <- makeFilter(c("expression_matrix_accession", "IN", paste(e_m_a, collapse=";")))
  exMat <- labkey.selectRows(baseUrl = labkey.url.base, folderPath = labkey.url.path, schemaName = "study",
                             colFilter=em_filter, queryName = "gene_expression_matrices", colNameOpt = "rname")
  rownames(exMat) <- exMat$biosample_accession
  #
  umat <- unique(exMat[, c("expression_matrix_accession", "download_link_download_link", "feature_mapping_file", "feature_mapping_file_link_download_link", "matrix_description_description")])
  links <- lapply(umat$download_link_download_link, URLdecode)
  urls <- lapply(links, getURL, .opts=list(netrc=TRUE))
  dts <- lapply(paste0("\t", urls), fread)
  names(dts) <- umat$expression_matrix_accession
  exprs_dfs <- lapply(dts, function(X){ data.frame(X[,2:ncol(X), with=FALSE], row.names=X[,1, with=FALSE][[1]])})

  fdata_links <- lapply(umat[[grep("feature.*download", colnames(umat))]], URLdecode)
  fdata_urls <- lapply(fdata_links, getURL, .opts=list(netrc=TRUE))
  fdata_dfs <- lapply(fdata_urls, function(X){ read.table(textConnection(X), header=TRUE, sep="\t")})
  names(fdata_dfs) <- umat$expression_matrix_accession

  lEset <- vector('list', length(e_m_a))
  names(lEset) <- e_m_a
  for(exm in 1:length(e_m_a)){
    exprs <- as.matrix(exprs_dfs[[exm]])
    phenodata <- exMat[ exMat$expression_matrix_accession==e_m_a[exm],]
    phenodata <- phenodata[ match(colnames(exprs), phenodata$biosample_accession),]
    phenodata <- as(phenodata, "AnnotatedDataFrame")
    fdata <- as(fdata_dfs[[exm]], "AnnotatedDataFrame")
    descr <- paste(umat[exm, "expression_matrix_accession"], ": ", umat[exm, "matrix_description_description"])
    eset <- ExpressionSet(assayData = exprs, phenoData = phenodata, featureData = fdata)
    eset@experimentData@title <- descr
    lEset[[exm]] <- eset  
  }
  return(lEset)
}
  
##
#' Standard limma analysis
#' 
#' A wrapper around limma functions to do a standard gene expression analysis.
#' Fit a linear model for each gene, computes some statistics and returns the result in a table (or list of tables).
#' By default the function tries a time course experiment.
#' 
#' @param eset An \code{ExpressionSet}. The expression set to analyse.
#' @param contrast A \code{character}. The characteristic used to make contrasts. It should be a valid column name in the phenoData of eset.
#' @param FDRthresh A \code{numeric}. The threshold used for filtering the False Discovery Rate.
#' 
#' @return A named \code{list} of \code{data.frame} of the same length as the number of contrasts.
#' 
#' @seealso \code{\link{read_exprs_mat}} \code{\link{write_gea}} \code{\link{lmFit}} \code{\link{topTable}}
#' 
#' @export
#' @docType methods
#' @rdname limma_standard
#' @aliases limma_standard
#' 
#' @importFrom limma lmFit eBayes makeContrasts topTable
##
limma_standard <- function(eset, contrast = "study_time_reported", FDRthresh = 0.1){
  f <- factor(eset[[contrast]])
  design <- model.matrix(~0+f+eset$subject_accession)
  fit <- lmFit(eset, design)
  contrasts <- makeContrasts(paste0("f", f[2:length(levels(f))], "-f", f[[1]]), levels=design)
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)
  tt <- list()
  for(cont in 1:ncol(contrasts)){
    tmp_tt <- topTable(fit2, coef=cont, number="Inf", p.value=FDRthresh)
    tmp_tt <- tmp_tt[, c("feature_id", "gene_symbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
    setnames(tmp_tt, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val"),
             c("log_fc", "ave_expr", "statistic", "p_value", "adj_p_val"))
    tt <- c(tt, list(tmp_tt))
  }
  names(tt) <- colnames(contrasts)
  return(tt)
}

##
#' Write gene expression analysis result to the database
#'  
#' This function writes the output of a standard limma analysis to the database and fill the appropriate tables.
#' 
#' @param topTable A \code{data.frame}, as returned by the \code{limma_standard} function,
#'  or a \code{list} of \code{data.frame}.
#' @param analysis_accession A \code{character} vector. The accession associated with the gene expression analysis.
#'  The vector should be of the same length as topTable if topTable is a list.
#' @param description A \code{character} vector. The description of the analysis. 
#' 
#' @return A \code{boolean}. TRUE if eveything went as expected.
#'
#' @seealso \code{\link{read_exprs_mat}} \code{\link{limma_standard}} \code{\link{labkey.importRows}}
#' 
#' @export
#' @docType methods
#' @rdname write_gea
#' @aliases write_gea
#' 
#' @importFrom Rlabkey labkey.importRows
##
write_gea <- function(topTable, analysis_accession, description){
  # topTable: a topTable or list of topTables
  # analysis_accession: a string
  if(is.list(topTable)){
    if(any(sapply(list(topTable, analysis_accession, description), length)!=length(topTable))){
      stop("The three arguments should be of the same length. Given: ",
           paste(sapply(list(topTable, analysis_accession, description), length), collapse=", "))
    }
  }
  if(any(!grepl("^GEA", analysis_accession))){
    stop("Invalid analysis_accession: ", paste(analysis_accession[which(!grepl("^GEA", analysis_accession))], collapse=", "))
  }  

  for(ind in topTable){
    #gene_expression_analysis_results
    toImport <- topTable[[ind]]
    toImport$analysis_accession <- analysis_accession[ind]

    setnames(toImpot, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val"),
             c("log_fc", "ave_expr", "statistic", "p_value", "adj_p_val"))
    labkey.importRows(baseUrl=labkey.url.base, folderPath=labkey.url.path, schemaName="lists",
                      queryName="gene_expression_analysis_result_wide", toImport=toImport)
    toImport <- 
    toImport$description <- description[ind]    
    labkey.importRows(baseUrl=labkey.url.base, folderPath=labkey.url.path, schemaName="lists",
                      queryName="gene_expression_analysis", toImport=toImport)
    
  }
  
}

# ##########################################################################
# ##########################################################################
# ##########################################################################
# f <- function(eset, subset){
#   if(missing(subset)){
#   } else{
#     eset <- subset(eset, subset=subset)
#   }
#   return(eset)
# }
#   
# 
# #Takes GE_files 
# # - Isn't it always the same matrix (0 args + opts)
# # - Do we put everything in GE_files?
# # - Do we do the querying in R? Prolly not
# # - Do we detect data type
# # - Do we normalize?
# # - Does this populate the tables with the new matrices as well? (e.g: 'list.expression_matrices')
# # - Always one matrix per cohort? It's not always convenient (SDY162)
# #return the expression matrice(s) and the feature mapping file
# get_meta <- function(normalize=TRUE, subset){
#   #query
#   ds_GE_files <- labkey.selectRows(baseUrl = labkey.url.base, folderPath = labkey.url.path, schemaName = "study",
#                                    queryName = "gene_expression_files", colNameOpt = "rname")
#   dt_GE_files <- data.table(ds_GE_files)
#   
#   #subset
#   if(missing(subset)){
#     call <- substitute(TRUE)
#   } else{
#     call <- substitute(subset)
#   }
#   rows <- eval(expr=call, envir=dt_GE_files, enclos=parent.frame())
#   dt_GE_files <- dt_GE_files[rows,]
#   
#   #split by cohort
#   for(arm in dt_GE_files$arm_name){
#     #The reading should be different when not executed from ISpace
#     affyBatch <- ReadAffy(filenames=paste(labkey.file.root, "rawdata",
#                           "gene_expression", pd_TIV_07$file_info_name, sep = "/"))
#     if(normalize){
#       eset <- rma(affyBatch)
#     }
#   }
# }
# 