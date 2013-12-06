




```r
library(ggplot2)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(Rlabkey)
library(Biobase)
library(pheatmap)
library(glmnet)
library(ISpaceR)
```


```r
mypalette <- rev(brewer.pal(name="PiYG", n=11))
brewer.qual <- "Paired"
IS_theme <- theme_bw(base_size=14)
```




```r
lub <- labkey.url.base <- "https://posey.fhcrc.org/"
lup <- labkey.url.path <- "/Emory/SDY61/"

param <-list(day = 7, GEA_acc = "GEA1", EM_acc = "EXPM0001", FDR_thresh = 0.02, FC_thresh = 0)
```



```r

GEM <- read_exprs_mat(param$EM_acc)[[1]]

GER_filter <- makeFilter(c("analysis_accession", "EQUAL", param$GEA_acc),
                         c("adj_p_val", "LESS_THAN", param$FDR_thresh))
ds_GER <- labkey.selectRows(baseUrl = labkey.url.base, folderPath = labkey.url.path, schemaName = "lists", queryName = "gene_expression_analysis_result_wide", colNameOpt="rname", colFilter=GER_filter)
GER <- data.table(ds_GER)

HAI_filter <- makeFilter(c("subject_accession", "IN", paste(GEM$subject_accession, collapse=";")))
ds_HAI <- labkey.selectRows(baseUrl=labkey.url.base, schemaName="study", folderPath=labkey.url.path, queryName ="hai", colNameOpt="rname", colFilter=HAI_filter)
HAI <- data.table(ds_HAI)
```



```r
#Filter-subset
#Filter out subjects with missing timepoints
#Filter selected probes
pd <- data.table(pData(GEM))
timepoints <- c(0, param$day)
pd <- pd[,keep:=(sum(study_time_reported%in%timepoints)==length(timepoints)), by="biosample_accession_name,subject_accession"]
pd <- pd[keep==TRUE]
pd <- pd[study_time_reported %in% timepoints]
probes_selected <- unique(GER$feature_id)
GES <- GEM[probes_selected, pd$biosample_accession] #sorted by subject, time
HAI <- HAI[subject_accession %in% GES$subject_accession & biosample_accession_name %in% GES$biosample_accession_name]
```



```r
#lFC
lFC <- exprs(GES)[, GES$study_time_reported == param$day] - exprs(GES)[, GES$study_time_reported == 0]
```



```r
#hai
#Class subjects as high/low responders
HAI <- HAI[,list(study_time_collected=study_time_collected, hai_response=value_reported/value_reported[study_time_collected==0]),by="virus_strain,biosample_accession_name,subject_accession"]
HAI <- HAI[study_time_collected==28]
HAI <- HAI[,list(hai_response=max(hai_response)),by="subject_accession"]
HAI <- HAI[,responder:=ifelse(hai_response<=2, "Low", ifelse(hai_response>=4, "High", "Medium"))]
```



```r
#glmnet
x <- t(lFC)
y <- log2(HAI$hai_response)
fit_hai <- glmnet(x, y, alpha=0.5)
cv_fit_hai <- cv.glmnet(x, y)
coef_hai <- predict(fit_hai, s=cv_fit_hai$lambda.min, type="coefficients")
```



```r
#Selected probes
selected_probes <- names(which(coef_hai[,1] > 0))
selected_probes <- fData(GES)[selected_probes[grep("Intercept", selected_probes, invert=TRUE)],]
#If the user selected less variables than observation, no selection is done
if(length(probes_selected) < nrow(HAI)){
  x.selected <- x
  selected_probes <- fData(GES)[probes_selected,]
} else{
  x.selected <- x[, colnames(x) %in% selected_probes$feature_id, drop=FALSE]
}

if(nrow(selected_probes) < 2){
  opts_chunk$set(eval=FALSE, cache=FALSE)
  stop("No probes were selected as predictive. You may try to change the filtering criteria.")
}
```



```r
#lasso
relasso <- lm(y~x.selected)
sum_relasso <- summary(relasso)
sum_relasso_coef <- sum_relasso$coefficients
predictor_table <- cbind(selected_probes, sum_relasso_coef[,c("t value","Pr(>|t|)")][-1,])
setnames(predictor_table, c("t value","Pr(>|t|)"), c("t-value","p-value"))
```



```r
#Display
#Table
predictor_table[,"gene_symbol"] <- paste0('<a href="http://www.genenames.org/cgi-bin/quick_search.pl?.cgifields=type&type=equal&num=50&search=', gsub(";","+",predictor_table[,"gene_symbol"]),'&submit=Submit" target="_blank">', predictor_table[,"gene_symbol"], '</a>')
```



```r
kable(predictor_table, digits=2, format="html", row.names=FALSE, table.attr="id=\"res_table\"")
```

<table id="res_table">
 <thead>
  <tr>
   <th> gene_symbol </th>
   <th> feature_id </th>
   <th> t-value </th>
   <th> p-value </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td> <a href="http://www.genenames.org/cgi-bin/quick_search.pl?.cgifields=type&type=equal&num=50&search=MCCC2&submit=Submit" target="_blank">MCCC2</a> </td>
   <td> 209623_at </td>
   <td>  0.39 </td>
   <td> 0.70 </td>
  </tr>
  <tr>
   <td> <a href="http://www.genenames.org/cgi-bin/quick_search.pl?.cgifields=type&type=equal&num=50&search=NA&submit=Submit" target="_blank">NA</a> </td>
   <td> 244447_at </td>
   <td>  1.69 </td>
   <td> 0.11 </td>
  </tr>
  <tr>
   <td> <a href="http://www.genenames.org/cgi-bin/quick_search.pl?.cgifields=type&type=equal&num=50&search=CLPTM1L&submit=Submit" target="_blank">CLPTM1L</a> </td>
   <td> 223020_at </td>
   <td> -0.02 </td>
   <td> 0.98 </td>
  </tr>
  <tr>
   <td> <a href="http://www.genenames.org/cgi-bin/quick_search.pl?.cgifields=type&type=equal&num=50&search=SPCS2+LOC653566&submit=Submit" target="_blank">SPCS2;LOC653566</a> </td>
   <td> 201240_s_at </td>
   <td>  0.42 </td>
   <td> 0.68 </td>
  </tr>
  <tr>
   <td> <a href="http://www.genenames.org/cgi-bin/quick_search.pl?.cgifields=type&type=equal&num=50&search=IGKC+IGK@&submit=Submit" target="_blank">IGKC;IGK@</a> </td>
   <td> 216576_x_at </td>
   <td>  0.74 </td>
   <td> 0.47 </td>
  </tr>
  <tr>
   <td> <a href="http://www.genenames.org/cgi-bin/quick_search.pl?.cgifields=type&type=equal&num=50&search=NA&submit=Submit" target="_blank">NA</a> </td>
   <td> 217179_x_at </td>
   <td>  0.51 </td>
   <td> 0.62 </td>
  </tr>
  <tr>
   <td> <a href="http://www.genenames.org/cgi-bin/quick_search.pl?.cgifields=type&type=equal&num=50&search=PRDX1&submit=Submit" target="_blank">PRDX1</a> </td>
   <td> 208680_at </td>
   <td>  1.30 </td>
   <td> 0.21 </td>
  </tr>
</tbody>
</table>

<script type="text/javascript" charset="utf-8">
    $(document).ready(function() {
        $('#res_table').dataTable();
    } );
</script>


```r
options(markdown.HTML.header = unlist(sapply(system.file("misc", c("vignette.css", "datatables.txt"), package = "knitr"), readLines)))
```


```r
#Heatmap
annotation <- data.frame(hai=log2(HAI$hai_response))
rownames(annotation) <- colnames(lFC)
mat <- lFC[ as.character(predictor_table$feature_id), order(annotation$hai)]
pheatmap(mat, dendrogram="none", cluster_cols=FALSE, cluster_rows=TRUE, 
         breaks = seq(-1, 1, length.out = length(mypalette) + 1), 
         show_rownames=TRUE, show_colnames=FALSE, scale="none",cluster_method="ward",
         cluster_distance="correlation",color=mypalette, annotation=annotation, 
         annotation_colors=list(hai=grey(10:0/10)))
```

<img src="figure/heatmap.png" title="plot of chunk heatmap" alt="plot of chunk heatmap" style="display: block; margin: auto;" />



```r
#Prediction on training set
x.test <- data.frame(t(lFC[ as.character(selected_probes$feature_id),]), check.names=FALSE)
hai_pred <- as.matrix(x.test)%*%as.vector((coef(relasso)[-1]))+coef(relasso)[1]
df <- data.frame(observed=HAI$hai_response, fitted=relasso$fitted.values)
ggplot(df, aes(x=log2(observed), y=fitted)) + geom_point() + 
  geom_smooth(method="lm") +IS_theme
```

<img src="figure/prediction.png" title="plot of chunk prediction" alt="plot of chunk prediction" style="display: block; margin: auto;" />

