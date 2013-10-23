1_ Create the eset
  - Argument:
    - expression_matrix_accession (Potentially more than one)
  - Read study.gene_expression_matrices with filtering on e_m_a
    - This will include all infos (PTL, TRT, ...)
  - Read the files:
    - expression_matrix
    - fdata
  - Create the expressionset

2_ Limma
  - Arguments:
    - contrast = "study_time_collected"
    - FDR thresh = 0.1
  - Value:
    - topTable or list of topTable

3_ Writing
  - Arguments:
    - topTables
    - analysis_accession
    - analysis_description
  - Fill the tables:
    - lists.gene_expression_analysis
    - lists.gene_expression_analysis_result_wide
