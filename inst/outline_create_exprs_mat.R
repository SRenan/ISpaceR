1_ Reading
  - Read the table
  - Detect whether it is Illumina or Affy
  - Read all files
    ? What if some files have different featuredata (SDY161)
  - Make ExpressionSet

2_ Analysis
  - Split the eset
  - Subset (e.g: PBMC & macro should not be normalized together)
  - Normalize using rma or lumi (optional)
  - Subset (e.g: PolyIC L is not needed)

3_ Writing
  - Write the matrices to files
  - Write the fData to files
  - Fill the matrix tables
    - lists.expression_matrix: accession, description
    - lists.expression_matrix_2_biosample: accession, file, biosample
  - Fill the reagent tables
    - lists.features_2_genes
    - lists.gene_expression_platforms: reagent_accession, feature_mapping_file
    - lists.gene_expression_reagent_2_matrix: reagent_accession, matrix_accession

