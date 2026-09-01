%% ParTI on BacDrop Ciprofloxacin Data
% Leveraging struxcture provided by ParTI toolbox (Uri Alon Lab)

% Load the numeric gene expression matrix
% Rows = cells (samples), Columns = genes
% Load expression matrix (CSV with header and row labels)
%geneExpression = readmatrix("path/to/parti_expression_matrix.csv", 'NumHeaderLines', 1);
% Load expression matrix as a table (includes barcodes and gene names)
T = readtable("path/to/ParTI_expression_matrix.csv");

% Get the size of the table
[rows, cols] = size(T);

% Display the shape
fprintf('The shape of the table is %d rows × %d columns.\n', rows, cols);

% Drop the first column (barcodes) and convert to numeric matrix
geneExpression = table2array(T(:, 2:end));

% Load gene names (as a list of strings)
geneNames = importdata("path/to/ParTI_gene_names.txt", ',');

% Load discrete labels 
[discrAttrNames, discrAttr] = read_enriched_csv("path/to/ParTI_discrete_attributes.csv", ',');

% there are no continuous attributes in this simplified version
contAttr = [];
contAttrNames = {};
GOcat2Genes = [];

% Clean up variable names (optional)
discrAttrNames = regexprep(discrAttrNames, '_', ' ');

%% Run ParTI with randomization controls for single p-value

% Full ParTI version with randomization controls
[arc, arcOrig, pc, errs, pval] = ParTI(geneExpression, ...
    1, 8, ...
    discrAttrNames, discrAttr, 0, ...
    contAttrNames, contAttr, ...
    GOcat2Genes, ...
    0.05, 'BacDropCipro_enrichmentAnalysis_0208_SISAL_balanced');

% Display the overall p-value
fprintf('Overall statistical significance p-value: %.6f\n', pval);
%% 

save('ParTI_output_0208_SISAL_balanced.mat', 'arc', 'arcOrig', 'pc', 'errs', 'pval')


