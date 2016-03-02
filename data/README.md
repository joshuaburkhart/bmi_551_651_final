# File Descriptions

### expression.txt
tab-delimited text file with expression values for 18,632 genes for each of the 39 cell lines.

### subtypes.txt
tab-delimited text file of subtypes (basal, luminal, claudin-low and normal-like) for each of 39 cell lines

### training_set_answers.txt
tab-delimited text file of the correct classification of 0 (non-responsive) or 1 (responsive) for each combination of 25 cell lines and 12 drugs.

### scoring_and_test_set_id_mappings.csv
comma-delimited text file of the id used by Kaggle for each of the cell line/drug combinations in the scoring set and test set. The first 108 values are the scoring set (9 cell lines and 12 drugs) and the last 60 are the final test set (5 cell lines 12 drugs). Scores on the final test set will not be shown until the competition is over.

### rand_sub_cont.csv
a sample submission file in the correct format with random predictions between 0 and 1. The calculation of the AUROC value summarizes the performance of these guesses at all thresholds between 0 and 1.

