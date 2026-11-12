# Ordinal CAD severity prediction with class-imbalance resampling

This repository contains the R scripts and Python notebook for a study of how class-imbalance resampling affects prediction of coronary artery disease (CAD) severity. The workflow uses clinical and genetic predictors to estimate three cumulative SYNTAX score thresholds and combines the estimates into four ordered severity classes.

## Study design

| Class | SYNTAX score | Severity group |
| --- | --- | --- |
| 0 | 0 | No CAD |
| 1 | Greater than 0 and less than 23 | Low severity |
| 2 | At least 23 and less than 32 | Intermediate severity |
| 3 | At least 32 | High severity |

The three binary tasks are CAD presence (`CADP`, score > 0), CAD severity (`CADS`, score ≥ 23), and high CAD severity (`CADHS`, score ≥ 32). Each task compares logistic regression, k-nearest neighbors, a radial-kernel support vector machine, and random forest. The training conditions are the original class distribution (`orig`), random undersampling (`rus`), random oversampling (`ros`), Edited Nearest Neighbours (`enn`), SMOTE (`smote`), and SMOTE followed by ENN (`smote_enn`).

The scripts evaluate held-out folds and include repeated 10-fold experiments. For the four-class evaluation, the pool-adjacent-violators algorithm makes the three cumulative probabilities nonincreasing. Differences between adjacent probabilities yield class probabilities, and the class with the highest probability is the prediction. Binary evaluation includes balanced accuracy and ROC AUC. Four-class evaluation includes macro balanced accuracy, macro recall, weighted kappa, and error in predicted class rank.

## Explore the saved results

The clinical and genetic data are restricted and are excluded from the repository. The checked-in R archives let you inspect results without access to those data. Clone the repository and run the examples from its root.

```bash
git clone https://github.com/dtrygoni/cad-severity-ordinal-resampling.git
cd cad-severity-ordinal-resampling
```

In R, load the four-class results and preview fold-level performance. This example uses base R.

```r
results <- new.env()
load("Results/SK_ready_fourClass_res_f.Rda", envir = results)
head(results$dataset[c("Algorithm", "Technique", "Fold",
                       "Macro.BalancedAccuracy", "Kappa2")])
```

| Saved files | Contents |
| --- | --- |
| [`Results/Binary_CAD_*.Rda`](Results/) | Fold-level evaluation results for the three binary tasks |
| [`Results/SK_ready_*.Rda`](Results/) | Binary and four-class metric tables formatted for Scott–Knott comparisons |
| [`Models/Models_*_whole.Rda`](Models/) | Saved models from the whole-dataset training branch. Fold-based evaluation expects separately generated `Models_*.Rda` files. |

To rerun the saved-result comparisons, install the packages loaded at the top of [`7.Scott_Knott_script.R`](7.Scott_Knott_script.R), create the output directory, and run the script from the repository root.

```bash
mkdir -p Results/ScottKnott
Rscript 7.Scott_Knott_script.R
```

## Rerun the experiments with authorized data

The preprocessing script expects `Data/Whole_dataset.Rda` to contain an object named `dataset` with `SYNTAX.SCORE` and the clinical and genetic columns referenced by the scripts. Use R for preprocessing, modeling, and analysis, and Python with Jupyter for resampling. The notebook imports NumPy, pandas, scikit-learn, imbalanced-learn, and scikit-learn-extra. Install the R packages listed in each script's `library()` calls and namespace references. Dependency versions are unpinned.

Create `Train_test_sets/` and, for each repeat from 1 to 10, the `Train_test_sets/`, `Models/`, and `Results/` subdirectories under `Repeats/Rep<N>/`. Start each stage from the repository root. Before running preprocessing in a fresh R session, move the benchmark setup and `record_preparation()` definition above the first fold loop, keeping the benchmark output writes at the end. In the model training repeat loop, resolve each repeat path from the repository root or restore the root directory after each iteration.

| Stage | Script | What it does |
| --- | --- | --- |
| 1. Prepare folds | [`1.Train_set_construction.R`](1.Train_set_construction.R) | Builds train/test CSVs for each threshold using SNP quality control, univariate filtering, missing-value imputation, scaling, and Boruta feature selection. |
| 2. Resample training data | [`2.Resampling_script.ipynb`](2.Resampling_script.ipynb) | Creates resampled training CSVs. Update the Windows-style paths for your system and point the repeated section at each repeat directory. The notebook names `Rep10` explicitly. |
| 3. Fit models | [`3.Model_training.R`](3.Model_training.R) | Trains the classifiers for each fold and resampling condition. Run the fold-training blocks for evaluation. The whole-dataset blocks require additional `Whole_set_*.csv` inputs. |
| 4. Evaluate | [`4.Evaluation_script_BinaryClassifiers.R`](4.Evaluation_script_BinaryClassifiers.R), [`5.Multiclass_evaluation.R`](5.Multiclass_evaluation.R) | Scores the binary models and combines their probabilities to evaluate the four ordered classes. |
| 5. Compare results | [`6.Test_set_preparation_for_SK.R`](6.Test_set_preparation_for_SK.R), [`7.Scott_Knott_script.R`](7.Scott_Knott_script.R) | Formats the metric tables and runs statistical comparisons. |

The supporting scripts generate [descriptive tables](8.Descriptive_table_construction.R), [study plots](9.Plots_MDPI.R), [SHAP plots](10.Shap_plots.R), and [additional analyses](11.Complementary_scripts.R). Check their input paths before running them. /`.

## License

The code is available under the [MIT License](LICENSE).
