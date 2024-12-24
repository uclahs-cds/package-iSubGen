## Generalized Feature Selection ################################################
#' Feature Selection Function
#'
#' @param group1.identifiers Character vector of IDs for group 1.
#' @param group2.identifiers Character vector of IDs for group 2.
#' @param feature.dataframe Data frame where rows represent features and columns represent samples.
#' @param selection.method String specifying the method for feature selection ('t.test', 'wilcox.test', 'randomforest', 'LASSO').
#' @return A character vector of significant feature names.
#' @examples
#' feature_selection(group1.identifiers, group2.identifiers, feature.dataframe, selection.method = 't.test')
extract.features <- function(group1.identifiers,
    group2.identifiers,
    feature.dataframe,
    selection.method = 't.test') {
    # Validate inputs
    if (!is.data.frame(feature.dataframe)) {
        stop('feature.dataframe must be a valid data frame where rows are features and columns are samples.')
    }
    if (!is.character(group1.identifiers) || !is.character(group2.identifiers)) {
        stop('group1.identifiers and group2.identifiers must be character vectors representing sample IDs.')
    }
    if (!selection.method %in% c('t.test', 'wilcox.test', 'randomforest', 'LASSO')) {
        stop('Unsupported selection.method. Choose from "t.test", "wilcox.test", "randomforest", or "LASSO".')
    }

    # Initialize the output
    significant.features <- NULL

    # Perform feature selection based on method
    if (selection.method %in% c('t.test', 'wilcox.test')) {
        # Compute p-values for each feature using the specified statistical test
        p.values <- apply(feature.dataframe, 1, function(feature.row) {
            if (var(feature.row[group1.identifiers]) == 0 || var(feature.row[group2.identifiers]) == 0) {
                return(1) # No variance in one of the groups, p-value = 1
            } else {
                if (selection.method == 't.test') {
                    return(t.test(feature.row[group1.identifiers], feature.row[group2.identifiers])$p.value)
                } else if (selection.method == 'wilcox.test') {
                    return(wilcox.test(feature.row[group1.identifiers], feature.row[group2.identifiers])$p.value)
                }
            }
        })

        # Adjust p-values for multiple testing
        adjusted.p.values <- p.adjust(p.values, method = 'BH')

        # Create a data frame with adjusted p-values and feature names
        p.value.dataframe <- data.frame(
            adjusted.p.value = adjusted.p.values,
            feature.name = rownames(feature.dataframe)
        )

        # Order features by adjusted p-value
        p.value.dataframe <- p.value.dataframe[order(p.value.dataframe$adjusted.p.value), ]

        # Select top 100 significant features
        significant.features <- p.value.dataframe$feature.name[p.value.dataframe$adjusted.p.value < 0.05][1:100]
    } else if (selection.method == 'randomforest') {
        # Transpose the feature data frame for random forest input
        transposed.dataframe <- data.frame(t(feature.dataframe))
        transposed.dataframe$group <- ifelse(rownames(transposed.dataframe) %in% group1.identifiers, 1, 0)

        # Fit a Random Forest model
        random.forest.model <- randomForest(
            x = transposed.dataframe[, colnames(transposed.dataframe) != 'group'],
            y = as.factor(transposed.dataframe$group),
            importance = TRUE
        )

        # Extract feature importance and sort by Gini importance
        importance.dataframe <- data.frame(importance(random.forest.model))
        importance.dataframe <- importance.dataframe[order(importance.dataframe$MeanDecreaseGini, decreasing = TRUE), ]

        # Select top 100 significant features
        significant.features <- rownames(importance.dataframe)[importance.dataframe$MeanDecreaseGini > 0.01][1:100]
    } else if (selection.method == 'LASSO') {
        # Transpose the feature data frame for LASSO input
        transposed.dataframe <- data.frame(t(feature.dataframe))
        transposed.dataframe$group <- ifelse(rownames(transposed.dataframe) %in% group1.identifiers, 1, 0)

        # Fit a LASSO model using cross-validation
        lasso.model <- cv.glmnet(
            x = as.matrix(transposed.dataframe[, colnames(transposed.dataframe) != 'group']),
            y = transposed.dataframe$group,
            alpha = 1
        )

        # Extract coefficients from the LASSO model
        coefficient.matrix <- as.matrix(coef(lasso.model, s = lasso.model$lambda.min))
        coefficient.dataframe <- data.frame(
            feature.name = rownames(coefficient.matrix),
            coefficient.value = as.numeric(coefficient.matrix)
        )

        # Filter features with non-zero coefficients
        coefficient.dataframe <- coefficient.dataframe[coefficient.dataframe$coefficient.value != 0, ]

        # Order by absolute coefficient value and select top 100 significant features
        coefficient.dataframe <- coefficient.dataframe[order(abs(coefficient.dataframe$coefficient.value), decreasing = TRUE), ]
        significant.features <- coefficient.dataframe$feature.name[coefficient.dataframe$coefficient.value > 0.01][1:100]
    }

    return(significant.features)
}