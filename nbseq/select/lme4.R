suppressPackageStartupMessages({
    library('tidyverse')
    library('lme4')
    library('lmerTest')
    library('parallel')
})


args <- commandArgs(TRUE)

input_path  <- args[[1]] # 'intermediate/cdr3/features/top_enrichment_expt/feature_table.csv'
out_path    <- args[[2]]
feature_col <- args[[3]] # 'CDR3ID'
ag_col      <- args[[4]] # 'FliC'
threads     <- as.integer(args[[5]])
# verbose <- as.logical(as.integer(args[[4]]))


df <- read_csv(input_path, show_col_types = FALSE)
columns <- c(feature_col, 'name', 'r', 'abundance', ag_col)
df_ag <- df[, columns, drop=FALSE] %>% rename('Ag'=all_of(ag_col)) %>% arrange(Ag)
# df_ag <- df_ag[order(df_ag[feature_col]),]

fit_model <- function(.feature, data) {
    df_feature <- data[data[feature_col] == .feature, ]

    tryCatch({
        model <- lmer("abundance ~ r * Ag + (r | name)", data=df_feature)
        s <- summary(model)
        pvalues <- s$coefficients[,'Pr(>|t|)']
        names(pvalues) <- paste0(names(pvalues), '_pvalue')
        coeffs <- s$coefficients[,'Estimate']
        val <- c(coeffs, pvalues, 'warnings'=paste(model@optinfo$conv$lme4$messages, collapse=';'))
        # print(model)
        return(val)
    },
    error = function(e) {
        # print(e)
        return(list(warnings=as.character(e)))
    })
}


features <- unique(df_ag[[feature_col]])
# TEMP
# features <- features[1:10]

if (threads == 1) {
    .apply <- lapply
} else {
    if (threads == -1) {
        threads <- future::availableCores()
    }
    .apply <- mclapply
    options('mc.cores' = threads)
    message("Using parallel processing with ", getOption('mc.cores'), " cores...")
}

out <- .apply(features,
            fit_model,
            data = df_ag)

out_df <- bind_rows(out) %>% add_column( "{feature_col}" := features, .before = 1 ) %>% write_csv(out_path)
