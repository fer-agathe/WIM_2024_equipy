library(dplyr)

setwd("/Users/agathefernandesmachado/Documents/PhD/WIM-concordia/equipy_freMPL/data")
df <- read.csv("SEER_clean_datav1210.csv")
str(df)
df %>% head(10)
summary(df)
nrow(df) ; length(unique(df$patient_ID))
nrow(which(is.na(df), arr.ind=TRUE))/nrow(df)*100
idx <- unique(which(is.na(df), arr.ind=TRUE)[,1])
length(idx)/nrow(df)*100 # delete 0.64%
df <- df[-idx,]
#unique(df$Cod_site)

# Verify statistics of the paper : https://arxiv.org/pdf/2307.13616.pdf

df <- df %>% filter(Metastasis == "M0" & Age_dx >= 18 & Age_dx <= 80 & Surv_mth != 0)
df_surv <- df %>% select(Age_dx, Yr_dx, Surv_mth, Vital_status_rec, Cod_site)
df_surv <- df_surv %>% mutate(Death_melanoma = as.numeric(
  Vital_status_rec == "Dead" & Cod_site == "Melanoma of the Skin")
)
sum(df_surv$Death_melanoma==1)/nrow(df_surv)*100 # 4.9%

new_table <- function(i) {
  new_df <- tibble()
  ind_surv <- df_surv[i,]
  Surv_yr <- ind_surv$Surv_mth/12
  max_duration <- ceiling(Surv_yr)
  for (j in 1:max_duration){
    duration <- j-1
    Age <- ind_surv$Age_dx + duration
    Year <- ind_surv$Yr_dx + duration
    d <- as.numeric(ind_surv$Death_melanoma == 1 & j == max_duration)
    exposure <- ifelse(j < max_duration, 1, 
                       ifelse(d == 0, Surv_yr-duration, 1))
    current_new_df <- tibble(
      i = i,
      Surv_yr = Surv_yr,
      max_duration = max_duration,
      j = j, 
      duration = duration,
      Age = Age,
      Year = Year, 
      d = d,
      exposure = exposure
    )
    new_df <- rbind(new_df, current_new_df)
  }
  new_df
}



#library(parallel)
#num_cores <- detectCores()
#print(num_cores)
# Example using mclapply
#result <- mclapply(1:nrow(df_surv), new_table, mc.cores = 6)
result <- lapply(1:nrow(df_surv), new_table)
combined_results <- do.call(rbind, result)

# suppress all individuals with survival duration > 1 when year of diagnosis > 2018
i_error <- combined_results$i[which(combined_results$Year > 2018)]
i_error <- unique(i_error)
combined_results <- combined_results[-which(combined_results$i %in% i_error),]
df <- df %>% mutate(i = seq_along(patient_ID))
new_tibble <- combined_results %>% left_join(df, by = "i")

str(new_tibble)
to_export <- new_tibble %>% select(d, exposure, duration, Age, Sex, reg_nod_pos, Mitoticrate, 
                                   Laterality, Ulceration, Site_rec_WHO08, Origin, Mar_stat, 
                                   Extent, Surg_primsite, Tumor, Positive_Node)

# Suppress missing origin for this study
to_export <- to_export %>% filter(!(Origin == "Missing"))

write.csv(to_export, file = 'data_mortality.csv', row.names = FALSE)

# Logistic regression
set.seed(123)

# Sample indices
indices <- sample(nrow(to_export), replace = FALSE)

# Train set size
train_size <- round(0.6 * nrow(to_export))

# Sélectionner les indices pour l'ensemble d'entraînement
indices_train <- indices[1:train_size]

to_export$d <- as.factor(d)

# Créer les ensembles d'entraînement et de test
df_train <- to_export[indices_train, ]
df_test <- to_export[-indices_train, ]

# Formula
str <- ""
for (col in colnames(df_train)){
  if (col != "exposure" & col != "d"){
    if (col == "duration"){
      str <- paste(str, col, sep = "")
    } else {
      str <- paste(str, "+", col, sep = "") 
    }
  }
}
#formula <- paste("offset(log(exposure)/(1-log(exposure))) +", str, sep="")
formula <- paste("d ~ ", str, sep="")
formula <- as.formula(formula)

# Fit logistic regression model with weights
lr_model <- glm(formula, data = df_train, family = binomial, weights=exposure)
summary(lr_model)

library(pROC)

# Example data
actual_labels <- df_train$d
predicted_scores <- lr_model$fitted.values

# Create a ROC curve
roc_curve <- roc(actual_labels, predicted_scores)

# Plot the ROC curve
plot(roc_curve, main = "ROC Curve", col = "blue", lwd = 2)

# Add diagonal reference line
abline(a = 0, b = 1, lty = 2, col = "red")

# Add AUC (Area Under the Curve) to the plot
auc_value <- auc(roc_curve)
legend("bottomright", legend = paste("AUC =", round(auc_value, 4)), col = "blue", lwd = 2)

# Set threshold
death_percentage <- sum(to_export$d == 1)/nrow(to_export)*100
f <- function(t){
  predicted_labels <- as.numeric(predicted_scores >= t)
  predicted_death <- sum(predicted_labels == 1)/length(predicted_labels)*100
  return(predicted_death - death_percentage)
}

# Use uniroot to find the root in the interval [a, b]
root <- uniroot(f, interval = c(0, 1))
optimal_threshold <- root$root

# Create a confusion matrix
predicted_labels <- as.numeric(predicted_scores >= optimal_threshold)
conf_matrix <- table(Actual = actual_labels, Predicted = predicted_labels)

# Compute acceptance rate (also known as positive predictive value or precision)
acceptance_rate <- conf_matrix[2, 2] / sum(conf_matrix[, 2])

# Compute true positive rate (sensitivity)
true_positive_rate <- conf_matrix[2, 2] / sum(conf_matrix[2, ])

# Compute false positive rate (specificity)
false_positive_rate <- conf_matrix[1, 2] / sum(conf_matrix[1, ])

# Print the results
cat("Acceptance Rate (Precision):", 1-acceptance_rate, "\n")
cat("True Positive Rate (Sensitivity):", 1-true_positive_rate, "\n")
cat("False Positive Rate (Specificity):", 1-false_positive_rate, "\n")

# Test set
# Predict on the test set
predictions <- exp(predict(lr_model, newdata = df_test))/(1+exp(predict(lr_model, newdata = df_test)))

# Export
df_to_export <- tibble(
  scores = predictions,
  y_true = df_test$d,
  y_pred = as.numeric(predictions >= optimal_threshold),
  gender = df_test$Sex,
  origin = df_test$Origin
)
write.csv(df_to_export, file = 'data_mortality_equipy.csv', row.names = FALSE)



