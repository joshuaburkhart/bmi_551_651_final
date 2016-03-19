# BMI 651: Tuning SVM

library(MASS)
library(plyr) #this must be loaded before dplyr
library(dplyr)
library(ggplot2)
library(broom)
library(knitr)
library(magrittr)
library(reshape2)
library(infotheo)
library(stats)
library(car)
library(e1071)
library(bnlearn)
library(elmNN)
library(gbm)
library(stats)
set.seed(2013)

setwd("~/SoftwareProjects/bmi_551_651_final/JoshBeta")

# gene id cols, cell line rows
expression <-
read.table(
"../data/expression.tsv",header = TRUE,sep = "\t",quote = "",row.names =
1,check.names = FALSE,stringsAsFactors = FALSE
)

# subtype cols, cell line rows
subtypes <-
read.table(
"../data/subtypes.txt",header = TRUE,sep = "\t",row.names = 1,check.names =
TRUE
)
subtypes[,1] <- subtypes[,1] %>% as.numeric()

for (subtype_num in 1:max(subtypes[,1]))
{
    zeros <- rep(0,nrow(subtypes))
    zeros[which(subtypes[,1] %>% as.numeric() == subtype_num)] <- 1
    subtypes <- cbind(subtypes,zeros)
}

# Kristen's celltypes
celltypes <- data.frame(
    t(
        read.table("../data/cell_type_classifications.txt",
                    header=TRUE,
                    sep="\t",
                    check.names=FALSE,
                    row.names=1)
    ),
    check.names = FALSE)

for (celltype in 1:ncol(celltypes))
{
    celltypes[,celltype] <- as.numeric(celltypes[,celltype])
}

# drug cols, cell line rows
training.classes <-
t(
read.table(
"../data/training_set_answers.txt",header = TRUE,sep = "\t",check.names = FALSE
)
)
training.classes[,1:ncol(training.classes)] <- training.classes[,1:ncol(training.classes)] %>% as.numeric()

# Usage, id, drug, cell line columns, values in rows
scoring_and_test_set_id_mappings <-
read.table(
"../data/scoring_and_test_set_id_mappings.csv",header = TRUE,sep = ",",check.names =
FALSE
)

sae <- expression

for (subtype_num in 2:(max(subtypes[,1]) + 1))
{
    # organize subtypes, expression, train, test, classes...
    subtype_df <- data.frame(t(subtypes[,subtype_num]),check.names = FALSE)
    colnames(subtype_df) <- colnames(expression)
    sae <- rbind(subtype_df,sae)
}

x.nrow_before_celltypes <- nrow(sae)

sae <- rbind(sae,celltypes)

x.nrow_after_celltypes <- nrow(sae)
x.celltype_rows <- ((x.nrow_before_celltypes+1):x.nrow_after_celltypes)

# split training data and kaggle holdout
training.features <- sae[,colnames(training.classes)]
kaggle.features <- sae[setdiff(colnames(sae),colnames(training.classes))]

print("Loading feature indeces...")
x.all_features <- vector("list",12)
#From Wilcoxon Rank Sum p < 0.01
x.all_features[[1]] <- c(
7883,6065,13969,833,1243,3980,12900,13998,586,2720,3013,3510,5108,12168,13639,14424,17449,17759,431,4544,9660,6097,8704,10815,16297
)
x.all_features[[2]] <- c(519,1444,4537,9768,16316,513,4471,4567,10857,16330,152,357,3284,5083,7418,10153,12351,15522,17562,2945,3255,3750,4220,6672,7625
)
x.all_features[[3]] <-c(3456,4537,13438,8035,17562,8097,17606,18270,14218,18066,988,2614,5732,12061,13022,14477,4145,6860,10177,10799,13509,15909,17631,1064,1334
)
x.all_features[[4]] <-c(4688,6564,12536,4215,4448,6999,10882,12247,3020,16649,9149,9319,14649,2137,3065,6901,6965,8772,15701,5782,6808,9256,10385,10415,16721
)
x.all_features[[5]] <-c(18094,4982,5821,17005,10094,16613,17747,7413,8279,3785,7078,10639,18486,14105,17523,11,1013,1087,2379,10737,11452,12066,12924,12948,14575
)
x.all_features[[6]] <-c(17669,16319,1193,7623,17254,1134,2299,3353,4488,6052,12638,13820,14580,298,3291,3830,9158,11888,12252,12572,16521,7513,8902,9950,11199
)
x.all_features[[7]] <-c(8042,17551,4141,18616,13179,353,12212,5644,9128,15294,3251,1247,2851,4169,11170,12598,13482,72,8159,11653,11984,14054,569,3843,5809
)
x.all_features[[8]] <-c(10222,2798,6566,3864,4616,13108,13934,14911,15066,16031,18252,738,14033,3181,6915,13946,17539,2252,4372,4816,7831,9554,17357,379,2666
)
x.all_features[[9]] <-c(484,6063,9634,15952,11969,16912,5084,12021,16365,17159,18392,3583,6076,7013,8068,10356,16127,16654,715,3288,4159,5684,9834,12461,13080
)
x.all_features[[10]] <-c(5869,9675,14833,18208,16006,2063,8710,12381,17581,6808,11230,13520,436,3877,5544,7395,10594,14890,793,9874,12808,15537,1744,2350,3169
)
x.all_features[[11]] <-c(9451,15958,17021,16230,16928,4561,9111,4292,13620,210,7762,16219,4798,9982,14355,18299,6013,7639,3442,3740,6104,8744,14725,16654,7013
)
x.all_features[[12]] <-c(3201,11654,2701,5520,5805,14662,15838,9131,9646,10273,17534,2832,8712,11696,12874,13836,7428,9162,9715,12695,3776,5889,6790,7221,8619
)

print("done.")

### Training and Kaggle Sets

kaggle.classes <- data.frame()
training.features_train <- data.frame()
training.classes_train <- numeric()
x.cv_err <- Inf
training.features_ncol <- numeric()
x.loop = NA
x.cv_boost <- numeric()
x.all_drug_min_errs <- numeric()
x.SPLIT <- 0.6

for (drug_idx in 1:nrow(training.classes))
{
    training.features_train <- training.features
    training.classes_train <- training.classes[drug_idx,]

    ### Cross Validation Loop

    # variables to store loop results
    x.cv_err <- Inf
    training.features_ncol <- round(x.SPLIT * ncol(training.features_train))
    x.cv_boost <- rep(1 / ncol(training.features_train),ncol(training.features_train))
    x.loop_models <- list(type=any)
    prev_avg <- NA

    for (loop in 1:1000)
    {
        ### Split Training & Validation sets with random sampling (Monte Carlo cross validation)

        # Cross Validation
        index <- sample(1:ncol(training.features_train),round(x.SPLIT * ncol(training.features_train)),replace=FALSE,prob=x.cv_boost)
        training.features_train_cv <- training.features_train[,index]
        training.classes_train_cv <- training.classes_train[index]

        x.fcols <- 1:ncol(training.features_train)
        negative_idxs <- match(sample(x.fcols[-index],round((1-x.SPLIT) * ncol(training.features_train)),replace=FALSE),x.fcols)

        training.features_validation_cv <- training.features_train[,negative_idxs]
        training.classes_validation_cv <- training.classes_train[negative_idxs]

        #now we can retain only our selected columns
        training.features_train_w <-
        training.features_train_cv[c(x.all_features[[drug_idx]],x.celltype_rows),]

        #add the class labels to the feature data frame
        training.train_svm_input <- as.matrix(t(training.features_train_w))
        #data.frame(t(training.features_train_w),check.names = FALSE)

        #filter features
        training.validation_svm_input <-
        data.frame(t(training.features_validation_cv[c(x.all_features[[drug_idx]],x.celltype_rows),]), check.names = FALSE)

        g.denom <- length(c(x.all_features[[drug_idx]],x.celltype_rows))

        #train model
        #x.model <- best.tune(svm,
        #    train.x = training.train_svm_input,
        #    train.y = training.classes_train_cv,
        #    validate.x = training.validation_svm_input,
        #    validate.y = training.classes_validation_cv,
        #    ranges = list(gamma = 2^seq(-3,0), cost = 2^seq(-4,0)),
        #    tunecontrol = tune.control(sampling = "fix"))
        x.model <- svm(
            x = training.train_svm_input,y = training.classes_train_cv,kernel="linear"
        )

        x.loop_models[[loop]] <- x.model

        ### Validation (on hold out data)

        x.prediction_cv <- predict(x.model,training.validation_svm_input)

        x.class_ag <-
        table(round(unlist(lapply(x.prediction_cv,function(x) {ifelse(x < 0,0,ifelse(x > 1,1,x))}))),training.classes_validation_cv) %>% classAgreement()

        #misclassification rate
        x.err_rate <- 1 - x.class_ag$diag

        print(c("loop",loop,"error rate:",x.err_rate))
        x.cv_err <- x.err_rate
        x.loop <- loop

        x.prediction_cv <- round(unlist(lapply(x.prediction_cv,function(x) {ifelse(x < 0,0,ifelse(x > 1,1,x))})))

        ### Boosting probabilities
        miscalled_idxs <- which(abs(round(x.prediction_cv) - training.classes_validation_cv)>0)
        x.cv_boost[(negative_idxs[miscalled_idxs])] <- abs(x.cv_boost[(negative_idxs[miscalled_idxs])] * 1.01)
    }
    x.all_drug_min_errs <- c(x.all_drug_min_errs,x.cv_err)

    print(c("drug:",drug_idx,"loop:",x.loop,'min drug err:', x.cv_err))#,"x.cv_boost:",x.cv_boost,"prev_avg:",prev_avg))

    ### Predict Kaggle Holdout

    #filter features
    kaggle.svm_input <-
    data.frame(t(kaggle.features[c(x.all_features[[drug_idx]],x.celltype_rows),]), check.names = FALSE)

    kaggle.predictions <- (predict(x.loop_models[[1]],kaggle.svm_input) / x.loop)
    for(loop_it in 2:x.loop)
    {
        kaggle.predictions <- (kaggle.predictions + predict(x.loop_models[[loop_it]],kaggle.svm_input) / x.loop)
    }

    #print(c("kaggle.predictions",kaggle.predictions))

    kaggle.classes <- rbind(kaggle.classes,round(unlist(
    lapply(kaggle.predictions,function(x) {
        ifelse(x < 0,0,ifelse(x > 1,1,x))
    })
    )))
}
rownames(kaggle.classes) <- rownames(training.classes)
colnames(kaggle.classes) <- colnames(kaggle.features)
kaggle.classes %>% str()
mapping <-
read.csv("../data/scoring_and_test_set_id_mappings.csv",check.names = FALSE)
kaggle.out <- data.frame()
for (row in 1:nrow(mapping))
{
    kaggle.out <-
    rbind(kaggle.out,c(mapping[row,"id"],kaggle.classes[(mapping[row,"drug"]),(mapping[row,"cellline"])]))
}
colnames(kaggle.out) <- c('id','value')

fptr <- file(description = "./tuningSVM_out.csv",'w')
write.csv(kaggle.out,fptr,row.names = FALSE, quote = FALSE)
close(fptr)

print(
c(
    "all drug min errs:",x.all_drug_min_errs,"avg_err:",sum(x.all_drug_min_errs) / length(x.all_drug_min_errs)
)
)
