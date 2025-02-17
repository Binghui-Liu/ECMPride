####################################################################################
# author: Binghui Liu
# creation date: 2020-02-12
# file description: Predict whether a specified human protein is ECM protein, 
#                   and calculate the corresponding probability.
#                   Run from the command line.
# version: R-3.6.1
####################################################################################

########## Receive command line arguments ##########
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
absolute.path <- strsplit(script.name, "\\\\ECMPride_parallel.R")[[1]][1]  # Get the unzip path for the ECMPride package
setwd(absolute.path) # Set the unzip path of the ECMPride package to the current working path

args.arg.name <- "--args"
where.args.begin <- grep(args.arg.name, initial.options)
where.proteins.to.be.predicted <- initial.options[where.args.begin + 1]
proteins.to.be.predicted <- read.csv(file = where.proteins.to.be.predicted) # Read the ID of the protein to be predicted

result.path <- dirname(where.proteins.to.be.predicted)

cat("\n\n\nThe prediction is in progress. Please Press Enter once and wait for the result...\n\n")

########## Load required packages ##########
suppressMessages(isInstalled <- require(randomForest))
if (!isInstalled) {
  install.packages("randomForest")
  suppressMessages(require(randomForest))
}

suppressMessages(isInstalled <- require(plyr))
if (!isInstalled) {
  install.packages("plyr")
  suppressMessages(require(plyr))
}

suppressMessages(isInstalled <- require(dplyr))
if (!isInstalled) {
  install.packages("dplyr")
  suppressMessages(require(dplyr))
}

suppressMessages(isInstalled <- require(xlsx))
if (!isInstalled) {
  install.packages("xlsx")
  suppressMessages(require(xlsx))
}

suppressMessages(isInstalled <- require(mRMRe))
if (!isInstalled) {
  install.packages("mRMRe")
  suppressMessages(require(mRMRe))
}

suppressMessages(isInstalled <- require(caret))
if (!isInstalled) {
  install.packages("caret")
  suppressMessages(require(caret))
}

suppressMessages(isInstalled <- require(parallel))
if (!isInstalled) {
  install.packages("parallel")
  suppressMessages(require(parallel))
}
########## load required function
GetPredictionResult <- function(x) {
  if (x[1] > x[2]) {
    result <- c(x[1], "ECM")
  } else {
    result <- c(x[2], "nonECM")
  }
}

GenerateSubModel = function(training.sub.models, test.data.feature)
{
  # Uses random forest to model and predict whether the tested protein was ECM protein
  #  
  # Args:  
  #   part.of.data.feature.ECM: An ECM training dataset containing the specific set of features.  
  #   part.of.data.feature.nonECM: A nonECM training dataset containing the specific set of features.  
  #   test.data.feature: The proteins to be predicted and their features set.  
  # 
  # Returns:  
  #   Protein prediction results: whether ECM protein or not, and the corresponding probability.
  #                               1 represents ECM, 2 represents nonECM.
  #   Training model: the model in this run that generate the prediction result
  library(randomForest)
  set.seed(1)
  # test.of.ECMPride.rf <- predict(training.sub.models, test.data.feature)
  # test.of.ECMPride.rf <- data.frame(test.of.ECMPride.rf)
  test.of.ECMPride.rf.prob <- predict(training.sub.models, test.data.feature, type = "prob")
  # test.of.ECMPride.rf.prob <- data.frame(test.of.ECMPride.rf.prob)
  # prediction.result.this.run <- cbind(test.of.ECMPride.rf, test.of.ECMPride.rf.prob)
  return(test.of.ECMPride.rf.prob)
}

########## Generate the dataset for training (part.of.data.feature)##########
load(file = "materials\\training_features\\training_sub_models.Rdata")
num.sub.model <- 99
cores <- detectCores(logical = TRUE)

########## Generate the dataset for testing (part.of.test.data.feature)##########
load(file = "materials\\full_human_proteins_features\\record_of_PP_human_sp.Rdata")
load(file = "materials\\full_human_proteins_features\\record_of_InterPro_human_sp.Rdata")
load(file = "materials\\full_human_proteins_features\\record_of_GreyPSSM_human_sp.Rdata")
##########load annotations of human proteins
annotation.human.sp <- read.csv(file = "materials\\full_human_proteins_features\\human_full_protein_annotation.csv")

record.of.PP.human.sp.to.be.predicted <- data.frame()
record.of.InterPro.human.sp.to.be.predicted <- data.frame()
record.of.GreyPSSM.human.sp.to.be.predicted <- data.frame()
record.of.annotation.to.be.predicted <- data.frame()

is.successful.prediction <- TRUE
unmatched.ID <- data.frame()

for (i in 1:dim(proteins.to.be.predicted)[1]) {
  temp <- which(as.character(record.of.PP.human.sp[,1]) == as.character(proteins.to.be.predicted[i,1]))
  if (length(temp)!=0) {
    record.of.PP.human.sp.to.be.predicted <- rbind(record.of.PP.human.sp.to.be.predicted, record.of.PP.human.sp[temp,])
    record.of.InterPro.human.sp.to.be.predicted <- rbind(record.of.InterPro.human.sp.to.be.predicted, record.of.InterPro.human.sp[temp,])
    record.of.GreyPSSM.human.sp.to.be.predicted <- rbind(record.of.GreyPSSM.human.sp.to.be.predicted, record.of.GreyPSSM.human.sp[temp,])
    temp.annotation <- which(as.character(annotation.human.sp[,1]) == as.character(proteins.to.be.predicted[i,1]))
    if (length(temp.annotation) != 0){
      record.of.annotation.to.be.predicted <- rbind(record.of.annotation.to.be.predicted, annotation.human.sp[temp.annotation, 2:6])
    } else {
      record.of.annotation.to.be.predicted <- rbind(record.of.annotation.to.be.predicted, annotation.human.sp[1, 2:6])
    }
  } else {
    is.successful.prediction <- FALSE
    unmatched.ID <- rbind(unmatched.ID, as.character(proteins.to.be.predicted[i,1]))
    names(unmatched.ID) <- "unmatched.ID"
  }
}
test.data.feature <- cbind(record.of.PP.human.sp.to.be.predicted,
                           record.of.InterPro.human.sp.to.be.predicted,
                           record.of.GreyPSSM.human.sp.to.be.predicted)
load("materials\\training_features\\IGR_after_selection_mRMR.Rdata")
feature.choose <- IGR[,2]
where.feature <- c()
for (m in 1:length(feature.choose))
{
  where.feature[m] <- which(names(test.data.feature) == feature.choose[m])
}

feature.collection <- test.data.feature[,where.feature]
part.of.test.data.feature <- data.frame(feature.collection=feature.collection)
########## Uses training.sub.models to predict the proteins in part.of.test.data.feature.
cl <- makeCluster(cores)
clusterExport(cl, c("GenerateSubModel","training.sub.models", "part.of.test.data.feature"))
len <- num.sub.model
prediction.result.this.fold.list <- parLapply(cl,  1:len, function(x) {GenerateSubModel(training.sub.models[[x]], part.of.test.data.feature)})
stopCluster(cl)

# prediction.result.submodel.df <- prediction.result.submodel.list[[1]][, 1]
prediction.result.submodel.df.prob <- prediction.result.this.fold.list[[1]][, 1:2]
for (i in 2:length(prediction.result.this.fold.list)) {
  # prediction.result.submodel.df <- cbind(prediction.result.submodel.df, prediction.result.submodel.list[[i]][, 1])
  prediction.result.submodel.df.prob <- cbind(prediction.result.submodel.df.prob, prediction.result.this.fold.list[[i]][, 1:2])
}
temp.ncol <- ncol(prediction.result.submodel.df.prob)
temp1.mean <- apply(as.matrix(prediction.result.submodel.df.prob[,seq(1, temp.ncol, 2)]),1,mean)
temp2.mean <- apply(as.matrix(prediction.result.submodel.df.prob[,seq(2, temp.ncol, 2)]),1,mean)
temp <- cbind(temp1.mean, temp2.mean)
prediction.result <- apply(temp, 1, GetPredictionResult)
prediction.result <- t(prediction.result)

ID <- record.of.PP.human.sp.to.be.predicted[,1]
prediction.result <- cbind(ID,
                           data.frame(prediction.result[,2]),
                           data.frame(as.numeric(prediction.result[,1])))
names(prediction.result) <- c("ID", "prediction.result", "prediction.result.prob")
prediction.result.annotation <- cbind(prediction.result, record.of.annotation.to.be.predicted)
names(prediction.result.annotation) <- c("ID", "prediction.result", "prediction.result.prob",
                                         "Gene.name", "label.human.protein.atlas", 
                                         "label.ExoCarta", "label.GO", "label.GO.term")
write.csv(prediction.result.annotation, file = "proteins_to_be_predicted\\prediction_result.csv", row.names = FALSE)

if (is.successful.prediction)
{
  cat("\n\n\nPrediction Succeed.\n")
  print(prediction.result)
  
  result.path <- file.path(result.path, "prediction_result.csv")
  cat(paste("You can find the full results in the file of '", result.path, "'",sep = ""),"\n")
} else {
  cat("\n\n\nSome proteins are not successfully predicted:\n")
  for (i in 1:dim(unmatched.ID)[1])
  {
    cat(as.character(unmatched.ID[i,1]), "\n")
  }
  
  cat("Possible reason: Not the normal UniProt ID.\n")
  cat("Please make sure all protein names have been converted into UniProt ID.\n")
  cat("\n\nThe other proteins are successfully predicted:\n")
  print(prediction.result)
  
  result.path <- file.path(result.path, "PredictionResult.csv")
  cat(paste("\nYou can find the full results in the file of '", result.path, "'",sep = ""),"\n")
}