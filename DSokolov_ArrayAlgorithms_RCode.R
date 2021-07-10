## R Code for "Grouping Algorithms for Informative Array Testing in Disease Surveillance" by David Sokolov ##

## purpose: to simulate three array assignment algorithms in an informative array group testing scheme and compare their performance across disease parameters ##

# load necessary libraries #

library(dplyr)
library(purrr)
library(tidyr)
library(genGT)
library(ggplot2)
library(RColorBrewer)

# create a collection of popSize individuals, each with an associated disease risk probability drawn from a beta distribution with desired mean risk (meanRisk) and standard deviation (stdDev) #

create_population <- function(popSize, meanRisk, stdDev) {
  patientID <- c(1:popSize)
  n <- (meanRisk*(1-meanRisk))/(stdDev^2)
  alpha <- meanRisk*n
  beta <- (1-meanRisk)*n
  risk <- rbeta(popSize,alpha,beta)
  population <- as.data.frame(cbind(patientID,risk), stringsAsFactors = FALSE)
  return(population)
}

# use each individual risk probability to assign each individual a disease status (0 = not infected, 1 = infected) #

create_disease_statuses <- function(population) {
  diseasedPopulation <- mutate(population, diseaseStatus = as.numeric(rbernoulli(nrow(population),as.numeric(risk))))
  return(diseasedPopulation)
}

# randomly partition population into array groups (of size 9,16,25, 36, 49, or 64) for testing- if population size (n) is not a multiple of array group size (k), set aside the highest risk n mod k individuals for individual testing (these are denoted by arrayGroup '0') #

partition_arrays_random <- function(population, groupSize) {
  if (groupSize %in% c(9,16,25,36,49,64) == FALSE) {
    return("Group size must be 9, 16, 25, 36, 49, or 64!")
  }
  else {
    individualTestGroup <- slice_max(population, as.numeric(risk), n=nrow(population) %% groupSize)
    arrayTestGroup <- setdiff(population, individualTestGroup)
    numberOfGroups <- nrow(arrayTestGroup) / groupSize
    groupIDs <- rep.int(1:numberOfGroups, groupSize)
    partitionedArrayPopulation <- mutate(arrayTestGroup, arrayGroup = sample(groupIDs, nrow(arrayTestGroup), replace = FALSE))
    individualTestGroup <- mutate(individualTestGroup, arrayGroup = 0)
    partitionedPopulation <- bind_rows(partitionedArrayPopulation,individualTestGroup)
    return(partitionedPopulation)
  }
  
}

# partition population into array groups (of size 9, 16, 25, or 36) using the 'dispersed risk' method- put riskiest individual in group 1, next riskiest individual into group 2, etc.- if population size (n) is not a multiple of array group size (k), set aside the highest risk n mod k individuals for individual testing (these are denoted by arrayGroup '0')#

partition_arrays_dispersedrisk <- function(population, groupSize) {
  if (groupSize %in% c(9,16,25,36,49,64) == FALSE) {
    return("Group size must be 9, 16, 25, 36, 49, or 64!")
  }
  else {
    individualTestGroup <- slice_max(population, as.numeric(risk), n=nrow(population) %% groupSize)
    arrayTestGroup <- setdiff(population, individualTestGroup)
    sortedArrayTestGroup <- arrange(arrayTestGroup, desc(risk))
    numberOfGroups <- nrow(arrayTestGroup) / groupSize
    groupIDs <- rep.int(1:numberOfGroups, groupSize)
    partitionedArrayPopulation <- mutate(sortedArrayTestGroup, arrayGroup = groupIDs)
    individualTestGroup <- mutate(individualTestGroup, arrayGroup = 0)
    partitionedPopulation <- bind_rows(partitionedArrayPopulation,individualTestGroup)
    return(partitionedPopulation)
  }
}

# partition population into array groups (of size 9, 16, 25, or 36) using the 'concentrated risk' method- fill group 1 with riskiest individuals, group 2 with next riskiest individuals, etc.- if population size (n) is not a multiple of array group size (k), set aside the highest risk n mod k individuals for individual testing (these are denoted by arrayGroup '0')#

partition_arrays_concentratedrisk <- function(population, groupSize) {
  if (groupSize %in% c(9,16,25,36,49,64) == FALSE) {
    return("Group size must be 9, 16, 25, 36, 49, or 64!")
  }
  else {
    individualTestGroup <- slice_max(population, as.numeric(risk), n=nrow(population) %% groupSize)
    arrayTestGroup <- setdiff(population, individualTestGroup)
    sortedArrayTestGroup <- arrange(arrayTestGroup, desc(risk))
    numberOfGroups <- nrow(arrayTestGroup) / groupSize
    groupIDs <- rep(1:numberOfGroups, each=groupSize)
    partitionedArrayPopulation <- mutate(sortedArrayTestGroup, arrayGroup = groupIDs)
    individualTestGroup <- mutate(individualTestGroup, arrayGroup = 0)
    partitionedPopulation <- bind_rows(partitionedArrayPopulation,individualTestGroup)
    return(partitionedPopulation)
  }
}

# given a group of individuals (of size 9, 16, 25, 36, 49, or 64), construct a testing array using the gradient method #

construct_testing_array_gradient <- function(testingGroup) {
  sortedTestingGroup <- arrange(testingGroup, desc(risk))
  sidelength <- sqrt(nrow(testingGroup))
  testingArray <- matrix(nrow=sidelength, ncol=sidelength)
  k <- 1
  for(j in 1:sidelength) {
    for(i in 1:sidelength){
      testingArray[i,j] <- sortedTestingGroup[k,]$diseaseStatus
      k <- k+1
    }
  }
  return(testingArray)
}

# given a group of individuals (of size 9, 16, 25, 36, 49, or 64), construct a testing array using the spiral method #

construct_testing_array_spiral <- function(testingGroup) {
  sortedTestingGroup <- arrange(testingGroup, desc(risk))
  sidelength <- sqrt(nrow(testingGroup))
  testingArray <- matrix(nrow=sidelength, ncol=sidelength)
  testingArray[1,1] <- sortedTestingGroup[1,]$diseaseStatus
  k <- 2
  for (i in 2:sidelength) {
    j <- 1
    while (j < i) {
      testingArray[i,j] <- sortedTestingGroup[k,]$diseaseStatus
      k <- k+1
      j <- j+1
    }
    while (j > 0) {
      testingArray[j,i] <- sortedTestingGroup[k,]$diseaseStatus
      k <- k+1
      j <- j-1
    }
  }
  return(testingArray)
}

# given a constructed array, "perform" array testing and calculate the number of tests used and the number of samples decoded in the first round of testing #

test_array <- function(testingArray) {
  arraySize <- nrow(testingArray)
  testsRound1 <- nrow(testingArray) + ncol(testingArray)
  testsRound2 <- 0
  numPosRows <- 0
  numPosCols <- 0
  decodedFirstRound <- 0
  for (j in 1:arraySize) {
    for (i in 1: arraySize) {
      if (1 %in% testingArray[i,] & 1 %in% testingArray[,j]) {
        numPosRows <- numPosRows + 1
        numPosCols <- numPosCols + 1
        testsRound2 <- testsRound2 + 1
      } 
      else {
        decodedFirstRound <- decodedFirstRound + 1
      }
    }
  }
  if (numPosRows == numPosCols & numPosRows == 1) {
    testsRound2 <- 0
    decodedFirstRound <- (arraySize)^2
  }
  totalTests <- testsRound1 + testsRound2
  return(c(totalTests, testsRound1, testsRound2, decodedFirstRound))
}

# Given a partitioned population and an array construction method, construct the arrays and test the population #

implementArrayTesting <- function(partitionedPopulation, constructionMethod) {
  totalTests <- 0
  totalDecodedFirstRound <- 0
  individualTests <- filter(partitionedPopulation, arrayGroup == 0)
  totalTests <- totalTests + nrow(individualTests)
  for (i in 1:max(partitionedPopulation$arrayGroup)) {
    arrayGroupTemp <- filter(partitionedPopulation, arrayGroup == i)
    if (constructionMethod == "gradient") {
      arrayTemp <- construct_testing_array_gradient(arrayGroupTemp)
      arrayTestResultsTemp <- test_array(arrayTemp, master_pooling = master_pooling)
      totalTests <- totalTests + arrayTestResultsTemp[1]
      totalDecodedFirstRound <- totalDecodedFirstRound + arrayTestResultsTemp[4]
    }
    if (constructionMethod == "spiral") {
      arrayTemp <- construct_testing_array_spiral(arrayGroupTemp)
      arrayTestResultsTemp <- test_array(arrayTemp, master_pooling = master_pooling)
      totalTests <- totalTests + arrayTestResultsTemp[1]
      totalDecodedFirstRound <- totalDecodedFirstRound + arrayTestResultsTemp[4]
    }
  }
  return(c(totalTests, totalDecodedFirstRound))
}

## Implementation of complete testing workflow: relative performance of three array assignment algorithms across disease prevalence range ##

## 10 sims each, prevs and sizes below, population size 5000 ##

testingWorkflow <- function(numSimulations, populationSize) {
  prevalences <- c(0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.4, 0.5)
  matrixSizes <- c(9, 16, 25, 36, 49, 64)
  nameVector <- c("Prevalence", "ArraySize", "TotalTestsRandomGradient", "TotalTestsRandomSpiral", "TotalTestsConcentratedGradient", "TotalTestsConcentratedSpiral", "TotalTestsDispersedGradient", "TotalTestsDispersedSpiral", "DecodedRandomGradient", "DecodedRandomSpiral", "DecodedConcentratedGradient", "DecodedConcentratedSpiral", "DecodedDispersedGradient", "DecodedDispersedSpiral")
  testingWorkflowPart1Data <- data.frame(t(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0)))
  colnames(testingWorkflowPart1Data) <- nameVector
  for (prev in prevalences) {
    for (matrixSize in matrixSizes) {
      
      totalTestsAcrossSimsRandomGradient <- c()
      totalTestsAcrossSimsRandomSpiral <- c()
      totalTestsAcrossSimsConcGradient <- c()
      totalTestsAcrossSimsConcSpiral <- c()
      totalTestsAcrossSimsDispGradient <- c()
      totalTestsAcrossSimsDispSpiral <- c()
      
      totalDecodedFirstRoundAcrossSimsRandomGradient <- c()
      totalDecodedFirstRoundAcrossSimsRandomSpiral <- c()
      totalDecodedFirstRoundAcrossSimsConcGradient <- c()
      totalDecodedFirstRoundAcrossSimsConcSpiral <- c()
      totalDecodedFirstRoundAcrossSimsDispGradient <- c()
      totalDecodedFirstRoundAcrossSimsDispSpiral <- c()
      
      for (sim in 1:numSimulations){
        
        population <- create_population(populationSize, prev, stdDev = prev/2)
        population <- create_disease_statuses(population)
        
        populationPartitionedRandom <- partition_arrays_random(population, matrixSize)
        populationPartitionedConcentrated <- partition_arrays_concentratedrisk(population, matrixSize)
        populationPartitionedDispersed <- partition_arrays_dispersedrisk(population, matrixSize)
        
        arrayTestResultsTempRandomGradient <- implementArrayTesting(populationPartitionedRandom, "gradient", master_pooling)
        arrayTestResultsTempRandomSpiral <- implementArrayTesting(populationPartitionedRandom, "spiral", master_pooling)
        arrayTestResultsTempConcGradient <- implementArrayTesting(populationPartitionedConcentrated, "gradient", master_pooling)
        arrayTestResultsTempConcSpiral <- implementArrayTesting(populationPartitionedConcentrated, "spiral", master_pooling)
        arrayTestResultsTempDispGradient <- implementArrayTesting(populationPartitionedDispersed, "gradient", master_pooling)
        arrayTestResultsTempDispSpiral <- implementArrayTesting(populationPartitionedDispersed, "spiral", master_pooling)
        
        totalTestsAcrossSimsRandomGradient <- append(totalTestsAcrossSimsRandomGradient, arrayTestResultsTempRandomGradient[1])
        totalTestsAcrossSimsRandomSpiral <- append(totalTestsAcrossSimsRandomSpiral, arrayTestResultsTempRandomSpiral[1])
        totalTestsAcrossSimsConcGradient <- append(totalTestsAcrossSimsConcGradient, arrayTestResultsTempConcGradient[1])
        totalTestsAcrossSimsConcSpiral <- append(totalTestsAcrossSimsConcSpiral, arrayTestResultsTempConcSpiral[1])
        totalTestsAcrossSimsDispGradient <- append(totalTestsAcrossSimsDispGradient, arrayTestResultsTempDispGradient[1])
        totalTestsAcrossSimsDispSpiral <- append(totalTestsAcrossSimsDispSpiral, arrayTestResultsTempDispSpiral[1])
        
        totalDecodedFirstRoundAcrossSimsRandomGradient <- append(totalDecodedFirstRoundAcrossSimsRandomGradient, arrayTestResultsTempRandomGradient[2])
        totalDecodedFirstRoundAcrossSimsRandomSpiral <- append(totalDecodedFirstRoundAcrossSimsRandomSpiral, arrayTestResultsTempRandomSpiral[2])
        totalDecodedFirstRoundAcrossSimsConcGradient <- append(totalDecodedFirstRoundAcrossSimsConcGradient, arrayTestResultsTempConcGradient[2])
        totalDecodedFirstRoundAcrossSimsConcSpiral <- append(totalDecodedFirstRoundAcrossSimsConcSpiral, arrayTestResultsTempConcSpiral[2])
        totalDecodedFirstRoundAcrossSimsDispGradient <- append(totalDecodedFirstRoundAcrossSimsDispGradient, arrayTestResultsTempDispGradient[2])
        totalDecodedFirstRoundAcrossSimsDispSpiral <- append(totalDecodedFirstRoundAcrossSimsDispSpiral, arrayTestResultsTempDispSpiral[2])
      }
      avgTotalTestsRandomGradient <- mean(totalTestsAcrossSimsRandomGradient)
      avgTotalTestsRandomSpiral <- mean(totalTestsAcrossSimsRandomSpiral)
      avgTotalTestsConcGradient <- mean(totalTestsAcrossSimsConcGradient)
      avgTotalTestsConcSpiral <- mean(totalTestsAcrossSimsConcSpiral)
      avgTotalTestsDispGradient <- mean(totalTestsAcrossSimsDispGradient)
      avgTotalTestsDispSpiral <- mean(totalTestsAcrossSimsDispSpiral)
      
      avgDecodedFirstRoundRandomGradient <- mean(totalDecodedFirstRoundAcrossSimsRandomGradient)
      avgDecodedFirstRoundRandomSpiral <- mean(totalDecodedFirstRoundAcrossSimsRandomGradient)
      avgDecodedFirstRoundConcGradient <- mean(totalDecodedFirstRoundAcrossSimsConcGradient)
      avgDecodedFirstRoundConcSpiral <- mean(totalDecodedFirstRoundAcrossSimsConcSpiral)
      avgDecodedFirstRoundDispGradient <- mean(totalDecodedFirstRoundAcrossSimsDispGradient)
      avgDecodedFirstRoundDispSpiral <- mean(totalDecodedFirstRoundAcrossSimsDispSpiral)
      
      tableData <- as.data.frame(t(c(prev, matrixSize, avgTotalTestsRandomGradient, avgTotalTestsRandomSpiral, avgTotalTestsConcGradient, avgTotalTestsConcSpiral, avgTotalTestsDispGradient, avgTotalTestsDispSpiral, avgDecodedFirstRoundRandomGradient, avgDecodedFirstRoundRandomSpiral, avgDecodedFirstRoundConcGradient, avgDecodedFirstRoundConcSpiral, avgDecodedFirstRoundDispGradient, avgDecodedFirstRoundDispSpiral)))
      colnames(tableData) <- nameVector
      
      testingWorkflowData <- rbind(testingWorkflowData, tableData)
      print(paste("prev", prev, "matrix size", matrixSize))
    }
  }
  return(testingWorkflowData)
}

## plots were created using dplyr/tidyr and ggplot2 (code not shown) ##
