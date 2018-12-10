
tutorial.data <- read.table("TutorialDatasetWideFormat.txt",header=TRUE,na.strings=".");

## Acknowledgement:
## The R code for the functions addWeights, convertToLong, addStages, and addReplicates
## was informed/motivated by the Tutorial provided by the d3 Lab, University of Michigan
## on Penn State The Methodoly Center Website
## Link: https://methodology.psu.edu/publications/materials-for-smart-longitudinal-analysis
## Please read the tutorial on this website for an in-depth review of SMART Data Analysis

# Helper Function #1 - addWeights(dataset = x)

# Requires: Column names for Randomization strictly has to be 'R'
# Modifies: x
# Effects: Reads in raw SMART dataset and assigns appropriate Weights column
addWeights <- function(dataset) {
    # Create column for the Weights
    dataset["Weights"] <- NA;
    
    # Two arm SMART, so R can only be 1 or 0
    dataset[(dataset$R == 1), ]$Weights <- 2;
    dataset[(dataset$R == 0), ]$Weights <- 4;
    
    print("Weights successfully added to SMART dataset")
    return(dataset)
}

tutorial.data <- addWeights(tutorial.data)

# Helper Function #2 - convertToLong(dataset = x)

# Requires: Column names must be appropriate as specified in help document (report)
# Modifies: x
# Effects: Reads in SMART dataset w/ Weights and converts to Long Format 
convertToLong <- function(dataset) {
    
    Time1 <- data.frame(cbind(dataset[ ,c("SubjectID", "A1", "R", "A2", "Weights")], 
                              Time = 1, Y = dataset[ ,"Y1"]))
    Time2 <- data.frame(cbind(dataset[ ,c("SubjectID", "A1", "R", "A2", "Weights")], 
                              Time = 2, Y = dataset[ ,"Y2"]))
    Time3 <- data.frame(cbind(dataset[ ,c("SubjectID", "A1", "R", "A2", "Weights")],
                              Time = 3, Y = dataset[ , "Y3"]))
    
    datasetLONG <- rbind(Time1, Time2, Time3)
    
    #Order the Dataset by SubjectID and Time
    print('Ordering dataset by SubjectID and Time')
    datasetLONG <- datasetLONG[order(datasetLONG["SubjectID"], datasetLONG["Time"]), ]
    
    print('SMART data converted to long format')
    
    return(datasetLONG)
}

tutorial.data <- convertToLong(tutorial.data)

# Helper Function #3 - addStages(dataset = x)

# Requires: Column names must be appropriate as specified in help document (report)
# Modifies: x
# Effects: Reads in long format SMART data and adds Stage 1 and Stage 2 columns to dataset

addStages <- function(dataset) {
    # S1 refers to Stage 1 treatment
    # S1 should be 1 by Time 2 and 3, 0 otherwise
    dataset["S1"] <- NA
    dataset$S1[which(dataset$Time == 1)] <- 0
    dataset$S1[which(dataset$Time == 2)] <- 1
    dataset$S1[which(dataset$Time == 3)] <- 1
    
    # S2 refers to Stage 2 treatment
    # S2 should be 1 by Time 3, 0 otherwise
    dataset["S2"] <- NA
    dataset$S2[which(dataset$Time == 1)] <- 0
    dataset$S2[which(dataset$Time == 2)] <- 0
    dataset$S2[which(dataset$Time == 3)] <- 1
    
    print("Successfully added Stages to SMART dataset")
    return(dataset)
}

tutorial.data <- addStages(tutorial.data)

# Helper Function #3 - addStages(dataset = x, n = # of times outcome measured)

# Requires: Column names must be appropriate as specified in help document (report)
# Modifies: x
# Effects: Reads in long format SMART data w/ Stages and replicates appropriate rows

addReplicates <- function(dataset, n = 3) {
    dataset["wave"] <- NA
    dataset$wave <- dataset$Time
    
    rowsReplicate <- dataset[which(dataset$R == 1), ]
    rowsNoReplicate <- dataset[which(dataset$R == 0), ]
    
    temp.pos1 <- rowsReplicate
    temp.pos1$A2 <- 1
    
    temp.neg1 <- rowsReplicate
    temp.neg1$A2 <- -1
    temp.neg1$wave <- temp.neg1$Time + n
    
    ## combine the dataset with appropriate replicates
    finaldata <- rbind(temp.pos1, temp.neg1, rowsNoReplicate)
    
    # Order by SubjectID and wave
    finaldata <- finaldata[order(finaldata$SubjectID, finaldata$wave), ]
    
    print("Successfully add replicated rows to SMART data")
    print("Data is now ready for Analysis")
    return(finaldata)
}

tutorial.data <- addReplicates(tutorial.data, 3)


