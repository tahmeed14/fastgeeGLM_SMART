## Function for adding Weights



add.Weights <- function(dataset, option = "known") {
    # Create column for the Weights
    dataset$Weights <- NA;
    
    if (option == "known") {
        dataset[(dataset$R == 1), ]$Weights <- 2;
        dataset[(dataset$R == 0), ]$Weights <- 4;
    }
    
    else if (option == 'estimate') {
        
    }
}