## Assign Stage Indicators
## This function also converts the Wide format dataset to Long

library(dplyr)

add.Stages(dataset) {
    
    subset.time.1 <- dataset %>% filter(Time == 1) %>% select(SubjectID, )
}