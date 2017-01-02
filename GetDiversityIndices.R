# Calculate diversity indices and return data frame of results
# Arguments: matrix=species abundance matrix w/species name as 
#                   column headers, sites as row headers, 
#                   and number of species occuring at site as
#                   cell values                 
GetDiversityIndices <- function(matrix) {

  require(vegan)

  N0 <- rowSums(matrix > 0)      #Species richness
  H <- diversity(matrix)         #Shannon entropy
  N1 <- exp(H)                   #Shannon diversity number
  N2 <- diversity(matrix, "inv") #Simpson diversity number
  J <- H/log(N0)                 #Pielou evenness
  E1 <- N1/N0                    #Shannon evenness (Hill's ratio)
  E2 <- N2/N0                    #Simpson evenness (Hill's ratio)

  #Put all indices in a data frame
  div <- data.frame(N0,H,N1,N2,E1,E2,J) 
    
  # Pass indices data frame back to user
  return(div)
}