# ==============================================================================
# Refactor Simulation Results
# ==============================================================================

## Convert the list of simulation results into a variable-wise structure.
## Each element of the returned list contains the corresponding result
## across all simulation replicates.
refactor <- function(result) {
  
  store.by.var <- lapply(
    names(result[[1]]),
    function(x) {
      sapply(result, `[[`, x)
    }
  )
  
  names(store.by.var) <- names(result[[1]])
  
  return(store.by.var)
}