# ==============================================================================
# RMST Plotting Functions
# ==============================================================================

## Plot one or more RMST curves on a common set of axes.
## The function automatically selects the plot title, y-axis limits, and legend
## placement based on the supplied curves.
plot.multiple.rmst <- function(
    titles, rmst.list, cols, main.title = ""
) {
  
  num.types <- length(rmst.list)
  
  # ----------------------------------------------------------------------------
  # Plot Configuration
  # ----------------------------------------------------------------------------
  
  ## Generate a default title when one is not explicitly supplied. A single
  ## curve uses its corresponding title, while multiple curves are presented
  ## as an RMST comparison.
  if (main.title == "") {
    main.title <- paste(titles[1], "RMST")
    
    if (num.types > 1) {
      main.title <- "RMST Comparison"
    }
  }
  
  ## Proceed only when the number of curves, labels, and colors agree and the
  ## number of curves is within the supported plotting range
  if (
    num.types == length(titles) &&
    num.types == length(cols) &&
    num.types > 0 &&
    num.types <= 5
  ) {
    
    # --------------------------------------------------------------------------
    # Determine Plot Limits and Legend Placement
    # --------------------------------------------------------------------------
    
    ## Use the final RMST value to determine an appropriate y-axis range and
    ## position the legend away from the region containing the plotted curves.
    avg.last.point <- mean(
      rmst.list[[1]][-nrow(rmst.list[[1]])]
    )
    
    legend.placement <- "topright"
    
    if (avg.last.point < -0.3) {
      ylim <- c(
        min(rmst.list[[1]][-nrow(rmst.list[[1]])]),
        0
      )
      
      legend.placement <- "bottomleft"
      
    } else if (avg.last.point > 0.3) {
      ylim <- c(
        0,
        max(rmst.list[[1]][-nrow(rmst.list[[1]])])
      )
      
    } else {
      ylim <- c(
        min(rmst.list[[1]][-nrow(rmst.list[[1]])]),
        max(rmst.list[[1]][-nrow(rmst.list[[1]])])
      )
    }
    
    
    ## Apply transparency to the plotted curves so that overlapping RMST
    ## trajectories remain visually distinguishable
    alpha.cols <- adjustcolor(
      cols,
      alpha.f = 0.05
    )
    
    # --------------------------------------------------------------------------
    # Draw RMST Curves
    # --------------------------------------------------------------------------
    
    ## Initialize the plot using the first RMST curve
    matplot(
      t.grid,
      rmst.list[[1]],
      type = "l",
      lty = 1,
      col = alpha.cols[1],
      xlab = "Time",
      ylab = "RMST",
      main = main.title
    )
    
    
    ## Overlay any additional RMST curves on the same axes
    if (length(rmst.list) > 1) {
      for (i in 2:num.types) {
        matplot(
          t.grid,
          rmst.list[[i]],
          type = "l",
          lty = 1,
          col = alpha.cols[i],
          add = TRUE
        )
      }
    }
    
    
    ## Add a horizontal reference line at zero RMST contrast
    abline(
      h = 0,
      col = "black",
      lwd = 3,
      lty = 2
    )
    
    # --------------------------------------------------------------------------
    # Add Legend
    # --------------------------------------------------------------------------
    
    legend(
      legend.placement,
      legend = titles,
      col = cols[1:num.types],
      lwd = 2,
      lty = 1
    )
    
  } else {
    
    ## Report an error when the supplied plotting arguments are inconsistent
    cat("error")
    
    if (length(titles) != num.types) {
      cat("mismatch")
    }
  }
}

# ==============================================================================
# Generate RMST Plots
# ==============================================================================

## Generate the complete set of RMST figures used to compare the estimated
## treatment effects from MS-Cox and case-base models against the marginal
## (true) RMST curve.
##
## When save.plots = TRUE, each figure is written to a separate PDF file.
plot.all.rmst <- function(
    true.rmst, cox.rmst, cb.rmst, save.plots = FALSE
) {
  
  # ----------------------------------------------------------------------------
  # Marginal RMST
  # ----------------------------------------------------------------------------
  
  ## Plot the marginal RMST curve alone to display the target causal estimand
  if (save.plots) {
    pdf(
      "true-rmst.pdf",
      width = 8,
      height = 5
    )
  }
  
  plot.multiple.rmst(
    titles = c("Marginal"),
    rmst.list = list(true.rmst),
    cols = "purple"
  )
  
  if (save.plots) {
    dev.off()
  }
  
  # ----------------------------------------------------------------------------
  # MS-Cox vs. Marginal RMST
  # ----------------------------------------------------------------------------
  
  ## Compare the MS-Cox RMST estimate with the marginal RMST target
  if (save.plots) {
    pdf(
      "cox-rmst.pdf",
      width = 8,
      height = 5
    )
  }
  
  plot.multiple.rmst(
    titles = c(
      "MS-Cox",
      "Marginal"
    ),
    rmst.list = list(
      cox.rmst,
      true.rmst
    ),
    cols = c(
      "red",
      "black"
    ),
    main.title = "MS-Cox RMST"
  )
  
  if (save.plots) {
    dev.off()
  }
  
  # ----------------------------------------------------------------------------
  # Case-Base vs. Marginal RMST
  # ----------------------------------------------------------------------------
  
  ## Compare the IPW case-base RMST estimate with the marginal RMST target
  if (save.plots) {
    pdf(
      "cb-rmst.pdf",
      width = 8,
      height = 5
    )
  }
  
  plot.multiple.rmst(
    titles = c(
      "Case Base",
      "Marginal"
    ),
    rmst.list = list(
      cb.rmst,
      true.rmst
    ),
    cols = c(
      "blue",
      "black"
    ),
    main.title = "Case Base RMST"
  )
  
  if (save.plots) {
    dev.off()
  }
  
  # ----------------------------------------------------------------------------
  # Overall RMST Comparison
  # ----------------------------------------------------------------------------
  
  ## Compare both estimators simultaneously against the marginal RMST target
  if (save.plots) {
    pdf(
      "all-rmst.pdf",
      width = 8,
      height = 5
    )
  }
  
  plot.multiple.rmst(
    titles = c(
      "MS-Cox",
      "Case Base",
      "Marginal"
    ),
    rmst.list = list(
      cox.rmst,
      cb.rmst,
      true.rmst
    ),
    cols = c(
      "red",
      "blue",
      "black"
    )
  )
  
  if (save.plots) {
    dev.off()
  }
}