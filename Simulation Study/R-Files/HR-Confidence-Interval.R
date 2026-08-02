# ==============================================================================
# Plot Hazard Ratio Confidence Intervals
# ==============================================================================

## Plot a random sample of simulation-specific hazard ratio estimates and their
## corresponding confidence intervals. The true HR and null value are included
## as reference lines to visually assess estimator accuracy and interval
## coverage.
plot.hr.ci <- function(
    true.value,
    title,
    ci.list,
    hr.list,
    colour = "red",
    coverage = 0,
    seed = 123
) {
  
  # ----------------------------------------------------------------------------
  # Select Simulation Replicates
  # ----------------------------------------------------------------------------
  
  ## Set the seed to ensure that the same simulation replicates are selected
  ## each time the plot is generated
  set.seed(seed)
  x = 1:100
  sampled.indices <- sample(1:length(hr.list), 100)
  
  # ----------------------------------------------------------------------------
  # Construct Confidence-Interval Plot
  # ----------------------------------------------------------------------------
  
  ## Plot the confidence intervals as shaded bands, with the corresponding
  ## hazard ratio estimates overlaid as individual simulation trajectories
  plot <- ggplot() +
    geom_ribbon(
      aes(
        x = x,
        ymin = ci.list[[1]][sampled.indices],
        ymax = ci.list[[2]][sampled.indices]
      ),
      fill = colour,
      alpha = 0.5,
      show.legend = FALSE
    ) +
    geom_line(
      aes(
        x = x,
        y = hr.list[sampled.indices],
        colour = title
      )
    ) +
    
    ## Add the true marginal HR as a horizontal reference line
    geom_hline(
      aes(
        yintercept = true.value,
        colour = "True HR",
        linetype = "True HR"
      )
    ) +
    
    ## Add the null HR of one as a reference for hypothesis testing
    geom_hline(
      aes(
        yintercept = 1,
        colour = "Null value",
        linetype = "Null value"
      )
    ) +
    
    # --------------------------------------------------------------------------
  # Configure Plot Appearance
  # --------------------------------------------------------------------------
  
  ## Assign distinct colours to the estimator, true HR, and null value
  scale_colour_manual(
    name = "",
    values = c(
      setNames(colour, title),
      "True HR" = "black",
      "Null value" = "purple"
    )
  ) +
    
    ## Use distinct line types to differentiate the estimator from the
    ## reference values
    scale_linetype_manual(
      name = "",
      values = c(
        title = "solid",
        "True HR" = "dashed",
        "Null value" = "dotted"
      )
    ) +
    
    ## Label the axes and include the empirical coverage in the plot title
    labs(
      x = "Simulation",
      y = "Hazard ratio",
      title = paste(
        title,
        "HR Confidence Intervals - Coverage:",
        coverage
      )
    ) +
    
    theme_bw()
  
  # ----------------------------------------------------------------------------
  # Return Plot
  # ----------------------------------------------------------------------------
  
  ## Return the ggplot object so that it can be displayed, modified, or combined
  ## with other plots
  return(plot)
}

# ==============================================================================
# Plot HR Confidence Intervals for Both Estimators
# ==============================================================================

## Generate and combine confidence-interval plots for the MS-Cox and Case-Base
## estimators, allowing direct visual comparison of their interval coverage and
## estimated hazard ratios across simulation replicates.
plot.all.hr.ci <- function(
    true.hr,
    cox.hr, cox.hr.stats,
    cb.hr, cb.hr.stats,
    seed = 123,
    save.plot = FALSE
) {
  
  # ----------------------------------------------------------------------------
  # Generate MS-Cox Confidence-Interval Plot
  # ----------------------------------------------------------------------------
  
  ## Extract the MS-Cox confidence intervals and estimated HRs for plotting
  cox.hr.ci.plot <- plot.hr.ci(
    true.hr,
    title = "MS-Cox",
    list(
      cox.hr.stats$ci[, "lower"],
      cox.hr.stats$ci[, "upper"]
    ),
    cox.hr,
    colour = "red",
    coverage = cox.hr.stats$summary["coverage"],
    seed = seed
  )
  
  # ----------------------------------------------------------------------------
  # Generate Case-Base Confidence-Interval Plot
  # ----------------------------------------------------------------------------
  
  ## Extract the Case-Base confidence intervals and estimated HRs for plotting
  cb.hr.ci.plot <- plot.hr.ci(
    true.hr,
    title = "Case Base",
    list(
      cb.hr.stats$ci[, "lower"],
      cb.hr.stats$ci[, "upper"]
    ),
    cb.hr,
    colour = "blue",
    coverage = cb.hr.stats$summary["coverage"],
    seed = seed
  )
  
  # ----------------------------------------------------------------------------
  # Combine Estimator-Specific Plots
  # ----------------------------------------------------------------------------
  
  ## Stack the two estimator-specific plots vertically and combine their legends
  combined_plot <- (cox.hr.ci.plot / cb.hr.ci.plot) +
    plot_layout(guides = "collect")
  
  # ----------------------------------------------------------------------------
  # Optionally Save Combined Plot
  # ----------------------------------------------------------------------------
  
  ## Save the combined figure as a PDF when requested
  if (save.plot) {
    ggsave(
      "hr-ci.pdf",
      combined_plot,
      width = 8,
      height = 8
    )
  }
}