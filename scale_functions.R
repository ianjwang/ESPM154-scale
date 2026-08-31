standardize <- function(x) as.numeric((x - mean(x)) / sd(x))

add_scale_bar <- function(xlim, ylim) {
  x_range <- diff(xlim)
  y_range <- diff(ylim)
  target_length <- 0.2 * x_range
  candidate_lengths <- c(1, 2, 5, 10, 20, 25, 50, 100)
  bar_length <- candidate_lengths[which.min(abs(candidate_lengths - target_length))]
  
  x0 <- xlim[1] + 0.08 * x_range
  y0 <- ylim[1] + 0.08 * y_range
  tick_height <- 0.02 * y_range
  
  segments(x0, y0, x0 + bar_length, y0, lwd = 3)
  segments(x0, y0 - tick_height, x0, y0 + tick_height, lwd = 2)
  segments(x0 + bar_length, y0 - tick_height,
           x0 + bar_length, y0 + tick_height, lwd = 2)
  text(x0 + bar_length / 2, y0 + 0.05 * y_range,
       labels = paste(bar_length, "km"), cex = 0.85)
}

rescale_landscape <- function(grain_km = 1) {
  if (!is.numeric(grain_km) || length(grain_km) != 1 ||
      is.na(grain_km) || grain_km <= 0) {
    stop("grain_km must be one positive number.")
  }
  data <- landscape
  data$grid_x <- floor((data$x - 1) / grain_km)
  data$grid_y <- floor((data$y - 1) / grain_km)
  out <- aggregate(
    cbind(genetic_diversity, temperature, precipitation) ~ grid_x + grid_y,
    data = data, FUN = mean
  )
  out$x_center <- out$grid_x * grain_km + grain_km / 2
  out$y_center <- out$grid_y * grain_km + grain_km / 2
  out$grain_km <- grain_km
  out
}

plot_landscape <- function(data) {
  if (!is.data.frame(data) ||
      !all(c("genetic_diversity", "temperature", "precipitation",
             "x_center", "y_center", "grain_km") %in% names(data))) {
    stop("data must be an object created by rescale_landscape().")
  }
  
  grain_km <- unique(data$grain_km)
  if (length(grain_km) != 1) {
    stop("data must contain only one spatial grain.")
  }
  
  add_color_bar <- function(values, colors, xlim, ylim) {
    value_range <- range(values, na.rm = TRUE)
    x_range <- diff(xlim)
    y_range <- diff(ylim)
    
    # Position the legend inside the right side of the map.
    bar_left <- xlim[2] - 0.075 * x_range
    bar_right <- xlim[2] - 0.045 * x_range
    bar_bottom <- ylim[1] + 0.20 * y_range
    bar_top <- ylim[1] + 0.80 * y_range
    
    # White backing makes the legend readable over any map color.
    rect(bar_left - 0.018 * x_range,
         bar_bottom - 0.045 * y_range,
         bar_right + 0.095 * x_range,
         bar_top + 0.045 * y_range,
         col = "white", border = "grey40")
    
    n <- length(colors)
    y_edges <- seq(bar_bottom, bar_top, length.out = n + 1)
    for (i in seq_len(n)) {
      rect(bar_left, y_edges[i], bar_right, y_edges[i + 1],
           col = colors[i], border = NA)
    }
    rect(bar_left, bar_bottom, bar_right, bar_top, border = "grey20")
    
    tick_values <- seq(value_range[1], value_range[2], length.out = 3)
    tick_y <- seq(bar_bottom, bar_top, length.out = 3)
    segments(bar_right, tick_y,
             bar_right + 0.010 * x_range, tick_y)
    text(bar_right + 0.018 * x_range, tick_y,
         labels = format(round(tick_values, 2), trim = TRUE),
         adj = 0, cex = 0.68)
  }
  
  plot_map <- function(values, main, palette) {
    map_colors <- hcl.colors(100, palette, rev = TRUE)
    value_range <- range(values, na.rm = TRUE)
    breaks <- seq(value_range[1], value_range[2], length.out = 101)
    color_index <- cut(values, breaks = breaks, labels = FALSE,
                       include.lowest = TRUE)
    
    half_cell <- grain_km / 2
    xlim <- c(min(data$x_center - half_cell), max(data$x_center + half_cell))
    ylim <- c(min(data$y_center - half_cell), max(data$y_center + half_cell))
    
    plot(data$x_center, data$y_center, type = "n", asp = 1, axes = FALSE,
         xlab = "", ylab = "", main = main, xlim = xlim, ylim = ylim)
    rect(data$x_center - half_cell, data$y_center - half_cell,
         data$x_center + half_cell, data$y_center + half_cell,
         col = map_colors[color_index], border = NA)
    box()
    
    add_scale_bar(xlim = xlim, ylim = ylim)
    add_color_bar(values = values, colors = map_colors,
                  xlim = xlim, ylim = ylim)
  }
  
  old_par <- par(mfrow = c(1, 3), mar = c(1.5, 1.5, 3.0, 0.7),
                 oma = c(0, 0, 2.2, 0))
  on.exit(par(old_par))
  
  plot_map(data$genetic_diversity, "Genetic diversity", "YlGnBu")
  plot_map(data$temperature, "Temperature", "Inferno")
  plot_map(data$precipitation, "Precipitation", "BluGrn")
  
  mtext(paste0("Spatial grain: ", grain_km, " km; n = ", nrow(data),
               " cells"), outer = TRUE, cex = 1.1, font = 2)
  
  invisible(data)
}

analyze_landscape <- function(data) {
  if (!is.data.frame(data) ||
      !all(c("genetic_diversity", "temperature", "precipitation",
             "x_center", "y_center", "grain_km") %in% names(data))) {
    stop("data must be an object created by rescale_landscape().")
  }
  
  grain_km <- unique(data$grain_km)
  if (length(grain_km) != 1) {
    stop("data must contain only one spatial grain.")
  }
  
  model_data <- data.frame(
    genetic_diversity = standardize(data$genetic_diversity),
    temperature = standardize(data$temperature),
    precipitation = standardize(data$precipitation)
  )
  model <- lm(genetic_diversity ~ temperature + precipitation,
              data = model_data)
  sm <- summary(model)
  statistics <- data.frame(
    grain_km = grain_km,
    n_cells = nrow(data),
    temperature_beta = unname(coef(model)["temperature"]),
    precipitation_beta = unname(coef(model)["precipitation"]),
    adjusted_r_squared = sm$adj.r.squared,
    temperature_p = sm$coefficients["temperature", "Pr(>|t|)"],
    precipitation_p = sm$coefficients["precipitation", "Pr(>|t|)"]
  )
  print(round(statistics, 3))
  
  old_par <- par(mfrow = c(1, 3), mar = c(4.2, 4.2, 3.2, 1),
                 oma = c(0, 0, 2.5, 0))
  on.exit(par(old_par))
  
  plot(data$precipitation, data$genetic_diversity, pch = 16,
       col = rgb(0.1, 0.4, 0.7, 0.45),
       xlab = "Mean precipitation", ylab = "Mean genetic diversity",
       main = "Diversity and precipitation")
  abline(lm(genetic_diversity ~ precipitation, data = data),
         col = "firebrick", lwd = 2)
  
  plot(data$temperature, data$genetic_diversity, pch = 16,
       col = rgb(0.1, 0.4, 0.7, 0.45),
       xlab = "Mean temperature", ylab = "Mean genetic diversity",
       main = "Diversity and temperature")
  abline(lm(genetic_diversity ~ temperature, data = data),
         col = "firebrick", lwd = 2)
  
  betas <- coef(model)[c("temperature", "precipitation")]
  barplot(betas, names.arg = c("Temperature", "Precipitation"),
          col = c("tomato", "steelblue"),
          ylim = range(c(-0.4, 0.9, betas)),
          ylab = "Standardized coefficient", main = "Multivariate model")
  abline(h = 0, col = "grey35")
  mtext(paste0("Adjusted R² = ", round(statistics$adjusted_r_squared, 2)),
        side = 3, line = 0.2, cex = 0.85)
  mtext(paste0("Spatial grain: ", grain_km, " km; n = ", nrow(data),
               " cells"), outer = TRUE, cex = 1.25, font = 2)
  
  invisible(list(data = data, model = model, statistics = statistics))
}

compare_scales <- function(grains_km = c(1, 5, 10, 20)) {
  if (!is.numeric(grains_km) || length(grains_km) == 0 ||
      any(is.na(grains_km)) || any(grains_km <= 0)) {
    stop("grains_km must contain one or more positive numbers.")
  }
  
  get_statistics <- function(g) {
    data <- rescale_landscape(grain_km = g)
    model_data <- data.frame(
      genetic_diversity = standardize(data$genetic_diversity),
      temperature = standardize(data$temperature),
      precipitation = standardize(data$precipitation)
    )
    model <- lm(genetic_diversity ~ temperature + precipitation,
                data = model_data)
    sm <- summary(model)
    data.frame(
      grain_km = g,
      n_cells = nrow(data),
      temperature_beta = unname(coef(model)["temperature"]),
      precipitation_beta = unname(coef(model)["precipitation"]),
      adjusted_r_squared = sm$adj.r.squared,
      temperature_p = sm$coefficients["temperature", "Pr(>|t|)"],
      precipitation_p = sm$coefficients["precipitation", "Pr(>|t|)"]
    )
  }
  
  results <- do.call(rbind, lapply(grains_km, get_statistics))
  print(round(results, 3))
  matplot(results$grain_km,
          results[, c("temperature_beta", "precipitation_beta")],
          type = "b", pch = c(16, 17), lty = 1, lwd = 2,
          col = c("tomato", "steelblue"),
          xlab = "Spatial grain (km)",
          ylab = "Standardized regression coefficient",
          main = "Environmental relationships across spatial grains")
  abline(h = 0, col = "grey40")
  legend("topleft", legend = c("Temperature", "Precipitation"),
         col = c("tomato", "steelblue"), pch = c(16, 17),
         lty = 1, lwd = 2, bty = "n")
  invisible(results)
}
