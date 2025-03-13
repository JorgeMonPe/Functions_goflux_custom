environment(click.peak_custom_CO2CH4) <- asNamespace("goFlux")
##Is the same function that click.peak2 but I added an argument (shoulder_used) in which you enter the shoulder that you used
#with goFlux() function and in the plots you will the limit of your incubaition and where start the shoulder
#Additionally, this one show both gases (CO2 and CH4) which could help you to decide the best strat and end time for your incubation
click.peak_custom_CO2CH4 <- function (ow.list, gastype = "CO2dry_ppm", sleep = 3, plot.lim = c(380, 
  1000), seq = NULL, warn.length = 60, save.plots = NULL, shoulder_used = 0) 
{
  if (!is.numeric(plot.lim) | length(plot.lim) != 2) {
    stop("'plot.lim' must be numeric and of length 2")
  }
  if (!is.numeric(warn.length)) {
    stop("'warn.length' must be of class numeric")
  }
  else {
    if (warn.length <= 0) 
      stop("'warn.length' must be greater than 0")
  }
  if (!is.null(save.plots)) {
    if (!is.character(save.plots)) 
      stop("'save.plot' must be a character string")
  }
  if (!is.null(seq)) 
    if (!is.numeric(seq)) 
      stop("'seq' must be of class numeric")
  if (is.null(seq)) 
    seq <- 1:length(ow.list)
  if (missing(ow.list)) 
    stop("'ow.list' is required")
  if (!is.list(ow.list)) 
    stop("'ow.list' must be of class list")
  for (ow in seq) {
    if (!is.data.frame(ow.list[[ow]])) {
      stop(paste("The ", toOrdinal(ow), " element in 'ow.list' is of class ", 
        class(ow.list[[ow]])[1], ". All elements in 'ow.list' must be of class data.frame", 
        sep = ""))
    }
  }
  for (ow in seq) {
    if (!any(grepl("\\<UniqueID\\>", names(ow.list[[ow]])))) {
      stop(paste("'UniqueID' is required and was not found in the", 
        toOrdinal(ow), "element of 'ow.list'"))
    }
  }
  for (ow in seq) {
    if (!any(grepl("\\<POSIX.time\\>", names(ow.list[[ow]])))) {
      stop(paste("'POSIX.time' is required and was not found in the", 
        toOrdinal(ow), "element of 'ow.list'"))
    }
    else if (!is.POSIXct(ow.list[[ow]]$POSIX.time)) {
      stop("'POSIX.time' in 'ow.list' must be of class POSIXct")
    }
  }
  if (is.null(gastype)) 
    stop("'gastype' is required")
  else {
    if (!is.character(gastype)) 
      stop("'gastype' must be a character string")
  }
  if (!any(grepl(paste("\\<", gastype, "\\>", sep = ""), c("CO2dry_ppm", 
    "COdry_ppb", "CH4dry_ppb", "N2Odry_ppb", "NH3dry_ppb", 
    "H2O_ppm")))) {
    stop(paste("'gastype' must be of class character and one of the following:", 
      "'CO2dry_ppm', 'COdry_ppm', 'CH4dry_ppb', 'N2Odry_ppb', 'NH3dry_ppb' or 'H2O_ppm'"))
  }
  for (ow in seq) {
    if (!any(grepl(paste("\\<", gastype, "\\>", sep = ""), 
      names(ow.list[[ow]])))) {
      stop("data frames in 'ow.list' must contain a column that matches 'gastype'")
    }
    else if (!is.numeric(Reduce("c", ow.list[[ow]][, gastype]))) {
      stop(paste("The column that matches 'gastype' in each data frame of", 
        "'ow.list' must be of class numeric"))
    }
  }
  if (!is.null(sleep)) {
    if (!is.numeric(sleep)) 
      stop("'sleep' must be of class numeric")
    if (sleep > 10) 
      stop("'sleep' must be shorter than 10 seconds")
    if (sleep < 0) 
      stop("'sleep' cannot be negative")
  }
  flag <- . <- POSIX.time <- identify.error <- rownum <- NULL
  sleeploop <- function(x) {
    p <- proc.time()
    Sys.sleep(x)
    proc.time() - p
  }
  default.par <- par(no.readonly = TRUE)
  on.exit(par(default.par))
  on.exit(Sys.unsetenv("TZ"))
  ow.list.name <- deparse(substitute(ow.list))
  ow.corr.ls <- list()
  plots.ls <- list()
  for (ow in seq) {
    flux.meas <- Reduce("c", ow.list[[ow]][, gastype])
    time.meas <- Reduce("c", ow.list[[ow]][, "POSIX.time"])
    ymax <- max(flux.meas, na.rm = TRUE)
    ymin <- min(flux.meas, na.rm = TRUE)
    ydiff <- ymax - ymin
    ylim <- c(ymin - ydiff * 0.05, ymax + ydiff * 0.05)
    ylim.min <- ifelse(ylim[1] < plot.lim[1], plot.lim[1], 
      ylim[1])
    ylim.max <- ifelse(ylim[2] > plot.lim[2], plot.lim[2], 
      ylim[2])
    #If CO2 is the target
    if (gastype=="CO2dry_ppm"){
      
      # Show timeseries for CO2 and CH4 to provide an overview of the incubation
      flux.meas_CO2 <- Reduce("c", ow.list[[ow]][, "CO2dry_ppm"])
      flux.meas_CH4 <- Reduce("c", ow.list[[ow]][, "CH4dry_ppb"])
      
      
      # Graph limits
      yaxis.limit.max_CO2 <- (max(flux.meas_CO2, na.rm = TRUE) + 0.01*max(flux.meas_CO2, na.rm = TRUE)) #%>% ifelse(. > 1000, 1000, .)
      yaxis.limit.min_CO2 <- (min(flux.meas_CO2, na.rm = TRUE) - 0.01*max(flux.meas_CO2, na.rm = TRUE))
      yaxis.limit.max_CH4 <- (max(flux.meas_CH4, na.rm = TRUE) + 0.01*max(flux.meas_CH4, na.rm = TRUE))
      yaxis.limit.min_CH4 <- (min(flux.meas_CH4, na.rm = TRUE) - 0.01*max(flux.meas_CH4, na.rm = TRUE))


    tryCatch({
      identify.error <- NULL
      dev.new(noRStudioGD = TRUE, width = 14, height = 8)
      plot(flux.meas_CH4 ~ time.meas, col="gray",
           xlab = "Time", ylab = "", xaxt = 'n', yaxt = 'n',   # Suppress x and y axes
           ylim = c(yaxis.limit.min_CH4, yaxis.limit.max_CH4))
      # Añadir el eje Y derecho
      axis(4, at = pretty(c(yaxis.limit.min_CH4, yaxis.limit.max_CH4))) # Eje derecho
      mtext(side=3, "¡¡¡ Focus on black dots, gray ones is CH4 and right axis is CH4 ppb !!!", col = "red")
      mtext("CH4dry [ppb]", side = 4) # Etiqueta del eje derecho
      # Añadir el segundo eje Y y nueva secuencia de puntos
      par(new = TRUE) # Permite superponer gráficos
      plot(flux.meas ~ time.meas, main = paste("Incubation:",unique(ow.list[[ow]]$UniqueID)), 
        xlab = "Time", ylab = gastype, xaxt = "n", ylim = c(ylim.min, 
          ylim.max))
      abline(v = c(min(time.meas)+shoulder_used, max(time.meas)-shoulder_used), col = "red", lty = "dashed")
      time.zone <- attr(time.meas, "tzone")
      Sys.setenv(TZ = time.zone)
      axis.POSIXct(1, at = seq(min(time.meas), max(time.meas), 
        by = "10 secs"), format = "%H:%M:%S")
      rownum <- identify(time.meas, flux.meas, pos = FALSE, 
        n = 2, plot = TRUE, atpen = FALSE, offset = 0.5, 
        tolerance = 0.25)
    }, error = function(e) {
      identify.error <<- e
    })
    }#----end if CO2 is the target
    #If CH4 is the target
    if (gastype=="CH4dry_ppb"){
      
      # Show timeseries for CO2 and CH4 to provide an overview of the incubation
      flux.meas_CO2 <- Reduce("c", ow.list[[ow]][, "CO2dry_ppm"])
      flux.meas_CH4 <- Reduce("c", ow.list[[ow]][, "CH4dry_ppb"])
      
      
      # Graph limits
      yaxis.limit.max_CO2 <- (max(flux.meas_CO2, na.rm = TRUE) + 0.01*max(flux.meas_CO2, na.rm = TRUE)) #%>% ifelse(. > 1000, 1000, .)
      yaxis.limit.min_CO2 <- (min(flux.meas_CO2, na.rm = TRUE) - 0.01*max(flux.meas_CO2, na.rm = TRUE))
      yaxis.limit.max_CH4 <- (max(flux.meas_CH4, na.rm = TRUE) + 0.01*max(flux.meas_CH4, na.rm = TRUE))
      yaxis.limit.min_CH4 <- (min(flux.meas_CH4, na.rm = TRUE) - 0.01*max(flux.meas_CH4, na.rm = TRUE))
      
      
      tryCatch({
        identify.error <- NULL
        dev.new(noRStudioGD = TRUE, width = 14, height = 8)
        plot(flux.meas_CO2 ~ time.meas, col="gray",
             xlab = "Time", ylab = "", xaxt = 'n', yaxt = 'n',   # Suppress x and y axes
             ylim = c(yaxis.limit.min_CO2, yaxis.limit.max_CO2))
        # Añadir el eje Y derecho
        axis(4, at = pretty(c(yaxis.limit.min_CO2, yaxis.limit.max_CO2))) # Eje derecho
        mtext(side=3, "¡¡¡ Focus on black dots, gray ones is CO2 and right axis is CO2 ppm !!!", col = "red")
        mtext("CO2dry [ppm]", side = 4) # Etiqueta del eje derecho
        # Añadir el segundo eje Y y nueva secuencia de puntos
        par(new = TRUE) # Permite superponer gráficos
        plot(flux.meas ~ time.meas, main = paste("Incubation:",unique(ow.list[[ow]]$UniqueID)), 
             xlab = "Time", ylab = gastype, xaxt = "n", ylim = c(ylim.min, 
                                                                 ylim.max))
        abline(v = c(min(time.meas)+shoulder_used, max(time.meas)-shoulder_used), col = "red", lty = "dashed")
        time.zone <- attr(time.meas, "tzone")
        Sys.setenv(TZ = time.zone)
        axis.POSIXct(1, at = seq(min(time.meas), max(time.meas), 
                                 by = "10 secs"), format = "%H:%M:%S")
        rownum <- identify(time.meas, flux.meas, pos = FALSE, 
                           n = 2, plot = TRUE, atpen = FALSE, offset = 0.5, 
                           tolerance = 0.25)
      }, error = function(e) {
        identify.error <<- e
      })
    }#----end if CH4 is the target
    if (isTRUE(identify.error[[1]] == "graphics device closed during call to locator or identify")) {
      dev.flush()
      dev.off()
      warning(paste(ow.list.name, "[[", ow, "]] UniqueID ", 
        unique(ow.list[[ow]]$UniqueID), ": ", identify.error[[1]], 
        sep = ""), call. = FALSE)
      if (!is.null(save.plots)) {
        dev.new()
        par(mfrow = c(1, 2))
        plot(flux.meas ~ time.meas, xlab = "Time", ylab = gastype, 
          xaxt = "n", ylim = c(ylim.min, ylim.max))
        time.zone <- attr(time.meas, "tzone")
        Sys.setenv(TZ = time.zone)
        axis.POSIXct(1, at = seq(min(time.meas), max(time.meas), 
          by = "10 secs"), format = "%H:%M:%S")
        mtext(line = -3, outer = T, cex = 1.5, paste(ow.list.name, 
          "[[", ow, "]] UniqueID ", unique(ow.list[[ow]]$UniqueID), 
          sep = ""))
        par(mar = c(0, 0, 0, 0))
        plot(x = 0:10, y = 0:10, ann = F, bty = "n", 
          type = "n", xaxt = "n", yaxt = "n")
        identify.message <- paste(ow.list.name, "[[", 
          ow, "]] UniqueID ", unique(ow.list[[ow]]$UniqueID), 
          ": ", identify.error[[1]], sep = "")
        mes <- strsplit(strwrap(identify.message, width = 60), 
          "\n")
        for (i in seq(along = mes)) text(mes[[i]], y = 6 - 
          i * 0.5, x = 5)
        plots.ls[[ow]] <- recordPlot()
        dev.off()
      }
    }
    else {
      dev.flush()
      dev.off()
      if (length(rownum) < 2) {
        rownum <- c(1, 1)
      }
      start.time_corr <- time.meas[rownum[1]]
      end.time_corr <- time.meas[rownum[2]]
      flux.flag <- which(between(time.meas, start.time_corr, 
        end.time_corr))
      ow.corr.ls[[ow]] <- ow.list[[ow]] %>% mutate(flag = if_else(row_number() %in% 
        flux.flag, 1, 0)) %>% mutate(Etime = as.numeric(POSIX.time - 
        start.time_corr, units = "secs")) %>% mutate(start.time_corr = start.time_corr, 
        end.time_corr = end.time_corr) %>% mutate(obs.length_corr = as.numeric(end.time_corr - 
        start.time_corr, units = "secs"))
      xmin <- min(round_any(ow.corr.ls[[ow]]$Etime, 30, 
        f = floor)) %>% if_else(. == 0, -15, .)
      xmax <- max(round_any(ow.corr.ls[[ow]]$Etime, 30, 
        f = ceiling)) %>% if_else(. == 0, 15, .)
      xmult <- (xmax + abs(xmin))/30
      ymax2 <- ow.corr.ls[[ow]] %>% filter(flag == 1) %>% 
        select(all_of(gastype)) %>% max()
      ymin2 <- ow.corr.ls[[ow]] %>% filter(flag == 1) %>% 
        select(all_of(gastype)) %>% min()
      ydiff2 <- ymax2 - ymin2
      ylim2 <- c(ymin2 - ydiff2, ymax2 + ydiff2)
      ylim.min2 <- ifelse(ylim2[1] < ylim.min, ylim.min, 
        ylim2[1])
      ylim.max2 <- ifelse(ylim2[2] > ylim.max, ylim.max, 
        ylim2[2])
      dev.new(noRStudioGD = TRUE, width = 14, height = 8)
      plot(flux.meas ~ ow.corr.ls[[ow]]$Etime, col = ow.corr.ls[[ow]]$flag + 
        1, main = paste(unique(ow.corr.ls[[ow]]$UniqueID)), 
        xlab = "Etime", ylab = gastype, xaxp = c(xmin, 
          xmax, xmult), ylim = c(ylim.min2, ylim.max2))
      if (!is.null(sleep)) {
        if (sleep > 0) 
          sleeploop(sleep)
        dev.off()
      }
      if (nrow(filter(ow.corr.ls[[ow]], flag == 1)) < 
        warn.length) {
        warning(paste(ow.list.name, "[[", ow, "]] UniqueID ", 
          unique(ow.corr.ls[[ow]]$UniqueID), ": Number of observations (", 
          nrow(filter(ow.corr.ls[[ow]], flag == 1)), 
          ") smaller than warn.length=", warn.length, 
          sep = ""), call. = FALSE)
      }
      else {
        message(paste(ow.list.name, "[[", ow, "]] UniqueID ", 
          unique(ow.corr.ls[[ow]]$UniqueID), ": Good window of observation", 
          sep = ""))
      }
      if (!is.null(save.plots)) {
        dev.new()
        par(mfrow = c(1, 2))
        plot(flux.meas ~ time.meas, xlab = "Time", ylab = gastype, 
          xaxt = "n", ylim = c(ylim.min, ylim.max))
        time.zone <- attr(time.meas, "tzone")
        Sys.setenv(TZ = time.zone)
        axis.POSIXct(1, at = seq(min(time.meas), max(time.meas), 
          by = "10 secs"), format = "%H:%M:%S")
        mtext(line = -3, outer = T, cex = 1.5, paste(ow.list.name, 
          "[[", ow, "]] UniqueID ", unique(ow.corr.ls[[ow]]$UniqueID), 
          sep = ""))
        plot(flux.meas ~ ow.corr.ls[[ow]]$Etime, col = ow.corr.ls[[ow]]$flag + 
          1, xlab = "Etime", ylab = gastype, xaxp = c(xmin, 
          xmax, xmult), ylim = c(ylim.min2, ylim.max2))
        plots.ls[[ow]] <- recordPlot()
        dev.off()
      }
    }
  }
  if (length(plots.ls) > 0) {
    pdf(file = save.plots, width = 11.6, height = 8.2)
    for (p in 1:length(plots.ls)) {
      if (!is.null(plots.ls[[p]])) 
        print(plots.ls[[p]])
    }
    dev.off()
  }
  ow.corr <- map_df(ow.corr.ls, ~as.data.frame(.x))
  return(ow.corr)
}
environment(click.peak_Jorge) <- environment(goFlux::click.peak2)
