plot.bathy2<-function (x, image = FALSE, bpal = NULL, land = FALSE, deepest.isobath, 
          shallowest.isobath, step, n = 20, lwd = 1, lty = 1, col = "black", 
          default.col = "white", drawlabels = FALSE, xlab = "Longitude", 
          ylab = "Latitude", asp = 1, ...) 
{
  mat <- x
  if (!missing("deepest.isobath")) {
    n.levels <- length(deepest.isobath)
  }
  else {
    n.levels <- 1
  }
  lon <- unique(as.numeric(rownames(mat)))
  lat <- unique(as.numeric(colnames(mat)))
  if (land == FALSE) 
    mat[mat > 0] <- 0
  if (image == FALSE) {
    if (n.levels == 1) {
      if (!missing("deepest.isobath") | !missing("shallowest.isobath") | 
          !missing("step")) {
        level.param <- seq(deepest.isobath, shallowest.isobath, 
                           by = step)
      }
      else {
        level.param <- pretty(range(mat, na.rm = TRUE), 
                              n = n)
      }
      contour(lon, lat, mat, levels = level.param, lwd = lwd, 
              lty = lty, col = col, drawlabels = drawlabels, 
              xlab = xlab, ylab = ylab, asp = asp, ...)
    }
    if (n.levels > 1) {
      level.param <- seq(deepest.isobath[1], shallowest.isobath[1], 
                         by = step[1])
      contour(lon, lat, mat, levels = level.param, lwd = lwd[1], 
              lty = lty[1], col = col[1], drawlabels = drawlabels[1], 
              xlab = xlab, ylab = ylab, asp = asp, ...)
      for (i in 2:n.levels) {
        level.param <- seq(deepest.isobath[i], shallowest.isobath[i], 
                           by = step[i])
        contour(lon, lat, mat, levels = level.param, 
                lwd = lwd[i], lty = lty[i], col = col[i], drawlabels = drawlabels[i], 
                add = TRUE)
      }
    }
  }
  if (image == TRUE) {
    if (is.null(bpal)) {
      ramp <- colorRampPalette(c("black","#245372", "#4871D9", 
                                 "#7D86A1", "cornflowerblue","aquamarine","cyan","white"))
      bpal <- ramp(5000)
    }
    if (is.list(bpal)) 
      bpal <- palette.bathy(mat, layers = bpal, land = land, 
                            default.col = default.col)
    if (n.levels == 1) {
      if (!missing("deepest.isobath") | !missing("shallowest.isobath") | 
          !missing("step")) {
        level.param <- seq(deepest.isobath, shallowest.isobath, 
                           by = step)
      }
      else {
        level.param <- pretty(range(mat, na.rm = TRUE), 
                              n = n)
      }
      image(lon, lat, mat, col = bpal, xlab = xlab, ylab = ylab, 
            asp = asp, ...)
      contour(lon, lat, mat, levels = level.param, lwd = lwd, 
              lty = lty, col = col, drawlabels = drawlabels, 
              add = TRUE)
    }
    if (n.levels > 1) {
      image(lon, lat, mat, col = bpal, xlab = xlab, ylab = ylab, 
            asp = asp, ...)
      level.param <- seq(deepest.isobath[1], shallowest.isobath[1], 
                         by = step[1])
      contour(lon, lat, mat, levels = level.param, lwd = lwd[1], 
              lty = lty[1], col = col[1], drawlabels = drawlabels[1], 
              add = TRUE)
      for (i in 2:n.levels) {
        level.param <- seq(deepest.isobath[i], shallowest.isobath[i], 
                           by = step[i])
        contour(lon, lat, mat, levels = level.param, 
                lwd = lwd[i], lty = lty[i], col = col[i], drawlabels = drawlabels[i], 
                add = TRUE)
      }
    }
  }
}
