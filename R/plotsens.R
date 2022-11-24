plotsens <- function(tau.res, coeff0, partialRsq = FALSE){
  if (partialRsq){
    # Interpolate points for the contour plot.
    fld <- with(tau.res, interp(x = pR2z, y = pR2t1, z = tau1))
    df <- melt(fld$z, na.rm = TRUE)
    names(df) <- c("x", "y", "tau1")
    df$pR2z <- fld$x[df$x]
    df$pR2t1 <- fld$y[df$y]

    ggplot(data = df, aes(x = pR2z, y = pR2t1)) +
      # geom_tile(aes(fill = tau1)) +
      stat_contour(aes(z = tau1, colour = after_stat(level))) +
      xlab("Partial R^2 of U with the treatment") +
      ylab("Partial R^2 of U with the response") +
      theme_bw() +
      annotate("text", x = 0, y = 0, label = coeff0) +
      geom_text_contour(aes(z = tau1), skip = 0, min.size = 0)
      # scale_fill_continuous(name = "tau1",
      #                      low = "darkblue", high = "white")
  }
  else {
    ggplot(data = tau.res, aes(x = zetaz, y = zetat1)) +
      stat_contour(aes(z = tau1, colour = after_stat(level))) +
      stat_contour(aes(z = t), colour = "red", breaks = c(-1.96,1.96)) +
      xlab("Coef. on U in model for treatment") +
      ylab("Coef. on U in model for response") +
      theme_bw() +
      annotate("text", x = 0, y = 0, label = coeff0) +
      geom_text_contour(aes(z = tau1), skip = 0, min.size = 0)
  }
}
