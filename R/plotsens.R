plotsens <- function(tau.res, zetaz, zetat, tau1, coeff0){
  g = ggplot(tau.res, aes_string(zetaz, zetat)) +
    stat_contour(aes_string(z = tau1, colour = "..level..")) +
    stat_contour(aes(z = t), colour = "red", breaks = c(-1.96,1.96)) +
    xlab("Coef. on U in model for treatment") +
    ylab("Coef. on U in model for response") +
    theme_bw() +
    annotate("text", x = 0, y = 0, label = coeff0)
  print(direct.label(g, method="bottom.pieces"))
}
