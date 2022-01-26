#' ccbiplot
#' @description The function adapts the \code{ggbiplot::ggbiplot()} function to an \code{acca}-object. 
#' @param pcobj           an object returned by acca:cc()
#' @param choices         which PCs to plot
#' @param shapes          which shapes to use when plotting observations
#' @param circle.prob     size of the circle in qchisq(circle.prob, df = 2)
#' @param arrows.color    color of the arrows 
#' @param scale           covariance biplot (scale = 1), form biplot (scale = 0). When scale = 1, the inner product between the variables approximates the covariance and the distance between the points approximates the Mahalanobis distance.
#' @param obs.scale       scale factor to apply to observations
#' @param var.scale       scale factor to apply to variables
#' @param pc.biplot       for compatibility with biplot.princomp()
#' @param groups          optional factor variable indicating the groups that the observations belong to. If provided the points will be colored according to groups
#' @param ellipse         draw a normal data ellipse for each group?
#' @param ellipse.prob    size of the ellipse in Normal probability
#' @param labels          optional vector of labels for the observations
#' @param labels.size     size of the text used for the labels
#' @param alpha           alpha transparency value for the points (0 = transparent, 1 = opaque)
#' @param circle          draw a correlation circle?
#' @param var.axes        draw arrows for the variables?
#' @param varname.size    size of the text for variable names
#' @param varname.adjust  adjustment factor the placement of the variable names, >= 1 means farther from the arrow
#' @param varname.abbrev  whether or not to abbreviate the variable names
#' @param ...             not used yet
#' @import ggplot2
#' @importFrom plyr ddply
#' @return It returns a \code{\link[ggplot2]{ggplot}} object.
#' @export
#' 
#' 
#' 
#' 
##@importFrom ggplot2 xvar yvar varname angle hjust
ccbiplot<-
  function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
            obs.scale = 1 - scale, var.scale = scale, groups = NULL, shapes = NULL,
            ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
            alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
            varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
            arrows.color=NULL,...) 
  {
    
    xvar<-yvar<-varname<-angle<-hjust<-NULL
    stopifnot(length(choices) == 2)
    report="percent"
    n_vars=NULL
   if (inherits(pcobj, "acca")) {
      nobs.factor <- sqrt(nrow(pcobj$scores$xscores))
      d <- pcobj$cor
      u <- sweep(pcobj$scores$xscores, 2, 1/(d * nobs.factor), FUN = "*")
      v <- rbind(pcobj$corr$corr.X.xscores,
                 pcobj$corr$corr.Y.xscores)
      report="correlation"
      n_vars=c(nrow(pcobj$xcoef),nrow(pcobj$ycoef))
    }  else {
      stop("An object of class acca is expected! Nothingh done")
      return(NULL)
    }
    if(is.null(n_vars)) n_vars=c(nrow(v),0)
    
    choices <- pmin(choices, ncol(u))
    df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                                FUN = "*"))
    # v <- sweep(v, 2, d^var.scale, FUN = "*")
    df.v <- as.data.frame(v[, choices])
    names(df.u) <- c("xvar", "yvar")
    names(df.v) <- names(df.u)
    if (pc.biplot) {
      df.u <- df.u * nobs.factor
    }
    r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
    v.scale <- rowSums(v^2)
    df.v <- r * df.v/sqrt(max(v.scale))
    if (obs.scale == 0) {
      u.axis.labs <- paste("standardized PC", choices, sep = "")
    }
    else {
      u.axis.labs <- paste("PC", choices, sep = "")
    }
    if(report=="percent"){
      u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                                100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
    } else if(report=="correlation"){
      u.axis.labs <- paste("CC", choices, sep = "")
      u.axis.labs <- paste(u.axis.labs, sprintf("(%0.3f correlation)", 
                                                pcobj$cor[choices]))
    }
    if (!is.null(labels)) {
      df.u$labels <- labels
    }
    if (!is.null(groups)) {
      df.u$groups <- groups
      df.u$shapes <- shapes 
    }
    if (varname.abbrev) {
      df.v$varname <- abbreviate(rownames(v))
    }
    else {
      df.v$varname <- rownames(v)
    }
    df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
    df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
    g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + 
      xlab(u.axis.labs[1]) + 
      ylab(u.axis.labs[2]) + coord_equal()+
      theme(panel.background = element_rect(fill='gray100'),
            legend.background = element_rect(fill='gray100'))
    if (var.axes) {
      if (circle) {
        theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                  length = 50))
        circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                               sin(theta))
        g <- g + geom_path(data = circle, color = "#005C08", 
                           size = 1/2, alpha = 1/3)
      }
      
      if(is.null(arrows.color)){
        arrows.color=c(rep("#832424",n_vars[1]),
                       rep("#005800",n_vars[2]))
      }
      g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
                                             xend = xvar, yend = yvar), 
                            arrow = arrow(length = unit(1/2, "picas")), 
                            color = arrows.color)
      # rep(("gray20"),4),
      # rep(muted("green"),6))
      # )
    }
    if (!is.null(df.u$labels)) {
      if (!is.null(df.u$groups)) {
        g <- g + geom_text(aes(label = labels, color = groups), 
                           size = labels.size)
      }
      else {
        g <- g + geom_text(aes(label = labels), size = labels.size)
      }
    }
    else {
      if (!is.null(df.u$groups)) {
        ###[TODO]  sizes
        g <- g + geom_point(aes(color = groups,
                                size=1,shape=shapes), alpha = alpha,stroke = 2)+
          theme(legend.title=element_blank())+
          guides(colour = guide_legend("title"),
                 size = "none",
                 shape = guide_legend("title")) 
      }
      else {
        g <- g + geom_point(alpha = alpha)
      }
    }
    if (!is.null(df.u$groups) && ellipse) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
      circle <- cbind(cos(theta), sin(theta))
      ell <- ddply(df.u, "groups", function(x) {
        if (nrow(x) <= 2) {
          return(NULL)
        }
        sigma <- var(cbind(x$xvar, x$yvar))
        mu <- c(mean(x$xvar), mean(x$yvar))
        ed <- sqrt(qchisq(ellipse.prob, df = 2))
        data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                         mu, FUN = "+"), groups = x$groups[1])
      })
      names(ell)[1:2] <- c("xvar", "yvar")
      g <- g + geom_path(data = ell, aes(color = groups, group = groups))
    }
    if (var.axes) {
      if(is.null(df.v$varname)) df.v$varname=1:nrow(df.v)

      g <- g + geom_text(data = df.v, aes(label = varname,
                                          x = xvar, y = yvar, 
                                          angle = angle, hjust = hjust), 
                         color = arrows.color,
                         # rep(muted("blue"),n_vars[2])), 
                         size = varname.size)#,min.segment.length=1,segment.alpha=.1)
    }
    return(g)
  }