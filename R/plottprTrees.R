#' @title Plot taper curve for an object of class \code{tprTrees}
#' @description creating a plot of the taper curve of a tree, over or under bark
#' @param x an object of class 'tprTrees'
#' @param bark either NULL or logical; if TRUE taper curve over bark is plotted,
#' if FALSE taper curve under bark is plotted; if NULL, both are plotted
#' @param col.bark color to be used for plot of bark, if plot of taper curve
#' over and under bark is requested
#' @param obs should observations (measured/observed diameters) be added to
#' the plot?
#' @param assort assortments produced by \code{tprAssortment(, value="merge")}
#' @param legend logical, if legend should be added
#' @param ... further arguments for \code{plot} and \code{points}
#' @details plots the taper curve of a tree. Either over bark or under bark, or
#' both. Elements design can partly be chosen. If assortments are given, these
#' are added to the plot. Doing that, the assortment bottom and top position is
#' indicated by a vertical line and mid-diameter is shown as a point with
#' vertical dashed line. N.B. the mid-diameter shown is under bark and rounded
#' downwards for 0.5 cm if mid-diameter < 20 and for 0.75 cm if bigger. Volume
#' is calculated using this diameter. Reason for that behaviour is that
#' assortment information with regard to diameter and volume reflects the legal
#' rules for roundwood assortments (german RVR).
#' Additionally, assortment names are indicated.
#' One can provide assortment names in a column of \code{assort} named
#' 'assortname', which will be used if available, otherwise the 'Sort'-column
#'  will be used. See Examples.
#' @return No return value, called for side effects
#' @return NULL
#' @examples
#' ## plotting the taper curve of a tree
#' oldpar <- par()
#' par(mfrow = c(1, 1))
#' tree <- tprTrees(spp=1L, Dm=40, Hm=1.3, H=35)
#' plot(tree, type = "l", las = 1, legend = TRUE)
#' plot(tree, bark = TRUE, las = 1)
#' plot(tree, bark = FALSE, las = 1, obs=TRUE) # obs incl. bark!!!
#' tree <- tprTrees(spp=c(1, 1), Dm = c(40, 35), Hm=c(1.3, 1.3), H = c(35, 30))
#' plot(tree, bark = FALSE, las = 1, legend = TRUE) # both trees are plotted
#' plot(tree, bark = TRUE, las = 1, legend = TRUE, obs=TRUE)
#'
#' tree <- tprTrees(spp=1L, Dm=c(40, 32), Hm=c(1.3, 10.5), H=35)
#' plot(tree, type = "l", las = 1, legend = TRUE, obs=TRUE)
#'
#' ## if monotonicity is not forced:
#' tree <- tprTrees(spp=3L, Dm=8, Hm=1.3, H=10)
#' plot(tree, type = "l", las = 1, obs=TRUE, mono=FALSE)
#' plot(tree, type = "l", las = 1, obs=TRUE, mono=TRUE) # default
#'
#' tree <- tprTrees(spp=c(1, 8), Dm = c(40, 40), Hm=c(1.3, 1.3), H = c(35, 35))
#' plot(tree, bark = NULL, las = 1, col.bark = "blue", legend = TRUE)
#' plot(tree, bark = NULL, las = 1, col.bark = "blue", legend = TRUE, obs = TRUE)
#' plot(tree[1, ], main = tprSpeciesCode(spp(tree[1, ]), out = "long"))
#' plot(tree[2, ], main = tprSpeciesCode(spp(tree[2, ]), out = "scientific"))
#' par(mfrow = c(1, 2))
#' plot(tree, bark = TRUE, las = 1)
#'
#' ## now add assortments into taper curve
#' par(mfrow = c(1, 1))
#' pars <- parSort(n=length(tree), Lxh=1, fixN=2, fixL=4, fixA=10)
#' ass <- tprAssortment(tree, pars=pars)
#' plot(tree, assort = ass)
#' plot(tree, bark = FALSE, assort = ass)
#' plot(tree, bark = FALSE, assort = ass, legend = TRUE)
#' plot(tree[1, ], assort = ass[ass$tree == 1, ], main = "first tree in subset")
#' plot(tree[2, ], assort = ass[ass$tree == 2, ], main = "second tree in subset")
#'
#' ## adding own assortment labels using column 'assortname'
#' ass$assortname <- ifelse(grepl("fix", ass$sort), paste0("Fix:", ass$length), ass$sort)
#' plot(tree, assort = ass)
#' par(oldpar)
#' @importFrom graphics abline lines plot points text
#' @export

plot.tprTrees <- function(x, bark = NULL, col.bark = NULL, obs = FALSE,
                          assort = NULL, legend = FALSE, ...) {
  ## catch call
  mc <- match.call()
  mc[[1]] <- as.name("plot")
  mc$bark <- NULL # remove from call to plot and points
  mc$col.bark <- NULL
  mc$legend <- NULL
  mc$assort <- NULL
  mc$obs <- NULL
  if(!is.null(mc$mono)){
    mono <- mc$mono
    mc$mono <- NULL
  } else {
    mono <- TRUE
  }



  ## for each tree do...
  for (i in seq(along = x)) {
    mci <- mc # for each tree a fresh call object to be manipulated
    if (is.null(mci$main)) {
      mci$main <- paste("spp=", spp(x[i]),
        ", D1=", round(D13(x[i]), 1),
        ", H=", round(Ht(x[i]), 1),
        sep = ""
      )
    }
    col <- ifelse(is.null(mci$col), "black", mci$col)
    if (is.null(mci$col)) mci$col <- col
    if (is.null(mci$xlab)) mci$xlab <- "H [m]"
    if (is.null(mci$ylab)) mci$ylab <- "D [cm]"
    if (is.null(mci$type)) mci$type <- "l"
    ## plotting
    plotx <- seq(0, Ht(x[i]), length.out = floor(Ht(x[i])*10))
    mci$x <- plotx
    ## if bark is NULL, use bark==TRUE
    ploty <- tprDiameterCpp(x[i], Hx = plotx, bark = ifelse(is.null(bark), TRUE, bark),
                            mono = mono)
    mci$y <- ploty
    eval(mci, parent.frame()) # plot the taper curve
    abline(v = 0, h = 0)

    if (is.null(bark)) {
      col.bark <- ifelse(is.null(col.bark), "grey", col.bark)
      mci[[1]] <- as.name("points")
      mci$col <- col.bark
      mci$y <- tprDiameterCpp(x[i], Hx = plotx, bark = FALSE, mono = mono)
      eval(mci, parent.frame()) # add lines for bark taper curve via 'points()'
      if (legend == TRUE) {
        if(isTRUE(obs)){
          lgndobs <- "observations"
        } else {
          lgndobs <- NULL
        }
        legend("topright",
          legend = c("taper curve over bark", "taper curve under bark",
                     lgndobs),
          col = c(col, col.bark, ifelse(isTRUE(obs), col, NA)),
          lty = c(1, 1, NA),
          pch = c(NA, NA, ifelse(isTRUE(obs), 3, NA))
        )
      }
    } else {
      if (legend == TRUE) {
        if(isTRUE(obs)){
          lgndobs <- "observations"
        } else {
          lgndobs <- NULL
        }
        lgnd <- ifelse(bark == TRUE, "taper curve over bark",
                       "taper curve under bark")
        lgnd <- c(lgnd, lgndobs)
        legend("topright",
          legend = lgnd,
          lty = c(1, NA),
          pch = c(NA, ifelse(isTRUE(obs), 3, NA)),
          col = col
          )
      }
    }
    if(isTRUE(obs)){
      ## add given observations into plot
      points(x=Hm(x[i]), y=Dm(x[i]), pch=3)
    }
    if (!is.null(assort)) {
      j <- unique(assort$tree)[i] # select correct assortments for actual tree
      asstmp <- assort[assort$tree == j & assort$length != 0, ]
      if (nrow(asstmp) > 0) {
        htmp <- unique(round(c(asstmp$height, asstmp$height + asstmp$length), 4))
        dtmp <- tprDiameterCpp(x[i], Hx = htmp,
          bark = ifelse(is.null(bark), TRUE, bark), mono = mono)
        invisible(sapply(seq(length(htmp)), function(a) {
          lines(x = rep(htmp[a], 2), y = c(0, dtmp[a]))
        }))
        invisible(sapply(seq(nrow(asstmp)), function(a) {
          lines(
            x = rep(asstmp$height[a] + asstmp$length[a] / 2, 2),
            y = c(0, asstmp$mdm[a]),
            lty = 2, col = "light grey"
          )
        }))
        points(x = asstmp$height + asstmp$length / 2, y = asstmp$mdm, pch = 16)
        if (!is.null(asstmp$assortname)) {
          asslabel <- asstmp$assortname
        } else {
          asslabel <- asstmp$sort
        }
        text(x = asstmp$height + asstmp$length / 2, y = 0, pos = 3, labels = asslabel)
      }
    }
  }
}
