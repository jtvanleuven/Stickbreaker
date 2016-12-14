#' Create a blank plot to add to
#'
#' @param mars Margins to set in \code{par(mar=mars)}. Default is c(0,0,0,0).
#' @return Nothing
#' @details Creates an blank plot with no margins and xlim=c(0,1) and ylim=c(0,1).
#' @export
#'
#'
empty.plot <- function(mars=c(0,0,0,0)){
  par(mar=mars)
  plot(1,1, type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), bty="n")
}



#' Plot a sector of a cirle
#'
#' @param a X offset
#' @param b Y offset
#' @param r Radius
#' @param p Proportion of circle to fill
#' @param col Fill color
#' @param den Desity of fill (it hatching). Default is NA.
#' @param angle Angle of hatching (if hatching). Default is NA.
#' @param ydistortion Distortion to make circles not streched by plotting
#' @param zero.angle Angle to offset filled sector
#' @param lwd Line width. Default is 1.
#' @param n.pts Number of points to discretize circle into. Default=100.
#' @return Nothing
#' @details Creates an blank plot with no margins and xlim=c(0,1) and ylim=c(0,1).
#' @export
#'
#'

plot.circle.color.sector <- function(a,b,r,p,col, den=NA, angle=NA, ydistortion=1, zero.angle=0, lwd=1, n.pts=100){

  t.v <- seq(0, 2*pi, length.out=n.pts)

  #a <- 1; b <- 1; r <- 0.05
  #p <- 0.33
  #col <- "deepskyblue"

  x <- a + r*cos(t.v)
  y <- b + (r*sin(t.v)*ydistortion)

  points(x,y, type="l")
  base.cut <- round(n.pts*zero.angle)+1
  pt.cut <- base.cut + round(n.pts*p)
  sector.x <- c(a, x[base.cut:pt.cut], a)
  sector.y <- c(b, y[base.cut:pt.cut], b)
  if (p > 0 && p < 1){
    polygon(sector.x, sector.y, col=col, density=den, angle=angle, border="black", lwd=lwd)
  } else if (p==1){
    polygon(sector.x[2:(length(x))], sector.y[2:(length(x))], col=col, density=den, angle=angle, border="black", lwd=lwd)
  }
}








#' Create multipanel figure that compares d esimators
#'
#' @param outpath The pathway (including file name) where to write figure to
#' @param error.table Generated in \code{analyze.stick.data.batch}.
#' @param seq.method.props Generated in \code{analyze.stick.data.batch}
#' @param muts.to.plot Mutation values to plot (values not indecies)
#' @param coes.to.plot Coefficient values to plot
#' @param sigs.to.plot Sigma values to plot
#' @param ylims Vector of upper y-axis limits for [1] Failure rate, [2] Bias and [3] rMSE
#' @param est.to.plot Vector of estimators to plot. Options: "MLE", "RDB.all", "RDB", "max", "seq" (case sensitive).
#' @param col.v Color vector corresponding to est.to.plot
#' @param angle.v Angle vector corresponding to est.to.plot
#' @param den.v Density vector corresponding to est.to.plot
#' @return Nothing
#' @details Generates figure from paper. This probably will not be included in the package. But I'm including it
#' here to keep everything together. Creates a pdf and a svg version.
#' @export

create.d.estimation.figure <- function(outpath, error.table, seq.method.props, mut.vals, coe.vals, sig.vals, ylims, est.to.plot, col.v, angle.v, den.v){

for (plot.i in 1:2){
  if (plot.i==1){
    pdf(file=paste(outpath, ".pdf", sep=""), width=10, height=12)
  } else if (plot.i ==2){
    svg(file=paste(outpath, ".svg", sep=""), width=10, height=12)
  }

lab.cex <- 1.5
layout(mat=matrix(nrow=6, ncol=5, data=c(seq(1, 5),6, seq(7,10),6,seq(11,14), 15, seq(16,19),15, seq(20,23),rep(24,5)), byrow=TRUE), widths=c(0.3,0.3, rep(3, 3)), heights=c(0.25, rep(3,4), 0.75)) -> l
#layout.show(l)
plot.sectors <- TRUE

empty.plot()  # 1
empty.plot()  # 2
empty.plot()  # 3
text(labels="Failure Rate", x=0.5, y=0.5, font=2, cex=lab.cex)
empty.plot()  # 4
text(labels="Bias", x=0.5, y=0.5, font=2, cex=lab.cex)
empty.plot()  # 5
text(labels="rMSE", x=0.5, y=0.5, font=2, cex=lab.cex)

mut.lab.v <- c("3 Mutations (8 Genotypes)", "4 Mutations (16 Genotypes)")

sig.lab.v <- c("sigma=0.02", "sigma=0.08")
hy.cols.v <- c("white", "black", "grey80")
#est.to.plot <- c("d.hat.MLE", "d.hat.median", "d.hat.median.pos", "d.hat.obs.max", "d.hat.hy")
#est.cols <- sapply(est.to.plot, function(x) which(colnames(fail.table)==x))
#est.to.plot <- c("MLE", "RDB.all", "RDB", "max", "seq")
props.to.plot <- c("MLE", "RDB", "max")
p.est.cols <- sapply(props.to.plot, function(x) which(colnames(seq.method.props)==paste("p.", x, sep="")))
legend.names.s <- c("MLE", "RDB.all", "RDB", "max", "seq")
legend.names.all <- c(expression(hat(d)[MLE]), expression(hat(d)[RDB[all]]), expression(hat(d)[RDB]), expression(hat(d)[Max]), expression(hat(d)[Seq]))
which.names <- unlist(sapply(est.to.plot, function(x) which(legend.names.s==x)))
legend.names <- legend.names.all[which.names]


for (mut.i in 1:length(muts.to.plot)){
  empty.plot()
  text(labels=mut.lab.v[mut.i], x=0.5, y=0.5,srt=90, font=2, cex=lab.cex)
  n.muts <- muts.to.plot[mut.i]
  mut.rows <- which(error.table$n.muts == muts.to.plot[mut.i])
  mut.rows2 <- which(seq.method.props$n.muts == muts.to.plot[mut.i])


  for (sig.i in 1:2){
    empty.plot()
    text(labels=bquote(sigma==.(sigs.to.plot[sig.i])), x=0.5, y=0.5,srt=90, font=2, cex=lab.cex)

    # --- Failure ---
    par(mar=c(5,4,1,1))
    yrng <- ylims[1]
    xrng <- 3
    yx.ratio <- (yrng/xrng)*1.25
    plot(1,1, type="n", xlim=c(0.5,3.5), ylim=c(0, yrng), xaxt="n", xlab="", ylab="Failure Rate", yaxt="n")
    axis(side=2, at=seq(0,yrng, by=0.1), labels=seq(0,yrng, by=0.1), cex=0.75)
    axis(side=1, at=seq(1,length(coes.to.plot)), labels=coes.to.plot, cex=0.75)
    abline(0,0)
    mtext(text=expression(u), side=1, line=2)
    sig.rows <- which(error.table$sigma==sigs.to.plot[sig.i])
    sig.rows2 <- which(seq.method.props$sigma==sigs.to.plot[sig.i])
    sig.mut.rows <- intersect(mut.rows, sig.rows)
    sig.mut.rows2 <- intersect(mut.rows2, sig.rows2)
    a <- 1/6

    for (coe.i in 1:length(coes.to.plot)){
      coes.rows <- which(error.table$coes==coes.to.plot[coe.i])
      coes.rows2 <- which(seq.method.props$coes==coes.to.plot[coe.i])
      plot.row <- intersect(coes.rows, sig.mut.rows)
      plot.row2 <- intersect(coes.rows2, sig.mut.rows2)
      subset.data <- error.table[plot.row,]
      row.order <- sapply(est.to.plot, function(x) which(subset.data$method==x))
      new.rows <- plot.row[row.order]

      x.offs <- c(-2.5, -1.5, -0.5, 0.5, 1.5,2.5)

      for (j in 1:length(est.to.plot)){
        rect(coe.i+x.offs[j]*a, 0, coe.i+x.offs[j+1]*a,error.table$fail.rate[new.rows[j]], col=col.v[j], angle=angle.v[j], density=den.v[j])
      }

      #rect(coe.i-2.5*a, 0, coe.i-1.5*a, error.table$fail.rate[new.rows[1]], col=col.v[1])
      #rect(coe.i-1.5*a, 0, coe.i-0.5*a, error.table$fail.rate[new.rows[2]], col=col.v[2])
      #rect(coe.i-0.5*a, 0, coe.i+0.5*a, error.table$fail.rate[new.rows[3]], col=col.v[3])
      #rect(coe.i+0.5*a, 0, coe.i+1.5*a, error.table$fail.rate[new.rows[4]], col=col.v[4])
      #rect(coe.i+1.5*a, 0, coe.i+2.5*a, error.table$fail.rate[new.rows[5]], col="black", angle=45, density=10)
    }

    # --- Bias ---
    par(mar=c(5,4,1,1))
    yrng <- ylims[2]
    xrng <- 3
    yx.ratio <- ((2*yrng)/xrng)*1.25
    plot(1,1, type="n", xlim=c(0.5,3.5), ylim=c(-yrng, yrng), xaxt="n", xlab="", ylab="Bias", yaxt="n")
    axis(side=2, at=seq(-yrng,yrng, by=0.2), labels=seq(-yrng,yrng,by=0.2), cex=0.75)
    axis(side=1, at=seq(1,length(coes.to.plot)), labels=coes.to.plot, cex=0.75)
    abline(0,0)
    mtext(text=expression(u), side=1, line=2)
    sig.rows <- which(error.table$sigma==sigs.to.plot[sig.i])
    sig.mut.rows <- intersect(mut.rows, sig.rows)

    for (coe.i in 1:length(coes.to.plot)){
      coe.rows <- which(error.table$coes==coes.to.plot[coe.i])
      coe.rows2 <- which(seq.method.props$coes==coes.to.plot[coe.i])
      plot.row <- intersect(coe.rows, sig.mut.rows)
      plot.row2 <- intersect(coe.rows2, sig.mut.rows2)
      subset.data <- error.table[plot.row,]
      row.order <- sapply(est.to.plot, function(x) which(subset.data$method==x))
      new.rows <- plot.row[row.order]

      for (j in 1:length(est.to.plot)){
        rect(coe.i+x.offs[j]*a, 0, coe.i+x.offs[j+1]*a,error.table$d.bias[new.rows[j]], col=col.v[j], angle=angle.v[j], density=den.v[j])
      }

      #rect(coe.i-2.5*a, 0, coe.i-1.5*a, error.table$d.bias[new.rows[1]], col=col.v[1])
      #rect(coe.i-1.5*a, 0, coe.i-0.5*a, error.table$d.bias[new.rows[2]], col=col.v[2])
      #rect(coe.i-0.5*a, 0, coe.i+0.5*a, error.table$d.bias[new.rows[3]], col=col.v[3])
      #rect(coe.i+0.5*a, 0, coe.i+1.5*a, error.table$d.bias[new.rows[4]], col=col.v[4])
      #rect(coe.i+1.5*a, 0, coe.i+2.5*a, error.table$d.bias[new.rows[5]], col="black", angle=45, density=10)

      rad <- 0.20
      if (plot.sectors==TRUE){
        plot.circle.color.sector(a=coe.i+2*a, b=yrng-0.15, r=rad, p=1, col="white", ydistortion =yx.ratio, zero.angle=0, lwd=0.5)
        zero.angle <- 0
        for (col.i in 1:length(p.est.cols)){
          #print(zero.angle)
          plot.circle.color.sector(a=coe.i+2*a, b=yrng-0.15, r=rad, p=seq.method.props[plot.row2, p.est.cols[col.i]], col=hy.cols.v[col.i], ydistortion=yx.ratio, zero.angle=zero.angle, lwd=0.5)
          zero.angle <- zero.angle + seq.method.props[plot.row2, p.est.cols[col.i]]
        }
      }
    }

    # --- rMSE ---
    par(mar=c(5,4,1,1))
    yrng <- ylims[3]
    plot(1,1, type="n", xlim=c(0.5,3.5), ylim=c(0,yrng), xaxt="n", xlab="", ylab="rMSE", yaxt="n")
    axis(side=2, at=seq(0,yrng, by=0.2), labels=seq(0,yrng,by=0.2), cex=0.75)
    axis(side=1, at=seq(1,length(coes.to.plot)), labels=coes.to.plot, cex=0.75)
    abline(0,0)
    mtext(text=expression(u), side=1, line=2)
    sig.rows <- which(error.table$sigma==sigs.to.plot[sig.i])
    sig.mut.rows <- intersect(mut.rows, sig.rows)

    for (coe.i in 1:length(coes.to.plot)){
      coe.rows <- which(error.table$coes==coes.to.plot[coe.i])
      plot.row <- intersect(coe.rows, sig.mut.rows)
      subset.data <- error.table[plot.row,]
      row.order <- sapply(est.to.plot, function(x) which(subset.data$method==x))
      new.rows <- plot.row[row.order]

      for (j in 1:length(est.to.plot)){
        rect(coe.i+x.offs[j]*a, 0, coe.i+x.offs[j+1]*a,error.table$d.rMSE[new.rows[j]], col=col.v[j], angle=angle.v[j], density=den.v[j])
      }

      #rect(coe.i-2.5*a, 0, coe.i-1.5*a, error.table$d.rMSE[new.rows[1]], col=col.v[1])
      #rect(coe.i-1.5*a, 0, coe.i-0.5*a, error.table$d.rMSE[new.rows[2]], col=col.v[2])
      #rect(coe.i-0.5*a, 0, coe.i+0.5*a, error.table$d.rMSE[new.rows[3]], col=col.v[3])
      #rect(coe.i+0.5*a, 0, coe.i+1.5*a, error.table$d.rMSE[new.rows[4]], col=col.v[4])
      #rect(coe.i+1.5*a, 0, coe.i+2.5*a, error.table$d.rMSE[new.rows[5]], col="black", angle=45, density=10)
    }
  }
}
empty.plot()
legend("center", 1, legend=c(legend.names),
       fill=col.v, col=col.v,
       ncol=length(est.to.plot), cex=2, angle=angle.v, density=den.v*2.5)
dev.off()
}
}








#' Create multipanel figure showing bias and rMSE when estimating stickbreaking coefficients
#'
#' @param outpath The pathway (including file name) where to write figure to
#' @param error.table Data frame with rMSE, bias and other summary statistics. Generated in \code{analyze.stick.data.batch}.
#' @param seq.method.props Generated in \code{analyze.stick.data.batch}
#' @param muts.to.plot Simulated mutation values to plot.
#' @param coes.to.plot Simulated coefficient values to plot.
#' @param sigs.to.plot Simulated sigma values to plot.
#' @return Nothing
#' @details Generates figure for paper. This probably will not be included in the package. But I'm including it
#' here to keep everything together. Creates a pdf and a svg version.
#' @export


create.stick.coe.estimation.figure <- function(outpath, error.table, muts.to.plot, coes.to.plot, sigs.to.plot){

  for (plot.i in 1:2){
  if (plot.i == 1){
    pdf(file=paste(outpath, ".pdf", sep=""), width=10, height=6.67)
  } else{
    svg(file=paste(outpath, ".svg", sep=""), width=10, height=6.67)
  }

  layout(mat=matrix(nrow=4, ncol=4, data=c(seq(1,12), rep(13, 4)), byrow=TRUE), widths=c(0.5,rep(3, 3)), heights=c(0.25, rep(3,2), 0.75)) -> l
  #layout.show(l)

  sub.table <- error.table[which(error.table$method=="seq"),]
  #offs <- c(-0.4, -0.2, 0, 0.2)
  offs <- c(-0.45, -0.15, +0.15)
  bar.wdth <- 0.3
  col.v <- c("black", "grey50", "white")

  lab.cex <- 1.5
  par(mar=c(5,4,1,1))
  empty.plot()
  empty.plot()
  text(labels="3 Mutations", x=0.5, y=0.5, font=2, cex=lab.cex)
  empty.plot()
  text(labels="4 Mutations", x=0.5, y=0.5, font=2, cex=lab.cex)
  empty.plot()
  text(labels="5 Mutations", x=0.5, y=0.5, font=2, cex=lab.cex)

  # --- Bias ---
  empty.plot()
  text(labels="Bias", x=0.5, y=0.5, srt=90, font=2, cex=lab.cex)
  par(mar=c(5,4,1,1))

  for (mut.i in 1:length(muts.to.plot)){
    plot(1,1, type="n", xlim=c(0.5,3.5), ylim=c(-0.12, 0.12), xaxt="n", xlab="", ylab="Bias")
    axis(side=1, at=seq(1,length(coes.to.plot)), labels=coes.to.plot, cex=0.75)
    abline(0,0)
    mtext(text=expression(u), side=1, line=2)


    mut.val <- muts.to.plot[mut.i]
    mut.rows <- which(sub.table$n.muts == mut.val)

    for (coe.i in 1:3){
      coe.val <- coes.to.plot[coe.i]
      coe.rows <- which(sub.table$coes==coe.val)
      mut.coe.rows <- intersect(mut.rows, coe.rows)
      vals <- sub.table$coe.bias[mut.coe.rows]
      for (sig.i in 1:length(sigs.to.plot)){
        sig.dex <- which(sub.table$sigma == sigs.to.plot[sig.i])
        rect(coe.i+offs[sig.dex], 0, coe.i+offs[sig.dex]+bar.wdth, vals[sig.dex], col=col.v[sig.dex])
      }
    }
  }

  # --- rMSE ---
  empty.plot()
  #text(labels=expression(sqrt(MSE)), x=0.5, y=0.5, srt=90, font=2, cex=lab.cex)
  text(labels="rMSE", x=0.5, y=0.5, srt=90, font=2, cex=lab.cex)

  coes.to.plot <- coe.vals[1:3]
  muts.to.plot <- mut.vals[1:3]
  sigs.to.plot <- sig.vals[1:3]
  #col.v <- c("black", "grey50", "white")
  #offs <- c(-0.4, -0.2, 0, 0.2)


  par(mar=c(5,4,1,1))

  for (mut.i in 1:length(muts.to.plot)){
    plot(1,1, type="n", xlim=c(0.5,3.5), ylim=c(0, 0.15), xaxt="n", xlab="", ylab="rMSE")

    axis(side=1, at=seq(1,length(coes.to.plot)), labels=coes.to.plot, cex=0.75)
    abline(0,0)
    mtext(text=expression(u), side=1, line=2)

    mut.val <- muts.to.plot[mut.i]
    mut.rows <- which(sub.table$n.muts == mut.val)

    for (coe.i in 1:3){
      coe.val <- coes.to.plot[coe.i]
      coe.rows <- which(sub.table$coes==coe.val)
      mut.coe.rows <- intersect(mut.rows, coe.rows)
      vals <- sub.table$coe.rMSE[mut.coe.rows]
      for (sig.i in 1:length(sigs.to.plot)){
        rect(coe.i+offs[sig.i], 0, coe.i+offs[sig.i]+bar.wdth, vals[sig.i], col=col.v[sig.i])
      }
    }
  }

  empty.plot()
  legend(x=0.4, y=0.9, legend=sig.vals, fill=col.v, cex=1.4, title=expression(sigma), ncol=3)

  dev.off()
}  #next plot.i

}  # end create stick.coe.estimation.figure




#' Create figure comparing coefficient estimates under three models
#'
#' @param error.table.list List of error tables from \code{analyze.stick.data.batch} and \code{analyze.mult.add.data.batch}
#' with names "stick", "mult" and "add" (in taht order).
#' @inheritParams create.stick.coe.estimation.figure



create.coe.estimation.figure <- function(outpath, error.table.list, muts.to.plot, coes.to.plot, sigs.to.plot){
  for (plot.i in 1:2){
    if (plot.i == 1){
      pdf(file=paste(outpath, ".pdf", sep=""), width=10, height=6.67)
    } else{
      svg(file=paste(outpath, ".svg", sep=""), width=10, height=6.67)
    }

    #layout(mat=matrix(nrow=3, ncol=3, data=c(seq(1,6), rep(7, 3)), byrow=TRUE), widths=c(3,3,3), heights=c(0.25, 3, 0.75)) -> l
    layout(mat=matrix(nrow=4, ncol=4, data=c(seq(1,12), rep(13, 4)), byrow=TRUE), widths=c(0.5,rep(3, 3)), heights=c(0.25, rep(3,2), 0.75)) -> l
    #layout.show(l)

    offs <- c(-0.45, -0.15, +0.15)
    bar.wdth <- 0.3
    col.v <- c("black", "grey50", "white")

    lab.cex <- 1.5
    par(mar=c(0,0,0,0))
    empty.plot()
    empty.plot()
    text(labels="3 Mutations", x=0.5, y=0.5, font=2, cex=lab.cex)
    empty.plot()
    text(labels="4 Mutations", x=0.5, y=0.5, font=2, cex=lab.cex)
    empty.plot()
    text(labels="5 Mutations", x=0.5, y=0.5, font=2, cex=lab.cex)


    for (sig.i in 1:2){
      sig.val <- sigs.to.plot[sig.i]
      empty.plot()
      text(labels=bquote(sigma==.(sigs.to.plot[sig.i])), x=0.5, y=0.5,srt=90, font=2, cex=lab.cex)

      par(mar=c(5,4,1,1))
      for (mut.i in 1:length(muts.to.plot)){
        plot(1,1, type="n", xlim=c(0.5,3.5), ylim=c(0, 0.15), xaxt="n", xlab="", ylab="rMSE")

        axis(side=1, at=seq(1,length(coes.to.plot)), labels=coes.to.plot, cex=0.75)
        abline(0,0)
        mtext(text=expression(u), side=1, line=2)
        mut.val <- muts.to.plot[mut.i]

        for (coe.i in 1:3){
          coe.val <- coes.to.plot[coe.i]
          stick.dex <- with(error.table.list$stick, which(n.muts==mut.val & coes==coe.val & sigma==sig.val))
          mult.dex <- with(error.table.list$mult, which(n.muts==mut.val & coes==coe.val & sigma==sig.val))
          add.dex <- with(error.table.list$add, which(n.muts==mut.val & coes==coe.val & sigma==sig.val))
          stick.val <- error.table.list$stick$coe.rMSE[stick.dex]
          mult.val <- error.table.list$mult$coe.rMSE[mult.dex]
          add.val <- error.table.list$add$coe.rMSE[add.dex]
          vals <- c(stick.val, mult.val, add.val)
          for (i in 1:3){
            rect(coe.i+offs[i], 0, coe.i+offs[i]+bar.wdth, vals[i], col=col.v[i])
          }
        }   # next coe.i
      } #next mut.i
    }  #next sig.i

    empty.plot()
    legend(x=0.3, y=0.9, legend=c("Stickbreaking", "Multiplicative", "Additive"), fill=col.v, cex=1.4, title="Model", ncol=3)

    dev.off()
  }  #next plot.i
}



#'  Create wireframe panel figure of posterior probabilities from analyzed simulated data
#'
#' @param posterior.list List of posterior probability dataframes generated in
#' \code{\link{summarize.posteriors.on.simulated.dataset}}
#' @param outpath File (with full path) where plot is written
#' @param coes.to.plot Coefficients to plot. If \code{NA} it plots all values in posterior.list.
#' @param sigs.to.plot Sigma values to plot. If \code{NA} it plots all values in posterior.list.
#' @param plot.type Specifies which field of posterior.list to plot on z-axis of wireframe.
#' Options are "mean.post" (mean posterior probabiltiy) or "true.unq" (proportion of time model is only accepted model).
#' Default is "mean.post".
#' @return Nothing. Creates plot (pdf and svg) written to \code{outpath}
#' @export

create.posterior.wireframe.panel <- function(posterior.list, outpath, coes.to.plot=NA, sigs.to.plot=NA, plot.type="mean.post"){
  grey.cols <- c(paste("grey", seq(0,50), sep=""), paste("grey", seq(60, 100, by=5)))
  model.coes <- c(expression(u), expression(s), expression(paste(Delta, w)))
  n.muts <- length(posterior.list)

  model.v <- c("stick", "mult", "add")
  coe.col <- which(colnames(posterior.list[[1]])=="coes")
  sig.col <- which(colnames(posterior.list[[1]])=="sigma")
  coe.v <- sort(unique(posterior.list[[1]]$coe))
  sig.v <- sort(unique(posterior.list[[1]]$sigma))
  if (is.na(coes.to.plot[1])){
    coe.v <- sort(unique(posterior.list[[1]]$coe))
  } else{
    coe.v <- coes.to.plot
  }
  coe.v <- rev(coe.v)
  if (is.na(sigs.to.plot[1])){
    sig.v <- sort(unique(posterior.list[[1]]$sigma))
  } else{
    sig.v <- sigs.to.plot
  }
  plot.list <- vector("list", length(posterior.list))
  for (i in 1:3){
    plot.list[[i]] <- vector("list",3)
  }

  scales.list <- vector("list", 4)
  names(scales.list) <- c("x", "y", "arrows", "col")
  scales.list$x$at <- sig.v
  scales.list$x$labels <- round(sig.v, 2)
  scales.list$y$at <- coe.v
  scales.list$y$labels <- coe.v
  scales.list$arrows <- FALSE
  scales.list$col <- "black"

  col.to.plot <- which(colnames(posterior.list[[1]])==plot.type)
  if (plot.type=="mean.post"){
    zlab <- "Mean posterior"
  } else if (plot.type == "true.unq"){
    zlab <- "Proportion"
  }

  for (mut.i in 1:length(posterior.list)){
    n.muts <- posterior.list[[mut.i]]$n.muts[1]
    for (mod.i in 1:3){
      model <- model.v[mod.i]
      subdata <- posterior.list[[mut.i]][which(posterior.list[[mut.i]]$model==model),]
      data.m <- matrix(nrow=length(coe.v), ncol=length(sig.v))
      for (coe.i in 1:length(coe.v)){
				for (sig.i in 1:length(sig.v)){
					the.row <- which(subdata$coes==coe.v[coe.i] & subdata$sigma==sig.v[sig.i])
					data.m[coe.i, sig.i] <- subdata[the.row, col.to.plot]
				}  #next sig.i
      }   #next coe.i
      d <- t(data.m)
      coe.type <- model.coes[mod.i]
      title <- paste(n.muts, " muts, ", model, sep="")
      plot.list[[mut.i]][[mod.i]] <- wireframe(d, shade=FALSE, light.source=c(0,0,10), col.regions=grey.cols, drape=TRUE, zlim=c(0, 1), row.values=sig.v, column.values=coe.v, col="black", zlab=list(zlab, rot=90), xlab=expression(sigma), ylab=coe.type, at=seq(0,1, 0.05), bty="n", colorkey=FALSE, bty="n", par.settings = list(axis.line = list(col = "transparent")), scales=scales.list, main=title, screen = list(z = -60, x = -60), )
      #scales=list(arrows=FALSE, tick.number=5, col="black")
    }
  } #next mut.i


  #layout(mat=matrix(nrow=4, ncol=4, data=c(seq(1,12), rep(13, 4)), byrow=TRUE), widths=c(0.5,rep(3, 3)), heights=c(0.25, rep(3,2), 0.75)) -> l
  #layout.show(l)
  #empty.plot()
  #empty.plot()
  #text(labels="Stickbreaking", x=0.5, y=0.5, font=2, cex=lab.cex)
  #empty.plot()
  #text(labels="Multiplicative", x=0.5, y=0.5, font=2, cex=lab.cex)
  #empty.plot()
  #text(labels="Additive", x=0.5, y=0.5, font=2, cex=lab.cex)


  wdth <- 14
  ht <- 14
  for (plot.i in 1:2){
    if (plot.i == 1){
      svg(file=paste(outpath, ".svg", sep=""), width=wdth, height=ht)
    } else if (plot.i == 2){
      pdf(file=paste(outpath, ".pdf", sep=""), width=wdth, height=ht)
    }

    ncol <- 3
    nrow <- 3
    for (mut.i in 1:length(posterior.list)){
      n.muts <- posterior.list[[mut.i]]$n.muts[1]
      #empty.plot()
      #text(label=paste(n.muts, " mutations", sep=""), x=0.5, y=0.5,srt=90, font=2, cex=lab.cex)
      for (mod.i in 1:3){
        print(plot.list[[mut.i]][[mod.i]], split=c(mut.i, mod.i, ncol, nrow), more=TRUE)
      }
    }  #next mut.i
    dev.off()
  } #next plot.i

  pdf(file=paste(outpath, "_scale.pdf", sep=""), width=wdth/2, height=ht/2)
    wireframe(d, shade=FALSE, light.source=c(0,0,10), col.regions=grey.cols, drape=TRUE, zlim=c(0, 1), row.values=sig.v, column.values=coe.v, col="black", zlab=list("Mean posterior", rot=90), xlab=expression(sigma), ylab=coe.type, scales=list(arrows=FALSE, tick.number=5, col="black"), at=seq(0,1, 0.05), bty="n", colorkey=TRUE, bty="n", par.settings = list(axis.line = list(col = "transparent")))
  dev.off()
}





#'  Create individual wireframe figures of posterior probabilities from analyzed simulated data
#'
#' @inheritParams create.posterior.wireframe.panel
#' @return Nothing. Creates pair of plots (pdf and svg) for each model/mutation
#' combination written to \code{outpath}
#' @export

create.posterior.wireframe.individual <- function(posterior.list, outpath, coes.to.plot=NA, sigs.to.plot=NA, plot.type="mean.post"){
  grey.cols <- c(paste("grey", seq(0,50), sep=""), paste("grey", seq(60, 100, by=5)))
  model.coes <- c(expression(u), expression(s), expression(paste(Delta, w)))
  n.muts <- length(posterior.list)

  model.v <- c("stick", "mult", "add")
  coe.col <- which(colnames(posterior.list[[1]])=="coes")
  sig.col <- which(colnames(posterior.list[[1]])=="sigma")
  if (is.na(coes.to.plot[1])){
    coe.v <- sort(unique(posterior.list[[1]]$coe))
  } else{
    coe.v <- coes.to.plot
  }
  if (is.na(sigs.to.plot[1])){
    sig.v <- sort(unique(posterior.list[[1]]$sigma))
  } else{
    sig.v <- sigs.to.plot
  }
  plot.list <- vector("list", length(posterior.list))
  for (i in 1:3){
    plot.list[[i]] <- vector("list",3)
  }

  scales.list <- vector("list", 4)
  names(scales.list) <- c("x", "y", "arrows", "col")
  scales.list$x$at <- sig.v
  scales.list$x$labels <- round(sig.v, 2)
  scales.list$y$at <- coe.v
  scales.list$y$labels <- coe.v
  scales.list$arrows <- FALSE
  scales.list$col <- "black"

  col.to.plot <- which(colnames(posterior.list[[1]])==plot.type)
  if (plot.type=="mean.post"){
    zlab <- "Mean posterior"
  } else if (plot.type == "true.unq"){
    zlab <- "Proportion"
  }

  for (mut.i in 1:length(posterior.list)){
    n.muts <- posterior.list[[mut.i]]$n.muts[1]
    for (mod.i in 1:3){
      model <- model.v[mod.i]
      subdata <- posterior.list[[mut.i]][which(posterior.list[[mut.i]]$model==model),]
      data.m <- matrix(nrow=length(coe.v), ncol=length(sig.v))
      for (coe.i in 1:length(coe.v)){
				for (sig.i in 1:length(sig.v)){
					the.row <- which(subdata$coes==coe.v[coe.i] & subdata$sigma==sig.v[sig.i])
					data.m[coe.i, sig.i] <- subdata[the.row, col.to.plot]
				}  #next sig.i
      }   #next coe.i
      d <- t(data.m)
      coe.type <- model.coes[mod.i]
      plot.list[[mut.i]][[mod.i]] <- wireframe(d, shade=FALSE, light.source=c(0,0,10), col.regions=grey.cols, drape=TRUE, zlim=c(0, 1), row.values=sig.v, column.values=coe.v, col="black", zlab=list(zlab, rot=90), xlab=expression(sigma), ylab=coe.type, at=seq(0,1, 0.05), bty="n", colorkey=FALSE, bty="n", par.settings = list(axis.line = list(col = "transparent")), scales=scales.list)
    }
  } #next mut.i



  wdth <- 5
  ht <- 5

  ncol <- 3
  nrow <- 3
  for (mut.i in 1:length(posterior.list)){
    n.muts <- posterior.list[[mut.i]]$n.muts[1]
    #empty.plot()
    #text(label=paste(n.muts, " mutations", sep=""), x=0.5, y=0.5,srt=90, font=2, cex=lab.cex)
    for (mod.i in 1:3){
      model <- model.v[mod.i]
      for (plot.i in 1:2){
        if (plot.i == 1){
          svg(file=paste(outpath,"_", model, "_", n.muts,"muts.svg", sep=""), width=wdth, height=ht)
        } else if (plot.i == 2){
          pdf(file=paste(outpath,"_", model, "_", n.muts,"muts.pdf", sep=""), width=wdth, height=ht)
        }
        print(plot.list[[mut.i]][[mod.i]])
        dev.off()
      } #next plot.i
    } #next mod.i
  }  #next mut.i


  pdf(file=paste(outpath, "_scale.pdf", sep=""), width=wdth/2, height=ht/2)
    wireframe(d, shade=FALSE, light.source=c(0,0,10), col.regions=grey.cols, drape=TRUE, zlim=c(0, 1), row.values=sig.v, column.values=coe.v, col="black", zlab=list("Mean posterior", rot=90), xlab=expression(sigma), ylab=coe.type, scales=list(arrows=FALSE, tick.number=5, col="black"), at=seq(0,1, 0.05), bty="n", colorkey=TRUE, bty="n", par.settings = list(axis.line = list(col = "transparent")))
  dev.off()
}



#' Create plot of observed vs model fit data
#'
#' @param outpath Path including filename to direct plot
#' @param fit.smry Summary of fit statistics from
#' @param preds Data frame with genotype rows (rownames are 0/1 strings), predicted fitness under each model
#' (column headers "stick", "mult" and "add") and observed fitness (column header "obs")
#' @param print.labels Print genotype labels. TRUE/FALSE. Default=TRUE.
#' @param pos.v Vector of which direction from to offest lables (pos argument in text function)
#' @param off.x Vector of offsets in x dimension
#' @param off.y Vector of offsets in y dimension
#' @param lims Vector of length 2 specifying plot limits (x and y get same limits)
#' @param cex.v Size of points for stickbreaking, multiplicative and additive models
#' @param col.v Color of points for stickbreaking, multiplicative and additive models
#' @param text.cex Size of the labels. Defaults=1.
#' @param xlab, ylab Labels for x and y axes. Default xlab="Observed Fitness", ylab="Model Predicted Fitness".
#' @return Nothing. Creates plot.
#' @export

create.model.fit.to.data.plot <- function(outpath, fit.smry, preds, print.labels=TRUE, pos.v, off.x, off.y, lims, cex.v=c(1,1,1), col.v, text.cex=0.75, xlab="Observed Fitness", ylab="Model Predicted Fitness"){
  for (plot.i in 1:2){
    if (plot.i ==1){
      file.out <- paste(outpath, ".svg", sep="")
      svg(file=file.out, width=5, height=5)
    } else if (plot.i ==2){
      file.out <- paste(outpath, ".pdf", sep="")
      pdf(file=file.out, width=5, height=5)
    }

    plot(x=preds$obs, y=preds$stick, ylim=lims, xlim=lims, ylab=ylab, xlab=xlab, pch=21, bg=col.v[1], cex=cex.v[1])
    abline(0,1, lty="dashed")
    points(x=preds$obs, y=preds$mult, pch=21, bg=col.v[2], cex=cex.v[2])
    points(x=preds$obs, y=preds$add, pch=21, bg=col.v[3], cex=cex.v[3])
    points(x=preds$obs, y=preds$stick, pch=21, bg=col.v[1], cex=cex.v[1])
    maxy <- apply(preds[,1:3], 1, function(x) max(x))
    miny <- apply(preds[,1:3], 1, function(x) min(x))
    minmax <- miny
    minmax[which(pos.v==3)] <- maxy[which(pos.v==3)]
    if (print.labels==TRUE){
      for (i in 1:length(preds$obs)){
        text(x=preds$obs[i]+off.x[i], y=minmax[i], labels=rownames(preds)[i], pos=pos.v[i], srt=90, offset=off.y[i], cex=text.cex)
      }
    }
    #legend("topleft", legend=c(as.expression(bquote("Stick (" ~ R^2==.(round(fit.smry$R2.stick,2)),")")), "Mult", "Add"), pch=21,pt.bg=col.v, bty="n")
    legend("topleft", pch=21, pt.bg=col.v, bty="n", legend=c(as.expression(bquote("Stick," ~ R^2==.(round(fit.smry$R2.stick,2)))), as.expression(bquote("Mult," ~ R^2==.(round(fit.smry$R2.mult,2)))), as.expression(bquote("Add," ~ R^2==.(round(fit.smry$R2.add,2))))))

    dev.off()
  }
}


#' Create 2 planel plot of observed vs model fit data
#'
#' @param outpath Path including filename to direct plot
#' @param smry.w.probs List with fit.smry statistics and posterior probabilties
#' @param preds Data frame with genotype rows (rownames are 0/1 strings), predicted fitness under each model
#' (column headers "stick", "mult" and "add") and observed fitness (column header "obs")
#' @param print.labels Print genotype labels. TRUE/FALSE. Default=TRUE.
#' @param pos.v Vector of which direction from to offest lables (pos argument in text function)
#' @param off.x Vector of offsets in x dimension
#' @param off.y Vector of offsets in y dimension
#' @param lims Vector of length 2 specifying plot limits (x and y get same limits)
#' @param cex.v Size of points for stickbreaking, multiplicative and additive models
#' @param col.v Color of points for stickbreaking, multiplicative and additive models
#' @param text.cex Size of the labels. Defaults=1.
#' @param xlab, ylab Labels for x and y axes. Default xlab="Observed Fitness", ylab="Model Predicted Fitness".
#' @param smry.w.probs.2 List with fit.smry statistics and posterior probabilties for second plot.
#' @param preds.2 Data frame for 2nd plot with genotype rows (rownames are 0/1 strings), predicted fitness under each model
#' (column headers "stick", "mult" and "add") and observed fitness (column header "obs")
#' @param print.labels.2 Print genotype labels in 2nd plot. TRUE/FALSE. Default=TRUE.
#' @param pos.v.2 Vector of which direction from to offest lables (pos argument in text function) in 2nd plot.
#' @param off.x.2 Vector of offsets in x dimension in 2nd plot.
#' @param off.y.2 Vector of offsets in y dimension in 2nd plot.
#' @param lims.2 Vector of length 2 specifying plot limits (x and y get same limits) in 2nd plot.
#' @param cex.v.2 Size of points for stickbreaking, multiplicative and additive models in 2nd plot.
#' @param text.cex.2 Size of the labels ijn 2nd plot. Defaults=1.
#' @param xlab.2, ylab.2 Labels for x and y axes in 2nd plot. Default xlab="Observed Fitness", ylab="Model Predicted Fitness".
#' @param leg.loc Location of each of the two legends. Default is "bottomright" for both.
#' @param leg.cex Text size in legend. Default = 1.
#' @param leg.bty bty argument for the legenes. Default = "n".
#' @param R2.digits Number of digits to round R2 values to when printing in legend. Default=3.
#' @param P.digits Number of digits to round P scores to when printing in legend. Default=2.
#' @return Nothing. Creates plot.
#' @export

create.2panel.model.fit.to.data.plot <- function(outpath, smry.w.probs, preds, print.labels=TRUE, pos.v, off.x, off.y, lims, cex.v=c(1,1,1), col.v, text.cex=0.75, xlab="Observed Fitness", ylab="Model Predicted Fitness", smry.w.probs.2, preds.2, print.labels.2=TRUE, pos.v.2, off.x.2, off.y.2, lims.2, cex.v.2=c(1,1,1), text.cex.2=0.75, xlab.2="Observed Fitness", ylab.2="Model Predicted Fitness", leg.loc=c("bottomright", "bottomright"), leg.cex=1, leg.bty="n", R2.digits=3, P.digits=2){
  for (plot.i in 1:2){
    if (plot.i ==1){
      file.out <- paste(outpath, ".svg", sep="")
      svg(file=file.out, width=10, height=5)
    } else if (plot.i ==2){
      file.out <- paste(outpath, ".pdf", sep="")
      pdf(file=file.out, width=10, height=5)
    }

    layout(mat=matrix(nrow=1, ncol=2, data=seq(1,2)), widths=c(5,5), heights=5)

    plot(x=preds$obs, y=preds$stick, ylim=lims, xlim=lims, ylab=ylab, xlab=xlab, pch=21, bg=col.v[1], cex=cex.v[1])
    abline(0,1, lty="dashed")
    points(x=preds$obs, y=preds$mult, pch=21, bg=col.v[2], cex=cex.v[2])
    points(x=preds$obs, y=preds$add, pch=21, bg=col.v[3], cex=cex.v[3])
    points(x=preds$obs, y=preds$stick, pch=21, bg=col.v[1], cex=cex.v[1])
    maxy <- apply(preds[,1:3], 1, function(x) max(x))
    miny <- apply(preds[,1:3], 1, function(x) min(x))
    minmax <- miny
    minmax[which(pos.v==3)] <- maxy[which(pos.v==3)]
    if (print.labels==TRUE){
      for (i in 1:length(preds$obs)){
        text(x=preds$obs[i]+off.x[i], y=minmax[i], labels=rownames(preds)[i], pos=pos.v[i], srt=90, offset=off.y[i], cex=text.cex)
      }
    }
    legend(leg.loc[1], pch=21, pt.bg=col.v, bty=leg.bty, legend=c(as.expression(bquote("Stick," ~ R^2==.(round(smry.w.probs$R2.stick,R2.digits)) ~ "," ~ P[post]==.(round(smry.w.probs$stick,2)))), as.expression(bquote("Mult," ~ R^2==.(round(smry.w.probs$R2.mult,2))~ "," ~ P[post]==.(round(smry.w.probs$mult,2)))), as.expression(bquote("Add," ~ R^2==.(round(smry.w.probs$R2.add,2))~ "," ~ P[post]==.(round(smry.w.probs$add,P.digits))))), cex=leg.cex)
    mtext(text="A", side=3, line=1, outer=FALSE, adj=0, cex=1.5, font=2)

    plot(x=preds.2$obs, y=preds.2$stick, ylim=lims.2, xlim=lims.2, ylab=ylab.2, xlab=xlab.2, pch=21, bg=col.v[1], cex=cex.v.2[1])
    abline(0,1, lty="dashed")
    points(x=preds.2$obs, y=preds.2$mult, pch=21, bg=col.v[2], cex=cex.v.2[2])
    points(x=preds.2$obs, y=preds.2$add, pch=21, bg=col.v[3], cex=cex.v.2[3])
    points(x=preds.2$obs, y=preds.2$stick, pch=21, bg=col.v[1], cex=cex.v.2[1])
    maxy <- apply(preds.2[,1:3], 1, function(x) max(x))
    miny <- apply(preds.2[,1:3], 1, function(x) min(x))
    minmax <- miny
    minmax[which(pos.v.2==3)] <- maxy[which(pos.v.2==3)]
    if (print.labels.2==TRUE){
      for (i in 1:length(preds.2$obs)){
        text(x=preds.2$obs[i]+off.x.2[i], y=minmax[i], labels=rownames(preds.2)[i], pos=pos.v.2[i], srt=90, offset=off.y.2[i], cex=text.cex.2)
      }
    }
    legend(leg.loc[2], pch=21, pt.bg=col.v, bty=leg.bty, legend=c(as.expression(bquote("Stick," ~ R^2==.(round(smry.w.probs.2$R2.stick,2)) ~ "," ~ P[post]==.(round(smry.w.probs.2$stick,2)))), as.expression(bquote("Mult," ~ R^2==.(round(smry.w.probs.2$R2.mult,2))~ "," ~ P[post]==.(round(smry.w.probs.2$mult,2)))), as.expression(bquote("Add," ~ R^2==.(round(smry.w.probs.2$R2.add,R2.digits))~ "," ~ P[post]==.(round(smry.w.probs.2$add,P.digits))))), cex=leg.cex)
    mtext(text="B", side=3, line=1, outer=FALSE, adj=0, cex=1.5, font=2)

    dev.off()
  }
}


