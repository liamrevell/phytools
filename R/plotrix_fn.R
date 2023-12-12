## functions from plotrix version 3.8-3 (since package was orphaned on CRAN)
## citation: J L (2006). "Plotrix: a package in the red light district of R." R-News, 6(4), 8-12.
## authors: Jim Lemon, Ben Bolker, Sander Oom, Eduardo Klein, Barry Rowlingson, Hadley Wickham, 
## Anupam Tyagi, Olivier Eterradossi, Gabor Grothendieck, Michael Toews, John Kane, Rolf Turner, 
## Carl Witthoft, Julian Stander, Thomas Petzoldt, Remko Duursma, Elisa Biancotto, Ofir Levy, 
## Christophe Dutang, Peter Solymos, Robby Engelmann, Michael Hecker, Felix Steinbeck, 
## Hans Borchers, Henrik Singmann, Ted Toal, Derek Ogle, Darshan Baral, Ulrike Groemping, 
## Bill Venables, The CRAN Team
## The following code is modified code from the plotrix R package version 3.8-3, which is licensed 
## GPLv3. This code therefore is also licensed under the terms of the GNU Public License, version 3.

## First attempt to load functions from plotrix package (if installed). If not installed, functions 
## are loaded from this source file.

GetYmult<-function() {
 if(dev.cur() == 1) {
  warning("No graphics device open.")
  ymult<-1
 }
 else {
  # get the plot aspect ratio
  xyasp<-par("pin")
  # get the plot coordinate ratio
  xycr<-diff(par("usr"))[c(1,3)]
  ymult<-xyasp[1]/xyasp[2]*xycr[2]/xycr[1]
 }
 return(ymult)
}

Arctext<-function(x,center=c(0,0),radius=1,start=NULL,middle=pi/2,end=NULL,
 stretch=1,clockwise=TRUE,cex=NULL, ...) {

 oldcex <- par("cex")
 # have to do this to get strwidth to work
 if(is.null(cex)) cex <- oldcex
 par(cex = cex)
 xvec <- strsplit(x, "")[[1]]
 lenx <- length(xvec)
 xwidths <- stretch * strwidth(xvec)
 charangles <- xwidths/radius
 # make really narrow characters wider
 changrang <- range(charangles)
 charangles[charangles < changrang[2]/2] <- changrang[2]/2
 if(!is.null(end)) {
  if(clockwise) start <- end + sum(charangles)
  else start <- end - sum(charangles)
 }
 if(is.null(start)) {
  if (clockwise) start <- middle + sum(charangles)/2
  else start <- middle - sum(charangles)/2
 }
 if(clockwise) {
  charstart <- c(start, start - cumsum(charangles)[-lenx])
  charpos <- charstart - charangles/2
 }
 else {
  charstart <- c(start, start + cumsum(charangles)[-lenx])
  charpos <- charstart + charangles/2
 }
 xylim <- par("usr")
 plotdim <- par("pin")
 ymult <- (xylim[4] - xylim[3])/(xylim[2] - xylim[1]) * plotdim[1]/plotdim[2]
 for(xchar in 1:lenx) {
  srt <- 180 * charpos[xchar]/pi - 90
  text(center[1] + radius * cos(charpos[xchar]), center[2] + 
  radius * sin(charpos[xchar]) * ymult, xvec[xchar], 
   adj = c(0.5, 0.5), srt = srt + 180 * (!clockwise),...)
 }
 par(cex = oldcex)
}

Draw.arc <- function(x=1, y=NULL, radius=1, angle1=deg1*pi/180,
 angle2=deg2*pi/180, deg1=0, deg2=45, n=0.05, col=NA, lwd=NA, ...) {

    if (all(is.na(col)))
        col <- par("col")
    if (all(is.na(lwd)))
        lwd <- par("lwd")
    xylim<-par("usr")
    ymult <- getYmult()
    devunits <- dev.size("px")
    draw.arc.0 <- function(x, y, radius, angle1, angle2, n, col, lwd, ...)
        {
        delta.angle <- (angle2 - angle1)
        if (n != as.integer(n))
            n <- as.integer(1+delta.angle/n) # Divide total angle by desired segment angle to get number of segments
        delta.angle <- delta.angle/n
        angleS <- angle1 + seq(0, length=n) * delta.angle
        angleE <- c(angleS[-1], angle2)
        # Move segment starts/ends so that segments overlap enough to make wide segments
        # not have an open slice in them.  The slice is open by delta.angle*half.lwd.user.
        # That subtends an angle of that/(radius+half.lwd.user) radians, from center.
        # Move segment endpoints by half of that, so together they equal that.
        if (n > 1)
            {
            half.lwd.user <- (lwd/2)*(xylim[2]-xylim[1])/devunits[1]
            adj.angle = delta.angle*half.lwd.user/(2*(radius+half.lwd.user))
            angleS[2:n] = angleS[2:n] - adj.angle
            angleE[1:(n-1)] = angleE[1:(n-1)] + adj.angle
            }
        p1x <- x + radius * cos(angleS)
        p1y <- y + radius * sin(angleS) * ymult
        p2x <- x + radius * cos(angleE)
        p2y <- y + radius * sin(angleE) * ymult
        segments(p1x, p1y, p2x, p2y, col=col, lwd=lwd, ...)
        }
    xy <- xy.coords(x, y)
    x <- xy$x
    y <- xy$y
    a1 <- pmin(angle1, angle2)
    a2 <- pmax(angle1, angle2)
    angle1 <- a1
    angle2 <- a2
    args <- data.frame(x, y, radius, angle1, angle2, n, col, lwd, stringsAsFactors=FALSE)
    for (i in 1:nrow(args))
        do.call("draw.arc.0", c(args[i, ], ...))
    invisible(args)
    }

Draw.circle<-function(x,y,radius,nv=100,border=NULL,col=NA,
 lty=1,density=NULL,angle=45,lwd = 1) {

 xylim<-par("usr")
 plotdim<-par("pin")
 ymult<-getYmult()
 angle.inc<-2*pi/nv
 angles<-seq(0,2*pi-angle.inc,by=angle.inc)
 if(length(col)<length(radius)) 
  col<-rep(col,length.out=length(radius))
 for(circle in 1:length(radius)) {
  xv<-cos(angles)*radius[circle]+x
  yv<-sin(angles)*radius[circle]*ymult+y
  polygon(xv,yv,border=border,col=col[circle],lty=lty,
   density=density,angle=angle,lwd=lwd)
 }
 invisible(list(x=xv,y=yv))
}

Draw.ellipse<-function(x,y,a=1,b=1,angle=0,segment=NULL,arc.only=TRUE, 
 deg=TRUE,nv=100,border=NULL,col=NA,lty=1,lwd=1,...) {

    if(is.null(segment)) {
     # set segment to full ellipse if not supplied
     if(deg) segment<-c(0,360)
     else segment<-c(0,2*pi)
    }
    ## workhorse internal function draw ellipse
    draw1ellipse <-
    function(x, y, a = 1, b = 1, angle = 0, segment=NULL, 
    arc.only=TRUE, nv = 100, deg = TRUE, border=NULL, col=NA, lty=1, lwd=1, ...)
    {
        # if input is in degrees
        if (deg) {
            angle <- angle * pi/180
            segment <- segment * pi/180
        }
        z <- seq(segment[1], segment[2], length = nv + 1)
        xx <- a * cos(z)
        yy <- b * sin(z)
        alpha <- xyangle(xx, yy, directed = TRUE, deg = FALSE)
        rad <- sqrt(xx^2 + yy^2)
        xp <- rad * cos(alpha + angle) + x
        yp <- rad * sin(alpha + angle) + y
        if (!arc.only) {
            xp <- c(x, xp, x)
            yp <- c(y, yp, y)
        }
        polygon(xp, yp, border=border, col=col, lty=lty, lwd=lwd, ...)
        invisible(NULL)
    }
    ## internal function for the internal function
    xyangle <-
    function(x, y, directed = FALSE, deg = TRUE)
    {
        if (missing(y)) {
            y <- x[,2]
            x <- x[,1]
        }
        out <- atan2(y, x)
        if (!directed)
            out <- out %% pi   
        if (deg) # if output is desired in degrees
            out <- out * 180 / pi
        out
    }
    if (missing(y)) {
        y <- x[,2]
        x <- x[,1]
    }
    n <- length(x)
    if (length(a) < n)
        a <- rep(a, n)[1:n]
    if (length(b) < n)
        b <- rep(b, n)[1:n]
    if (length(angle) < n)
        angle <- rep(angle, n)[1:n]
    if (length(col) < n)
        col <- rep(col, n)[1:n]
    if (length(border) < n)
        border <- rep(border, n)[1:n]
    if (length(nv) < n)
        nv <- rep(nv, n)[1:n]
    if(n==1)
     draw1ellipse(x,y,a,b,angle=angle,segment=segment,
      arc.only=arc.only,deg=deg,nv=nv,col=col,border=border,
      lty=lty,lwd=lwd,...)
    else {
     if (length(segment) < 2*n)
        segment <- matrix(rep(segment,n), n, 2, byrow=TRUE)
     lapply(1:n, function(i) draw1ellipse(x[i], y[i], a[i], b[i], 
      angle=angle[i], segment=segment[i,], arc.only=arc.only, deg=deg, 
      nv=nv[i], col=col[i], border=border[i],
      lty=lty, lwd=lwd, ...))
    }
    invisible(NULL)
}

Textbox<-function(x,y,textlist,justify=c("l","c","r"),cex=1,leading=0.5,
 box=TRUE,adj=c(0,0),font=NULL,vfont=NULL,col=NULL,border=NULL,fill=NA,
 density=NULL,angle=45,lty=par("lty"),lwd=par("lwd"),margin=0) {

 if (length(margin) == 1) margin<-rep(margin, 4)
 else if (length(margin) == 2) margin<-rep(margin, 2)
 saveAdj<-par(adj = 0)
 textstr<-paste(textlist, collapse=" ")
 words<-strsplit(textstr," ")[[1]]
 line.height<-strheight("hy",cex=cex,font=font,vfont=vfont) * (1+leading)
 if(margin[2] > 0) x[1]<-x[1] + margin[2]
 if(margin[3] > 0) y<-y-margin[3]
 if(margin[4] > 0) x[2]<-x[2]-margin[4]
 if(x[1] >= x[2]) x[2]<-x[1]+diff(par("usr")[1:2])*0.1
 x.len<-diff(x)
 y.pos<-y
 x.pos<-x[1]
 adj2<-c(0,1)
 if(justify[1] == "c") {
  x.pos<-x.pos + x.len/2
  adj2[1]<-0.5
 }
 else {
  if(justify[1] == "r") {
   x.pos<-x.pos + x.len
   adj2[1]<-1
  }
 }
 curword<-1
 curline<-1
 txtline<-""
 while(curword <= length(words)) {
  txtline[curline]<- ""
  txtline[curline]<-paste(txtline[curline],words[curword])
  curword<-curword + 1
  while(strwidth(paste(txtline[curline],words[curword]),cex=cex, 
   font=font,vfont=vfont) < x.len && !is.na(words[curword])) {
   txtline[curline]<-paste(txtline[curline],words[curword])
   curword<-curword+1
  }
  curline<-curline+1
  y.pos[curline]<-y.pos[curline-1]-line.height
 }
 if(box) {
  xbox<-x
  ybox<-c(y.pos[curline],y)
  ybox[1]<-ybox[1]-abs(margin[1])
  xbox[1]<-xbox[1]-abs(margin[2])
  ybox[2]<-ybox[2]+abs(margin[3])
  xbox[2]<-xbox[2]+abs(margin[4])
  rect(xbox[1],ybox[1],xbox[2],ybox[2],border=border, 
   col=fill,density=density,angle=angle,lty=lty,lwd=lwd)
 }
 text(x.pos,y.pos[1:curline-1],txtline,adj=adj+adj2,cex=cex,
  col=col,font=font,vfont=vfont)
 par(saveAdj)
 return(y.pos)
}

arctext<-if(.check.pkg("plotrix")) plotrix::arctext else Arctext
draw.arc<-if(.check.pkg("plotrix")) plotrix::draw.arc else Draw.arc
draw.circle<-if(.check.pkg("plotrix")) plotrix::draw.circle else Draw.circle
draw.ellipse<-if(.check.pkg("plotrix")) plotrix::draw.ellipse else Draw.ellipse
getYmult<-if(.check.pkg("plotrix")) plotrix::getYmult else GetYmult
textbox<-if(.check.pkg("plotrix")) plotrix::textbox else Textbox
