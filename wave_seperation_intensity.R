install.packages("signal")
install.packages("pracma")


fsg721 <- function(x) {
  # 1st derivative with SG filter
  #2nd order polynomial
  C = c(0.107143, 0.071429, 0.035714)
  B = integer(7)
  for (i in 1:3) {
    B[i] = C[i]
  }
  B[4] = 0.0
  for (i in 5:7) {
    B[i] = -C[8 - i]
  }
  A = c(1, 0)
  s = length(x)
  dx = signal::filter(B, A, x)
  dx = c(dx[7], dx[7], dx[7], dx[7:s], dx[s], dx[s], dx[s])
}

Tintersect <- function(wf, plot = FALSE) {
  
  k1 <- wf - min(wf[1:which.max(wf)])
  xvar <- (1:length(k1) - 1)
  spl <- smooth.spline(k1 ~ xvar)
  newx <- which.max(diff(k1))
  pred0 <- predict(spl, x = newx, deriv = 0)
  pred1 <- predict(spl, x = newx, deriv = 1)
  yint <- pred0$y - (pred1$y * newx)
  xint <- (-yint / pred1$y)
  
  if(isTRUE(plot)) {
    plot(xvar, k1, ylim=c(min(k1)-20, max(k1)))
    abline(h=min(k1), col="red", lty=3)
    lines(spl, col="red") 
    lines(xvar, yint + pred1$y*xvar, col="green", lwd=2)
    points(pred0,col="red", pch=8, lwd=2) 
    points(xint, 0, col="red", pch=8, lwd=2) 
    abline(v=xint, lty=3)
  }
  
  return(xint)
  
}

dicrotic <- function(pw, plot = FALSE) {
  
  # Get derivatives
  dp1 <- fsg721(pw)
  dp2 <- fsg721(fsg721(pw))
  dp3 <- fsg721(fsg721(fsg721(pw)))
  
  
  # FIND DICROTIC DEPRESSION ------------------------------------------------
  
  # Isolate notch area with 2nd and 3rd derivatives
  nni <- which.min(dp1)
  end <- length(pw)
  
  # End index without potential perturbation at end diastole  
  end2 <- end * .9
  
  # Dicrotic notch from local dp2 max
  dic <- which.max(dp2[nni:end2]) + nni - 1
  
  # plot(pw, type="l", lwd=2)
  # par(new=T)
  # plot(dp2, type='o',col="grey")
  # abline(v = dic, h = 0)
  
  
  # FIND DICROTIC PEAK ------------------------------------------------------
  
  end3 <- ((end - dic) * .6) + dic # 60% of diastolic duration
  #abline(v=end3, lty=2, col=2)
  
  if(sum(dp2[dic:end3] < 0) < 1) {
    dia <- 9999
  } else {
    dia <- which.min(dp2[dic:end3]) + dic - 1
  }
  
  # plot(pw, type="l", lwd=2)
  # par(new=T)
  # plot(dp2, type='o',col="grey")
  # abline(v = c(dic, dia), h = 0)
  
  
  # PLOTS -------------------------------------------------------------------
  
  if(isTRUE(plot)) {
    plot(pw, type = "l", lwd=2, ylab="BP (mmHg)")
    abline(v=c(dic, dia), col="grey", lty=3, lwd=2)
    mtext(c("Ed", "P3"), side = 3, at = c(dic,dia))
  }
  
  return(data.frame(dicrotic_notch = dic, 
                    dicrotic_peak = dia))
  
}

low.pass <- function(y, fq, do.plot = FALSE) {
  
  # Second order low pass filter
  # Removes high frequency components above fq
  # y = a numeric vector
  # fq = a numeric vector giving frequency or period of the filter.
  
  if (any(is.na(y))) stop("y contains NA")
  
  # n = a numeric value giving the order of the filter. 
  # Larger numbers create steeper fall off.
  n = 4
  
  if (any(fq>1)) {
    f <- 1/fq
    p <- fq
  } else {
    p <- 1/fq
    f <- fq
  }
  
  # sort f in case it's passed in backwards
  f <- sort(f)
  
  filt <- signal::butter(n = n,
                         W = f * 2,
                         type = "low",
                         plane = "z")
  
  # remove mean
  yAvg <- mean(y)
  y <- y - yAvg
  
  # pad the data to twice the max period
  pad <- max(p) * 2
  ny <- length(y)
  
  # pad the data
  yPad <- c(y[pad:1], y, y[ny:(ny - pad)])
  
  # run the filter
  yFilt <- signal::filtfilt(filt, yPad)
  
  # unpad the filtered data
  yFilt <- yFilt[(pad + 1):(ny + pad)]
  
  # return with mean added back in
  filt.sig <- yFilt + yAvg
  
  if(isTRUE(do.plot)){
    # plot results
    plot(filt.sig,
         type = "l",
         lwd = 2)
  }
  
  # return filtered signal
  return(filt.sig)
}

RootSpline1 <- function (x, y, y0 = 0, verbose = TRUE) {
  
  if (is.unsorted(x)) {
    ind <- order(x)
    x <- x[ind]; y <- y[ind]
  }
  z <- y - y0
  ## which piecewise linear segment crosses zero?
  k <- which(z[-1] * z[-length(z)] <= 0)
  ## analytical root finding
  xr <- x[k] - z[k] * (x[k + 1] - x[k]) / (z[k + 1] - z[k])
  ## make a plot?
  if (verbose) {
    plot(x, y, "l"); abline(h = y0, lty = 2)
    points(xr, rep.int(y0, length(xr)))
  }
  ## return roots
  xr
}

pwa <- function(pw, filt = FALSE, plot = FALSE) {
  
  # Low pass waveform
  if (isTRUE(filt)) {
    pw <- low.pass(pw, 10, do.plot = F)
  }
  
  # Create derivatives
  d1 <- fsg721(pw)
  d2 <- fsg721(fsg721(pw))
  d3 <- fsg721(fsg721(fsg721(pw)))
  d4 <- fsg721(fsg721(fsg721(fsg721(pw))))
  
  # Some additional calcs
  end <- length(pw)
  foot <- Tintersect(pw) - 1
  if(foot < 1) {foot <- 1}
  maxpi <- which.max(pw)
  notchdat <- dicrotic(pw)
  notch <- notchdat$dicrotic_notch
  notchpeak <- notchdat$dicrotic_peak
  
  # plot(pw,type="o")
  # abline(v=c(foot, notch, notchpeak))
  
  # Create time
  time <- (0:(length(pw)-1)) / 200
  X <- 1:length(pw)
  
  # get zero crossing of 4th derivative
  zero_cross <- RootSpline1(X[round(foot):end], d4[round(foot):end], 
                            verbose = F)
  
  # index of p1
  p1i <- zero_cross[2]
  
  # plot
  # plot(pw[foot:end], type='l', lwd=2, col='grey',yaxt='n',ylab="")
  # par(new=T)
  # plot(d4,type="l",xaxt='n',lwd=2,ylab="4th derivative")
  # abline(h=0,v=p1i,col=2,lwd=1.5)
  
  # Find p2 from 3rd derivative
  # when looking for p2 after sbp, the 3rd derivative method is more robust than
  # the 4th derivative method.
  p2i <- which.min(d3[maxpi:(notch - 5)]) + maxpi
  
  
  # Depending type of pressure waveform p1 or p2 will approx equal max p
  # Find which is closest to max P
  distp1 <- abs(maxpi - p1i)
  distp2 <- abs(maxpi - p2i)
  
  # Which is closer to max P, p1 | p2
  if(distp1 > distp2) {
    p2i <- which.max(pw)
  } else if(distp2 > distp1) {
    p1i <- which.max(pw)
  }
  
  # Calculate augmentation index
  ap <- (pw[p2i] - pw[p1i])
  pp <- (pw[maxpi] - pw[foot])
  aix <- (ap / pp) * 100
  
  # Determine Murgo waveform type
  type <- NA
  if(aix > 12) {
    type <- "A"
  } else  if (aix <= 12 & aix >= 0) {
    type <- "B"
  } else if (aix < 0) {
    type <- "C"
  }
  
  # plot(pw,type="o")
  # abline(v=c(p1i,p2i), col=2:3)
  
  # Find dp/dt max
  dpdt.max <- which.max(d1)
  
  # Find inflection after p1 (local max of 2nd derivative per SphygmoCor)
  p1i2 <- 9999
  if(p1i < maxpi-3) {
    p1i2 <- which.max(d2[p1i:maxpi]) + (p1i - 1)
  }
  
  # Find p1 from 1st derivative (local min as per Kelly et al. 10.1161/01.CIR.80.6.1652)
  # This method was simplified to the 4th derivative method but works fine here
  # p1i1 <- 9999
  # if(p1i < maxpi-3) {
  #   p1i1 <- which.min(d1[dpdt.max:p1i2]) + (dpdt.max - 1)
  # }
  
  # plot(pw,type="o"); abline(v=c(p1i, p1i1, p1i2), col=2)
  # par(new=T)
  # plot(d1,type="o"); abline(h=0)
  
  # Plot results
  if (isTRUE(plot)) {
    
    plot(time, pw,
         type = 'l',
         lwd = 3,
         ylab = "Pressure (mmHg)",
         xlab = "Time (s)")
    grid(NULL,NULL, lty = 3, col = "lightgrey") 
    legend("topright", type, bty = 'n')
    
    points(x = c(time[foot], 
                 time[p1i], 
                 time[p1i2], 
                 time[p2i], 
                 time[notch], 
                 time[notchpeak]),
           y = c(pw[foot], 
                 pw[p1i], 
                 pw[p1i2], 
                 pw[p2i], 
                 pw[notch], 
                 pw[notchpeak]),
           pch = "|",
           col = 2,
           lwd = 3,
           cex = 1.7)
    
  }
  
  
  df <- data.frame(
    # Index of values
    MaxP_index = maxpi,
    Foot_index = foot,
    P1_index = p1i,
    P2_index = p2i,
    Ed_index = notch,
    P3_index = notchpeak,
    #P1x_index = p1i1,
    DpDt_index = dpdt.max,
    # Values in unit seconds
    MaxP_sec = time[maxpi],
    Foot_sec = time[foot],
    P1_sec = time[p1i],
    P2_sec = time[p2i],
    Ed_sec = time[notch],
    P3_sec = time[notchpeak],
    #P1x_sec = time[p1i1],
    DpDt_sec = time[dpdt.max],
    # Values in unit mmHg
    MaxP_mmhg = pw[maxpi],
    Foot_mmhg = pw[foot],
    P1_mmhg = pw[p1i],
    P2_mmhg = pw[p2i],
    Ed_mmhg = pw[notch],
    P3_mmhg = pw[notchpeak],
    #P1x_mmhg = pw[p1i1],
    DpDt_mmhg = pw[dpdt.max],
    AP_mmHg = ap,
    AIX = aix,
    Type = type
  )
  
  # round values in df
  num_cols <- unlist(lapply(df, is.numeric)) # Identify numeric cols
  df[num_cols] <-  round(df[num_cols], 3)    # round numeric cols
  
  # # print values to console
  # for(i in 1:length(df)){
  #   print(paste0(names(df[i]),": ", df[1,i]), quote = F)
  # }
  
  return(df)
  
}

MattLag <- function(x, k) {
  
  if (k > 0) {
    d1 <- tail(x, k)
    d2 <- c(rep(NA, k), x)[1:length(x)]
    d2[is.na(d2)] <- d1
    return(d2)
  } else {
    b1 <- head(x, abs(k))
    b2 <- c(x[(-k + 1):length(x)], rep(NA, -k))
    b2[is.na(b2)] <- b1
    return(b2)
  }
}


wswi <- function(pressure, 
                 flow,
                 align = FALSE, 
                 plot = FALSE) {
  
  # pressure = averaged pressure waveform (numeric)
  # flow = averaged flow waveform (numeric)
  # lowpass = apply 10Hz low pass to pressure waveform (T/F)
  # align = time align pressure and flow waveforms (T/F)
  # plot = plot results (T/F)
  
  # Set up ------------------------------------------------------------------
  
  # Load data
  pw <- pressure[!is.na(pressure)]
  fw <- flow[!is.na(flow)]
  
  # Check p & u signals the same length
  if (length(pw) > length(fw)) {
    pw = pw[1:length(fw)]
  } else if (length(pw) < length(fw)) {
    fw = fw[1:length(pw)]
  }
  
  # Remove diastole
  pensavg <- pw - min(pw[1:which.max(pw)])
  uensavg <- fw - min(fw[1:which.max(fw)])
  
  sr <- 200 # sample rate
  
  # Make time variable
  tens <- (0:(length(pensavg) - 1)) / sr
  
  # Get waveform indices ----------------------------------------------------
  
  # Get PWA parameters
  params <- pwa(pw)
  
  Ifoot <- params$Foot_index
  Inotch <- params$Ed_index 
  Idiapeak <- params$P3_index
  Ip1 <- params$P1_index
  Ip2 <- params$P2_index
  AugInx <- params$AIX
  

  # Align PQ ----------------------------------------------------------------
  
  # Find foot with intersecting tangents
  xint.P <- Tintersect(pensavg)
  xint.U <- Tintersect(uensavg)
  
  if (isTRUE(align)) {
    maxlag <- round(xint.P - xint.U)
    
    if (maxlag != 0) {
      uensavg <-
        MattLag(uensavg, maxlag) # Lag flow so it aligns with pressure
    }
  }
  
  
  # Estimate Zc -------------------------------------------------------------

  # Isolate upstroke
  Imq <- which.max(uensavg)
  
  # Max of linear fit
  lmmax <- round(0.5*(Imq-round(Ifoot)))
  
  # The slope of P ~ U in early systole = blood density (rho) * compliance (c)
  linea <- 
    lm(pensavg[round(Ifoot):(lmmax+round(Ifoot))] ~ uensavg[round(Ifoot):(lmmax+round(Ifoot))])
  
  rho <- 1050                      # Blood density
  rhoc <- linea$coefficients[[2]]  # Slope <- Zc
  c <- rhoc / rho                  # c = wave speed
  

  # Wave intensity ----------------------------------------------------------

  # Matt method
  # p = filter order
  # n = filter length (must be odd)
  # m = return mth derivative
  # ts = time scaling factor
  dpe <- signal::sgolayfilt(
    pensavg,
    p = 2,
    n = 11,
    m = 1,
    ts = 1)
  
  due <- signal::sgolayfilt(
    uensavg,
    p = 2,
    n = 11,
    m = 1,
    ts = 1)
  
  die <- dpe * due
  
  # plot(die, type="l", lwd=2)
  
  # Rhoc.ss is rhoc but calculated via the sum of squares method
  rhocSS <- sqrt(sum(dpe[!is.na(dpe)] ^ 2, rm.na = TRUE) / sum(due[!is.na(due)] ^ 2, rm.na = TRUE))
  cSS <- rhocSS / rho
  
  # Calculate forward & backward WI
  dpep <- (dpe + rhoc * due) / 2 # Forward pressure difference eg. (P+Zc*U)/2
  dpem <- (dpe - rhoc * due) / 2 # Backward pressure difference eg. (P-Zc*U)/2
  duem <- -dpem / rhoc           # Backward velocity difference
  duep <- dpep / rhoc            # Forward velocity difference
  diep <- dpep * duep            # Forward wave intensity
  diem <- dpem * duem            # Backward wave intensity
  
  
  # Save variables for finding peaks and nadirs of WI
  TImaxp <- which.max(pensavg)                 # index of SBP
  lsys <- which.min(dpe[TImaxp:length(dpe)])   # min between max(p):end
  lsys <- round(lsys) + 4                      # round + 5 samples for margin of error
  lsys <- lsys + TImaxp
  
  
  # Find max WI as per Alun
  w1 <- pracma::findpeaks(die[1:lsys], npeaks = 1, minpeakheight = 0.7 *
                            max(die))
  w1val <- w1[1]         # WI @ w1 (peak)
  w1loc <- w1[2]         # index of w1
  
  # Find min in systole
  w0 <- pracma::findpeaks(-die[w1loc:lsys],
                          npeaks = 1,
                          minpeakheight = 0.7 * max(-die[w1loc:lsys]))
  
  if (is.null(w0)) {
    isw0 = TRUE
    #print(isw0)
    w0loca = which.min(die[w1loc:lsys])
    w0loc = w0loca + w1loc - 1
    w0val = die[w0loc]
  } else {
    isw0 = FALSE
    #print(isw0)
    w0val = w0[1]
    w0loca = w0[2]
    w0loc = w1loc + w0loca - 1
  }
  
  # Find second peak >5% of max WI
  minpeak <- max(die) / 20    # so peaks <5% of w1 ignored for w2
  w2 <- pracma::findpeaks(die[Inotch:1], npeaks = 1, minpeakheight = minpeak)
  
  if (is.null(w2)) {
    isw2 = TRUE
    w2val = NA    # WI @ w2
    w2loc = NA    # index of w2
  } else {
    isw2 = FALSE
    w2val = w2[1] # WI @ w2
    w2loc = Inotch - w2[2] + 1 # index of w2
  }
  
  # area is 1.06447 * height * width
  w1area <- 1.06447 * w1val * length(die[w1[3]:w1[4]][die[w1[3]:w1[4]] > 0])
  
  #w0 area
  w0srt <- w0[3] + w1loc
  w0end <- w0[4] + w1loc
  if (isTRUE(isw0)) {
    w0area <- 0
  } else {
    w0area <- 1.06447 * w0val * length(die[w0srt:w0end][die[w0srt:w0end] < 0])
  }
  
  #w2 area
  w2srt <- Inotch - w2[3] + 1
  w2end <- Inotch - w2[4] + 1
  if (isTRUE(isw2)) {
    w2area <- 0
  } else {
    w2area <- 1.06447 * w2val * length(die[w2srt:w2end][die[w2srt:w2end] > 0])
  }
  
  wri <- w0area / w1area
  
  
  # Wave separation ---------------------------------------------------------
  
  pep <- (pensavg + uensavg * rhoc) / 2 # Pf
  pepmax <- max(pep)                    # amplitude
  pept <- tens[which.max(pep)]          # time
  
  pem <- (pensavg - uensavg * rhoc) / 2 # Pb
  pemmax <- max(pem)                    # amplitude
  pemt <- tens[which.max(pem)]          # time
  
  
  
  # PLOTS -------------------------------------------------------------------
  
  if(isTRUE(plot)) {
    
    par(mfrow = c(2, 2),
        mar = c(3.5, 3.5, .5, .5),
        mgp = c(2, 1, 0))
    
    # Plot Zc
    plot(uensavg,
         pensavg,
         type = "o",
         ylab = "Pressure",
         xlab = "Flow",
         pch = 21,
         bg = "cornflowerblue",
         col = 'grey40',
         lwd = 0.9,
         ylim = c(min(pensavg)-5, max(pensavg)+5),
         xlim = c(min(uensavg)-10, max(uensavg)+10),
         cex = 1.2,
         frame.plot = F)
    grid()
    clip(uensavg[round(Ifoot)],
         uensavg[lmmax + round(Ifoot)],
         -1000,
         1000)
    abline(linea, 
           col = "firebrick", 
           lwd = 3)
    
    
    # Plot P and U
    plot(tens,
         pensavg,
         type = "l",
         lwd = 2,
         ylab = "BP & Flow",
         xlab = "Time (s)",
         frame.plot = F)
    grid()
    par(new = T)
    plot(
      tens,
      uensavg,
      type = "l",
      lwd = 2,
      col = "firebrick",
      xaxt='n', 
      yaxt='n', 
      ann=FALSE,
      #ylim = c(-0.1, max(P + 5)),
      frame.plot = F)
    abline(v = tens[Inotch], lty = 3, lwd = 2)
    
    
    # Plot Pf & Pb
    plot(tens,
         pensavg,
         type = "l",
         lwd = 0,
         col = "grey20",
         lty = 1,
         xlab = "Time (s)",
         ylab = "Pressure (mmHg)",
         xlim = c(0, round(max(tens), 1)),
         ylim = c(trunc(min(pensavg))-5, ceiling(max(pensavg))+5),
         frame.plot = F,
         las = 2,
         axes = F); grid()
    axis(1, at=seq(0, round(max(tens), 1), 0.1),
         cex.axis=.8)
    axis(2, at = seq( trunc(min(pensavg))-5, ceiling(max(pensavg))+5, 5),
         cex.axis=.8)
    lines(tens, pensavg, lwd = 2)
    polygon(tens,
            pep,
            col = adjustcolor("orange", alpha.f = 0.5), 
            border = F)
    polygon(tens,
            pem,
            col = adjustcolor("cornflowerblue", alpha.f = 0.7), 
            border = F)
    lines(tens, pep,
          lwd=1.5,
          col="grey30")
    lines(tens, pem,
          lwd=1.5,
          col="grey30")
    
    # Plot WI
    plot(tens,
         die,
         type = "l",
         lwd = 0,
         xlim = c(0, round(max(tens), 1)),
         ylab = "",
         xlab = "",
         frame.plot = F,
         las = 2,
         axes = F)
    
    polygon(tens,
            die,
            col = adjustcolor("cornflowerblue", alpha.f = 0.9), 
            border = F)
    
    #abline(h=0)
    
    points(x = c(tens[w1loc],
                 tens[w0loc],
                 tens[w2loc]),
           y = c(die[w1loc],
                 die[w0loc],
                 die[w2loc]),
           pch = c(25,24,25),
           bg = "firebrick")
    
    par(new=T)
    
    plot(tens,
         pw,
         type = "l",
         lwd = 0,
         col = "grey20",
         lty = 1,
         xlab = "Time (s)",
         ylab = "Pressure (mmHg)",
         xlim = c(0, round(max(tens), 1)),
         ylim = c(trunc(min(pw)), ceiling(max(pw))),
         frame.plot = F,
         las = 2,
         axes = F); grid()
    
    axis(1, at=seq(0, round(max(tens), 1), 0.1),
         cex.axis=.8)
    
    axis(2, at= seq( trunc(min(pw)), ceiling(max(pw)), 5),
         cex.axis=.8)
    
    lines(tens, pw, lwd = 2)
    
  }
  
  
  # Results -----------------------------------------------------------------
  
  # Save raw waveform
  df1 <-
    data.frame(
      time = tens,
      p = pensavg,
      u = uensavg,
      wi = die,
      wiforward = diep,
      wibackward = diem,
      pf = pep,
      pb = pem)
  
  # Save parameters
  df2 <-
    data.frame(
      sbp = pw[TImaxp],
      Tsbp = tens[TImaxp],
      dbp = pw[Ifoot],
      Tdbp = tens[Ifoot],
      p1 = pw[Ip1],
      Tp1 = tens[Ip1],
      p2 = pw[Ip2],
      Tp2 = tens[Ip2],
      ed = pw[Inotch],
      Ted = tens[Inotch],
      p3 = pw[Idiapeak],
      Tp3 = tens[Idiapeak],
      map = mean(pw),
      Aix = AugInx,
      qmax = fw[Imq],
      Tqmax = tens[Imq],
      qfoot = fw[Ifoot],
      qmean = mean(fw),
      # Wave intensity
      wi1 = w1val,
      Twi1 = tens[w1loc],
      Awi1 = w1area,
      wi0 = w0val,
      Twi0 = tens[w0loc],
      Awi0 = w0area,
      wi2 = w2val,
      Twi2 = tens[w2loc],
      Awi2 = w2area,
      wrm = wri,
      zc = rhoc,
      zcss = rhocSS,
      comp = c,
      compss = cSS,
      # Wave separation
      pfmax = pepmax,
      Tpfmax = pept,
      Apf = sum(pep) / sr,
      pbmax = pemmax,
      Tpbmax = pemt,
      Apb = sum(pem) / sr,
      rm = pemmax / pepmax,
      ri = pemmax / (pemmax + pepmax)
    )
  
  # round values in df
  num_cols <- unlist(lapply(df2, is.numeric)) # Identify numeric cols
  df2[num_cols] <-  round(df2[num_cols], 3)   # round numeric cols
  
  # print values to console
  for(i in 1:length(df2)){
    print(paste0(names(df2[i]),": ", df2[1,i]), quote = F)
  }
  
  return(list(df1, df2))
  
} 
