WaveIntensity4 <- function(pressure, flow, align = F) {
    
  #  Wave Intensity 4.0
  #  Matthew K. Armstrong (matthew-k-armstrong@uiowa.edu)
  #  GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007
  #  Copyright 2021 Matthew K. Armstrong (matthew.armstrong@utas.edu.au)
  #  Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
  #  Everyone is permitted to copy and distribute verbatim copies
  #  of this license document, but changing it is not allowed.
  #
  #  This program is distributed in the hope that it will be useful,
  #  but WITHOUT ANY WARRANTY; without even the implied warranty of
  #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  #  GNU General Public License for more details.
  #  http://www.gnu.org/licenses/gpl.html
    

  # Set up ------------------------------------------------------------------
    
  # Required functions
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
      
    return(dx)
  }
  
  Tintersect <- function(wf) {
    k1 <- wf - min(wf[1:which.max(wf)])
    xvar <- (1:length(k1) - 1)
    spl <- smooth.spline(k1 ~ xvar)
    newx <- which.max(diff(k1))
    pred0 <- predict(spl, x = newx, deriv = 0)
    pred1 <- predict(spl, x = newx, deriv = 1)
    yint <- pred0$y - (pred1$y * newx)
    xint <- (-yint / pred1$y)
    
    return(xint)
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
  
  # Load data
  sr <- 200
  p <- pressure[!is.na(pressure)]
  u <- flow[!is.na(flow)]
  
  # Make p & u signals the same length
  if (length(p) > length(u)) {
    p = p[1:length(u)]
  } else if (length(p) < length(u)) {
    u = u[1:length(p)]
  }
  
  # Remove diastole
  pensavg <- p - min(p[1:which.max(p)])
  uensavg <- u - min(u[1:which.max(u)])
  
  # Make time variable
  N = length(pensavg)
  tens = (0:(N - 1)) / sr # Time
  
  # Find foot with intersecting tangents
  xint.P <- Tintersect(pensavg)
  xint.U <- Tintersect(uensavg)
  
  xt.P <- tens[xint.P]
  xt.U <- tens[xint.U]
  
  # Isolate upstroke
  Imq = which.max(uensavg)
    
  # Max of linear fit
  lmmax  <- round(0.4*(Imq-round(xint.U)))
  
 
  # Align PQ ----------------------------------------------------------------
  
  if (isTRUE(align)) {
    maxlag <- round(xint.P - xint.U)
    
    if (maxlag != 0) {
      uensavg <-
        MattLag(uensavg, maxlag) # Lag flow so it aligns with pressure
    }
  }
    
    
  # WAVE INTENSITY ----------------------------------------------------------
    
  # Matt method
  # p = filter order
  # n = filter length (must be odd)
  # m = return mth derivative
  # ts = time scaling factor
  dpe = signal::sgolayfilt(
    pensavg,
    p = 2,
    n = 11,
    m = 1,
    ts = 1
  )
    
  due = signal::sgolayfilt(
    uensavg,
    p = 2,
    n = 11,
    m = 1,
    ts = 1
  )
  die = dpe * due
    
    
  # Estimate Zc & local complience ------------------------------------------
  
  # The slope of P ~ U in early systole = blood density (rho) * compliance (c)
  linea <- 
    lm(pensavg[round(xint.U):(lmmax+round(xint.U))] ~ uensavg[round(xint.U):(lmmax+round(xint.U))])
  rho = 1050                       # Blood density
  rhoc = linea$coefficients[[2]]   # Slope = Zc
  c = rhoc / rho                   # c = wave speed
  
  # Rhoc.ss is rhoc but calculated via the sum of squares method
  rhocSS = sqrt(sum(dpe[!is.na(dpe)] ^ 2, rm.na = TRUE) / sum(due[!is.na(due)] ^ 2, rm.na = TRUE))
  cSS = rhocSS / rho
  
    
  # Calculate forward & backward WI -----------------------------------------
  
  dpep = (dpe + rhoc * due) / 2 # Forward pressure difference eg. (P+Zc*U)/2
  dpem = (dpe - rhoc * due) / 2 # Backward pressure difference eg. (P-Zc*U)/2
  duem = -dpem / rhoc           # Backward velocity difference
  duep = dpep / rhoc            # Forward velocity difference
  diep = dpep * duep            # Forward wave intensity
  diem = dpem * duem            # Backward wave intensity
  
    
  # Wave separation ---------------------------------------------------------
  
  pep = (pensavg + uensavg * rhoc) / 2 # Pf
  pept = which.max(pep) / sr # time
  
  pem = (pensavg - uensavg * rhoc) / 2 # Pb
  pemt = which.max(pem) / sr           # time
  
  
  # Save Variables ----------------------------------------------------------
  
  TImaxp = which.max(pensavg)                 # time of max P
  
  lsys = which.min(dpe[TImaxp:length(dpe)])
  lsys = round(lsys) + 4                      # round + 4 samples for margin of error
  lsys = lsys + TImaxp
  
  
  # Find diacrotic notch ----------------------------------------------------
  
  dp <- fsg721(pensavg)
  nn <- which.min(dp)
  end <- length(p)
  
  dpdias <- fsg721(fsg721(p[nn:end]))
  fsect <- which.max(dpdias)
  dnotch <- (fsect + nn) + 1
  
  notch.t <- tens[dnotch]
  notch.p <- pensavg[dnotch]
    
    
  # Wave intensity variables ------------------------------------------------
  
  # Find max WI as per Alun
  w1 = pracma::findpeaks(die[1:lsys], npeaks = 1, minpeakheight = 0.7 *
                           max(die))
  w1val = w1[1]         # WI @ w1 (peak)
  w1loc = w1[2]         # index of w1
  
  # Find min in systole
  w0 = pracma::findpeaks(-die[w1loc:lsys],
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
  minpeak = max(die) / 20    # so peaks <5% of w1 ignored for w2
  w2 = pracma::findpeaks(die[dnotch:1], npeaks = 1, minpeakheight = minpeak)
  
  if (is.null(w2)) {
    isw2 = TRUE
    w2val = NA                  # WI @ w2
    w2loc = NA                  # index of w2
  } else {
    isw2 = FALSE
    w2val = w2[1]               # WI @ w2
    w2loc = dnotch - w2[2] + 1  # index of w2
  }
  
  
  # area is 1.06447 * height * width
  w1area = 1.06447 * w1val * length(die[w1[3]:w1[4]][die[w1[3]:w1[4]] > 0])
  
  #w0 area
  w0srt = w0[3] + w1loc
  w0end = w0[4] + w1loc
  if (isTRUE(isw0)) {
    w0area = 0
  } else {
    w0area = 1.06447 * w0val * length(die[w0srt:w0end][die[w0srt:w0end] < 0])
  }
    
  #w2 area
  w2srt = dnotch - w2[3] + 1
  w2end = dnotch - w2[4] + 1
  if (isTRUE(isw2)) {
    w2area = 0
  } else {
    w2area = 1.06447 * w2val * length(die[w2srt:w2end][die[w2srt:w2end] > 0])
  }
  
  wri = w0area / w1area
  

  # clean up calculated variables
  tmaxp = round(TImaxp / sr, 4)
  wi1 = round(w1val, 4)
  wit1 = round(w1loc / sr, 4)
  wia1 = round(w1area, 4)
  wi0 = round(w0val, 4)
  wit0 = round(w0loc / sr, 4)
  wia0 = round(w0area, 4)
  wi2 = round(w2val, 4)
  wit2 = round(w2loc / sr, 4)
  wia2 = round(w2area, 4)
  wrm = round(wri, 4)
  zc = round(rhoc, 4)
  zcss = round(rhocSS, 4)
  comp = round(c, 4)
  compss = round(cSS, 4)
  pfmax = round(max(pep), 4)
  pft = round(pept, 4)
  pbmax = round(max(pem), 4)
  pbt = round(pemt, 4)
  intpf = round(sum(pep) / sr, 4)
  intpb = round(sum(pem) / sr, 4)
  rm = round(pbmax / pfmax, 4)
  ri = round(pbmax / (pbmax + pfmax), 4)
  
    
  # PLOTS -------------------------------------------------------------------
  
  par(
    mfrow = c(2, 2),
    mar = c(3.5, 3.5, .5, .5),
    mgp = c(2, 1, 0))
  
  # Plot P and U
  plot(
    tens,
    pensavg,
    type = "l",
    lwd = 2,
    ylab = "Pressure & flow (mmHg)",
    xlab = "Time (s)",
    frame.plot = F)
  
  grid()
  
  par(new = TRUE)
  
  plot(
    tens,
    uensavg,
    type = "l",
    lwd = 2,
    col = "firebrick",
    xaxt = 'n',
    yaxt = 'n',
    ann = FALSE,
    frame.plot = F)
  
  abline(
    v = c(xt.P, notch.t),
    col = 'black',
    lwd = 2,
    lty = 3)
  
  # Plot Zc
  plot(
    uensavg,
    pensavg,
    ylab = "Pressure (mmHg)",
    xlab = "Flow velocity (ms)",
    pch = 21,
    bg = 'lightgrey',
    col = 'grey40',
    lwd = 0.9,
    cex = 1.2,
    frame.plot = F)
  
  grid()
  
  clip(uensavg[round(xint.U)],
       uensavg[lmmax + round(xint.U)],
       -1000,
       1000)
  
  abline(linea,
         col = "firebrick",
         lwd = 3)
  
  # Plot WI
  plot(
    tens,
    diem,
    type = "l",
    lwd = 2,
    col = "firebrick",
    lty = 1,
    xlab = "Time (s)",
    ylab = expression(paste("Wave Intensity (W/m" ^ "2" * ")")),
    ylim = c(min(diem, na.rm = T), max(die, na.rm = T)),
    frame.plot = F)
  
  grid()
  
  lines(tens,
        diep,
        col = "dodgerblue3",
        lwd = 2,
        lty = 1)
  
  lines(tens,
        die,
        col = 1,
        lwd = 2) 
  
  abline(v = tens[w1loc],
         col = "aquamarine4",
         lty = 3,
         lwd = 2)
  
  abline(v = tens[w0loc],
         col = "dodgerblue3",
         lty = 3,
         lwd = 2)
  
  abline(v = tens[w2loc],
         col = "firebrick",
         lty = 3,
         lwd = 2)
  
  # Plot Pf & Pb
  plot(
    tens,
    pensavg,
    type = "l",
    lwd = 2,
    ylab = "Pb & Pf (mmHg)",
    xlab = "Time (s)",
    ylim = c(min(pem), max(pensavg + 2)),
    frame.plot = F)
  
  grid()
  
  lines(tens,
        pep,
        col = "dodgerblue3",
        lwd = 2)
  
  lines(tens,
        pem,
        col = "firebrick",
        lwd = 2)
  
  
  # Results -----------------------------------------------------------------
  
  print(c('Sbp time........', tmaxp), quote = F)
  print(c('W1..............', wi1), quote = F)
  print(c('W1 time.........', wit1), quote = F)
  print(c('W1 area.........', wia1), quote = F)
  print(c('W0..............', wi0), quote = F)
  print(c('W0 time.........', wit0), quote = F)
  print(c('W0 area.........', wia0), quote = F)
  print(c('W2..............', wi2), quote = F)
  print(c('W2 time.........', wit2), quote = F)
  print(c('W2 area.........', wia2), quote = F)
  print(c('Reflect mag.....', wrm), quote = F)
  print(c('Zc PU...........', zc), quote = F)
  print(c('Zc SS...........', zcss), quote = F)
  print(c('Compliance PU...', comp), quote = F)
  print(c('Compliance SS...', compss), quote = F)
  print(c('Pb max..........', pbmax), quote = F)
  print(c('Pb integral.....', intpb), quote = F)
  print(c('Pb time.........', pbt), quote = F)
  print(c('Pf max..........', pfmax), quote = F)
  print(c('Pf integral.....', intpf), quote = F)
  print(c('Pf time.........', pft), quote = F)
  print(c('Reflct mag......', rm), quote = F)
  print(c('Reflct index....', ri), quote = F)
  
  
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
    
  df2 <-
    data.frame(
      tmaxp,
      wi1,
      wit1,
      wia1,
      wi0,
      wit0,
      wia0,
      wi2,
      wit2,
      wia2,
      wrm,
      zc,
      zcss,
      comp,
      compss,
      pbmax,
      intpb,
      pbt,
      pfmax,
      intpf,
      pft,
      rm,
      ri)
  
  Final_data <- list(df1, df2)
  
  return(Final_data)
  
}
