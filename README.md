# Wave intensity/separation analysis in R

The R script herein separates the blood pressure waveform into its forward and backward traveling components using simultaneously measured blood pressure and flow waves sampled at 200hz.

To load the function:
```R
devtools::source_url("https://raw.githubusercontent.com/mkarmstrong/WaveIntensity/main/WaveAnalyses.R")
```


To run the function:
```R
ans <- WaveAnalyses(mydata$pressure, 
                    mydata$flow, 
                    lowpass = FALSE, # apply a lowpass filter to the pressure wave?
                    align = FALSE,   # time align the pressure & flow wave?
                    plot = TRUE)     # Plot results?
```

Set `align = TRUE` if presure and flow curves require time alignment.
<br/><br/>

If `plot = TRUE` Results will be sent to the plot tab for inspection:

![alt text](WI_plot.png)

<br/><br/>

The function returns the following values:

**Variable**      | **Description**
------------------|-------------------------
sbp               | Systolic BP
Tsbp              | Time to systolic BP
dbp               | Diastolic BP
Tdbp              | Time to diastolic BP
p1                | BP at p1
Tp1               | Time to p1
p2                | BP at p2
Tp2               | Time to p2
ed                | BP at incisura
Ted               | Ejection duration
p3                | BP at incisura peak
Tp3               | Time to incisura peak
Aix               | Augmentation index
qmax              | Peak flow
Tqmax             | Time to peak flow
qfoot             | Flow at end diastole
qmean             | Mean flow
wi1               | Peak forward compresion wave
Twi1              | Time to peak FCW
Awi1              | Area under FCW
wi0               | Peak backward compresion wave
Twi0              | Time to peak BCW
Awi0              | Area under BCW
wi2               | Peak forward decompresion wave
Twi2              | Time to peak FDW
Awi2              | Area under FDW
wrm               | Wave reflection magnitude
zc                | Characteristic impedance (pu loop)
zcss              | Characteristic impedance (sum of squares)
comp              | Arterial complience (pu loop)
compss            | Arterial complience (sum of squares)
pfmax             | Peak forward BP
Tpfmax            | Time to peak Pf
Apf               | Area under Pf
pbmax             | Peak backward BP         
Tpbmax            | Time to peak Pb
Apb               | Area under Pb
rm                | Reflection magnitude
ri                | Resistance index

<br/><br/>

Enjoy!
