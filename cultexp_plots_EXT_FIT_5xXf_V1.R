##################################################################################
### Script for running, plotting, and fitting the EXT model to BAC- and BAC+ data#
### Tobias Reiner Vonnahme, UiT the Arctic University of Norway ##################
### v1.0 09.11.2020 ##############################################################
### Replication script to: https://doi.org/10.5194/bg-2020-314 ###################
### Data of this script under https://doi.org/10.18710/VA4IU9 ####################
##################################################################################
### Dependendy: Geider98v3.R #####################################################
###             FME package  #####################################################
###             deSolve package ##################################################
##################################################################################

# Input: t, x (Cf, Cr, Chl, N, r, DOC, resdoc, pcho, DON, TEPC, bactc, bactn, DIN)

####################################################################################
## 1) Parameter definition #########################################################
####################################################################################

#parameters START
#parameters and variables
#___________________________________________________________
# C   - Carbon content (starting condition, modelled)
# Chl - Chlorophyll content (starting condition, modelled)
# N   - Nitrogen content (starting condition, modelled)
# pc  - C-specific rate of photosynthesis (calculated)
# Pcmax - maximum value of pc at given temperature (calculated)
# Pcref - Value of PCmax at a given reference temperature (calculated)

SiPS2 <- 0.6#Photosynthesis reduction after Si limitation


Pcref <-   0.8
# cost - Cost of biosynthesis (parameter)
cost <- 1  
# Vcn  - max N uptake (calculated)
# Rc   - carbon based maintenance metabolic rate constant (parameter)
Rc <- 0.0958  # Rc of the G98 model (part wil become xf in the model)
xf <- 0.05463
xf2 <- xf * 5
#Rc <- RcG98 - xf
# Rchl and RN - Chl and N based maintenance metabolic rates(calculated)
# pchl - carbon specific light limited photosynt. rate (calculated)
# thetac - Chl : C ratio in the cell (calculated)
# thetan - Chl : N ratio in the cell (Calculated)
# thetanmax - max thetan (parameter)
thetanmax <- 1.7  
# Q    - N : C ratio in the cell (calculated)
# Qmin - minimum value of Q (parameter)
Qmin <- 0.05 
# Qmax - maximum value of Q (parameter)
Qmax <- 0.3 
# Alfchl - Chl specific linitial C assimilation rate (parameter)
Alfchl <- 0.04879
# I   - incident scalar irradiance (variable, model input)
I <- 100 # (measured)
# my  - specific growth rate (calculated)
# DIN - inorganic n in medium (variable, model input) in uM
# n   - shape factor describing dependence of Vcmax on Q
n <- 2.57249
# kn  - half saturation constant for nitrogen uptake (parameter)
kn <- 1 
kn2 <- 4.5372
nh4thres <- 0.4863#1.12#1.12 #1.12

Vmax <- 0.1 
ks <-7.6 
smin <- 1.82  

rem <-  5.391# remineralisation rate of excreted don (correction for bact C)
remd <-  0.000655 # remineralisation rate of background don (correction for bact C)
mubact <-0.8##0.04#0.8#0.04 # bacterial max growth rate
bactmax <- 0.7# 15/1000# 4e07/1e10#15/1000 # bacterial capacity
paramMod<-c(cost=cost, Rc=Rc, thetanmax=thetanmax, Qmin=Qmin, Qmax=Qmax, 
            Alfchl=Alfchl, I=I, n=n, kn=kn, Pcref=Pcref, xf=xf,
            ks=ks, Vmax=Vmax, smin=smin, rem=rem, remd=remd,
            mubact=mubact, bactmax=bactmax,
            kn2=kn2, nh4thres=nh4thres, SiPS=SiPS) 
pars<-paramMod
# parameter and variables END

tspan <- seq(1, 16, 0.01) # time of model run in days (16)

require(deSolve) # load package for solving differential equations
source('geider98v6.R') # load model function

########################################################################
### 2) Model fits to experiments #######################################
########################################################################
SiPS <-0.6

paramMod<-c(cost=cost, Rc=Rc, thetanmax=thetanmax, Qmin=Qmin, Qmax=Qmax, 
            Alfchl=Alfchl, I=I, n=n, kn=kn, Pcref=Pcref, xf=xf,
            ks=ks, Vmax=Vmax, smin=smin, rem=rem, remd=remd,
            mubact=mubact, bactmax=bactmax,
            kn2=kn2, nh4thres=nh4thres, SiPS=SiPS) 
pars<-paramMod


# Model fit to axenic experiment 
# Input: t, x (Cf, Cr, Chl, N, r, DOC, resdoc, pcho, DON, TEPC, bactc, bactn, DIN)
# define start variables
C0 <- 1.3   
Chl0 <-0.062 
N0 <-0.0892  
DIN0<-59.5
nh40 <- 5.5
cell0 <- 7.980
Si0 <- 19.3
PSi0 <- 7980  * 104 *10e-07 *3.9   
don0 <- 2800/6.625  
bact0 <- 0

# organize start variables into input vector
x0 <- c(C0, Chl0, N0, DIN0, PSi0, Si0, bact0, 0, 0, don0, nh40, 0)  
x <- x0

# run the model
out_axenic <- ode(x0, times = tspan, func = geider98v2, parms = paramMod, method= "rk")


# Model fit to bacteria enriched experiment
# define start variable 
C0 <- 1.3  # UNIT is mg L-1!!!
Chl0 <-0.062
N0 <-0.0892 
DIN0<-69.5 
nh40 <- 21.5 
cell0 <- 7.980
Si0 <- 22.15
PSi0 <- 7980  * 104 *10e-07 *3.9  #
don0 <- 2800/6.625
bact0 <- 35600 * 20e-09  # 35600 cells ml-1, 20fg C cell-10.18/1000 * 1e-12 fg mg-1 * 1e06 m3 ml-1   #0.18/1000

x0 <- c(C0, Chl0, N0, DIN0, PSi0, Si0, bact0, 0, 0, don0, nh40, 0) # 
x <- x0

out_bact <- ode(x0, times = tspan, func = geider98v2, parms = paramMod, method= "rk")

plot(tbactexp, bactexp * 20e-09)
lines(tspan,out_bact[,8])


################################
# Model runs with EPS aggregation
source('geider98v6_play2.R')

 SiPS <- SiPS2
paramMod<-c(cost=cost, Rc=Rc, thetanmax=thetanmax, Qmin=Qmin, Qmax=Qmax, 
            Alfchl=Alfchl, I=I, n=n, kn=kn, Pcref=Pcref, xf=xf,
            ks=ks, Vmax=Vmax, smin=smin, rem=rem, remd=remd,
            mubact=mubact, bactmax=bactmax,
            kn2=kn2, nh4thres=nh4thres, SiPS= SiPS, xf2 <- xf2) 
pars<-paramMod

# axenic
# Model fit to axenic experiment 
# Input: t, x (Cf, Cr, Chl, N, r, DOC, resdoc, pcho, DON, TEPC, bactc, bactn, DIN)
# define start variables
C0 <- 1.3   
Chl0 <-0.062 
N0 <-0.0892  
DIN0<-59.5
nh40 <- 5.5
cell0 <- 7.980
Si0 <- 19.3
PSi0 <- 7980  * 104 *10e-07 *3.9   
don0 <- 2800/6.625  
bact0 <- 0

source('geider98v6_play2.R')

x0 <- c(C0, Chl0, N0, DIN0, PSi0, Si0, bact0, 0, 0, don0, nh40, 0) # 
x <- x0

out_axenic_noEX <- ode(x0, times = tspan, func = geider98v2, parms = paramMod, method= "rk")

# with bact
# Model fit to bacteria enriched experiment
# define start variable 
C0 <- 1.3  # UNIT is mg L-1!!!
Chl0 <-0.062
N0 <-0.0892 
DIN0<-69.5 
nh40 <- 21.5 
cell0 <- 7.980
Si0 <- 22.15
PSi0 <- 7980  * 104 *10e-07 *3.9  #
don0 <- 2800/6.625
bact0 <- 35600 * 20e-09  # 35600 cells ml-1, 20fg C cell-10.18/1000 * 1e-12 fg mg-1 * 1e06 m3 ml-1   #0.18/1000


x0 <- c(C0, Chl0, N0, DIN0, PSi0, Si0, bact0, 0, 0, don0, nh40, 0) # 
x <- x0

out_bact_noEX <- ode(x0, times = tspan, func = geider98v2, parms = paramMod, method= "rk")



##########################




############################################################################################
### 3) Define measured values in the experimnents ##########################################
###########################################################################################

#### WITH BACTERIA
texp<-c(1,2,3,4,5,6,7,8,9,10,11,14,15,16)

cexp<-c(1314,1385,1477,2177,2531,3407,5066,5214,6208,4228,7177,8233,7030,8484)/1000
cexpmin<-c(977,1070,1332,1615,2119,2856,4580,4879,5542,4099,5862,8048,6563,6843)/1000
cexpmax<-c(2222,2003,1789,2190,2740,3502,5169,5782,6487,4925,7760,9524,8039,8955)/1000

nexp<-c(140,154,165,308,379,569,791,949,1089,728,1282,1563,1390,1538)/1000
nexpmin<-c(76,94,87,203,312,438,680,887,952,702,1093,1415,1365,1384)/1000
nexpmax<-c(187,238,202,330,412,603,900,1030,1180,876,1376,1790,1568,1744)/1000

chlexp<-c(63,70,108,330,250,370,550,596,617,670,615,560,604,418)/1000
chlexpmin<-c(57,65,104,165,231,363,534,147,585,651,567,518,588,410)/1000
chlexpmax<-c(67,74,116,349,256,370,567,602,629,681,653,598,622,446)/1000

tdinexp<-c(0,2,4,6,8,9,11,14,15)
dinexp<-c(70.3,69.2,49.3,37.0,37.6,34.1,16.5,4.7,0.0)
dinexpmin<-c(37.6,30.1,46.0,29.9,34.5,33.4,13.9,0.0,0.0)
dinexpmax<-c(90.2,69.9,50.3,48.1,41.8,34.8,19.1,17.5,0.1)

nh4exp<-c(19.65,20.5,11.2,7.7,16.3,14.5,12.2,12.8,13.3)
nh4expmin<-c(18.8,11.8,8.5,7.3,15.8,13.8,12.2,12.7,13.1)
nh4expmax<-c(20.5, 20.5,11.2,8.1,16.3,14.5,12.2,12.9,13.5)

tnh4exp<-c(0,
           1,
           2,
           3,
           4,
           5,
           6,
           7,
           8,
           10,
           11,
           15)
nh4expFull <-c(19.7,
               21.5,
               16.15,
               13.05,
               9.9,
               6.5,
               7.7,
               11.2,
               16.05,
               12,
               12.2,
               13.3)
nh4expminFull<-c(18.8,
                 19.0,
                 11.8,
                 11.9,
                 8.5,
                 6.3,
                 7.3,
                 11.2,
                 15.8,
                 11.7,
                 12.2,
                 13.1)
nh4expmaxFull<-c(20.5,
                 23.5,
                 20.5,
                 14.2,
                 11.2,
                 7.1,
                 8.1,
                 11.2,
                 16.3,
                 12.6,
                 12.2,
                 13.5)



po4exp<-c(29.1,18.2,14.2,16.3,27.9,38.4,33.9,47.0,11.3)
po4expmin<-c(16.3,16.8,14.0,15.1,16.6,36.7,31.7,2.1,4.1)
po4expmax<-c(52.4,18.5,14.7,16.6,43.7,40.2,36.1,53.3,18.5)

siexp<-c(23.5, 20.8,10.9,4.5,2.2,0.5,0.6,1.8,0.5)
siexpmin<-c(12.6, 9.2,10.0,3.8,1.9,0.3,0.5,0.3,0.4)
siexpmax<-c(29.5,21.4,11.6,5.0,4.0,0.7,0.8,2.2,0.7)

algexp<-c(6.9E+03,1.2E+04,1.5E+04,3.2E+04,3.8E+04,6.6E+04,8.7E+04,1.3E+05,
          1.3E+05 *1.3,7.6E+04*1.3,9.2E+04*1.3,7.2E+04*1.3,8.9E+04*1.3,7.2E+04*1.3)
algexpmin<-c(6.0E+03,1.1E+04,1.5E+04,2.8E+04,3.6E+04,6.4E+04,8.4E+04,1.2E+05,
             1.2E+05*1.3,6.1E+04*1.3,8.0E+04*1.3,5.7E+04*1.3,8.3E+04*1.3,7.0E+04*1.3)
algexpmax<-c(7.3E+03,1.3E+04,1.9E+04,3.4E+04,4.3E+04,6.7E+04,9.1E+04,1.4E+05,
             1.3E+05*1.3,9.2E+04*1.3,1.0E+05*1.3,8.7E+04*1.3,9.7E+04*1.3,8.2E+04*1.3)

tdonexp<-c(3,5,7,9,11,14,15)
donexp<-c(4.157, 3.145, 3.481, 3.419,3.523, 4.097, 3.223)
donexpmin<-c(3.864, 1.268, 2.208, 3.393, 3.421, 4.037, 2.745)
donexpmax<-c(4.530, 3.521, 4.514, 3.445, 3.626, 4.037, 3.552)

tbactexp<-c(1,2,3,4,5,6,7,8,9,10,11,14)
bactexp<-c(3.56E+04,2.69E+04,4.64E+04,2.63E+05,1.11E+06,2.1E+06,3.56E+06,
           8.1E+06,8E+06,2.7E+07,3E+07,7E+07)
bactexpmin<-c(2.13E+04,1.72E+04,2.55E+04,2.46E+05,6.78E+05,1.61E+06,3.11E+06,
              7E+06,4.06E+06,2.4E+07,9E+06,7E+07)
bactexpmax<-c(5.30E+04,3.96E+04,5.86E+04,2.75E+05,1.15E+06,2.50E+06,4.71E+06,
              8.97E+06,1E+07,2.8E+07,4E+07,7E+07)


###############################
### withOUT Bacteria #########

#time in days
texpout<- c(1,2,3,4,5,6,7,8,9,10,11,14,15,16)

#Total particluate Carbon ug/L (median, minimum, maximum)
cexpout<-c(1303.4,1419.4,1736.2,1752.0,1918.3,2505.3,2893.4,3841.3,4918.5,4870.2,4475.5,4126.4,4463.4,6078.9)/1000
cexpminout<-c(975.7,1127.7,1305.6,1324.7,1528.6,2334.2,2459.8,3710.3,4300.1,3969.3,4380.1,3946.5,4159.5,3481.9)/1000
cexpmaxout<-c(1450.4,2346.7,1848.3,1843.7,2164.8,2516.9,4304.5,4489.8,5222.4,5045.3,4489.3,4325.2,5708.1,6116.0)/1000
cerrout<- abs(scale(cexpmax - cexpmin))

#Total particulate Nitrogen ug/L (median, minimum, maximum)
nexpout<-c(89.2,136.2,194.9,182.1,242.8,381.9,473.6,604.4,739.6,808.2,514.7,422.2,471.5,601.1)/1000
nexpminout<-c(74.9,86.1,145.2,120.4,141.4,329.9,378.3,582.4,632.9,512.1,505.0,374.3,373.5,402.4)/1000
nexpmaxout<-c(139.5,235.7,203.7,186.8,279.7,408.6,491.6,691.4,764.7,1193.6,517.9,454.2,616.1,612.8)/1000
nerrout<- abs(scale(nexpmax -cexpmin))

#Total Chlorophyll ug/L (median, minimum, maximum)
chlexpout<-c(62.3,76.2,107.2,142.2,206.3,286.5,472.1,580.1,577.8,569.6,422.5,387.0,237.8,219.2)/1000
chlexpminout<-c(51.5,75.4,100.8,40.4,190.0,233.5,432.5,559.9,553.3,565.4,374.1,341.6,200.7,201.9)/1000
chlexpmaxout<-c(81.1,134.9,113.7,146.4,206.4,313.7,511.8,580.4,594.4,573.8,470.9,432.5,274.1,233.2)/1000
chlerrout<- abs(scale(chlexpmax - chlexpmin))

# algae cell counts (cells per ml)
algexpout<-c(8.0E+03,1.2E+04,1.6E+04,2.4E+04,3.2E+04,5.4E+04,6.9E+04,9.0E+04,
             1.1E+05,9.2E+04,9.9E+04,5.5E+04,9.1E+04,7.2E+04)
algexpminout<-c(6.6E+03,9.5E+03,1.5E+04,2.0E+04,2.1E+04,5.3E+04,6.7E+04,8.8E+04,
                1.0E+05,8.4E+04,9.1E+04,4.3E+04,8.8E+04,7.0E+04)
algexpmaxout<-c(8.8E+03,1.3E+04,1.6E+04,2.6E+04,3.8E+04,5.6E+04,7.1E+04,9.2E+04,
                1.3E+05,9.9E+04,1.1E+05,6.6E+04,9.5E+04,8.7E+04)

# bacteria cell counts (cells per ml)
tbactexpout<-c(1,2,3,4,5,6,7,8,9,10,11,14)
bactexpout<-c(2.7E+03,1.2E+03,1.5E+03,6.7E+03,2.3E+03,4.7E+03,1.5E+04,3.8E+03,
              2.5E+04,1.4E+05,4.6E+04,3.0E+06)
bactexpminout<-c(1.3E+03,1.0E+03,1.5E+03,3.8E+03,1.0E+04,2.7E+03,7.6E+03,3.5E+03,
                 1.6E+04,4.0E+04,2.9E+04,3.0E+06)
bactexpmaxout<-c(4.1E+03,2.9E+03,1.5E+03,4.2E+04,4.2E+03,6.7E+03,2.1E+04,4.3E+03,
                 5.3E+04,2.3E+05,6.4E+04,3.0E+06)

tdonexpout<-c(3,5,7,9,11,14,15)
donexpout<-c(1.838, 2.028, 1.424, 1.612, 1.7265 ,1.743, 1.772)
donexpminout<-c(1.661, 1.334, 0.987, 1.162, 1.573, 1.571, 1.398)
donexpmaxout<-c(2.075, 3.007, 2.030, 1.889, 1.880,1.926, 2.133)


# measublue NO3 and NO2 (uM) 
tdinexpout<- c(0,2,4,6,8,9,11,14,15)
dinexpout<-c(70.3,49.1,46.8,26.9,2.5,0.5,0.3,0.2,0.7)
dinexpminout<-c(37.6,33.4,40.1,18.1,1.9,0.4,0.1,0.1,0.0)
dinexpmaxout<-c(90.2,63.7,47.7,28.5,3.4,0.7,1.0,0.3,1.9)


#measublue nh4 (uM) (max >100 are removed, too high filtration pressure)
nh4expout<-c(19.65,4.4,4.3,2.6,5.8,3.25,2.1,15.2,8.0)
nh4expminout<-c(18.8,2.6,2.3,2.3,4.8,2.6,1.5,15.2,7.1)
nh4expmaxout<-c(20.5,7.7, 8.0, 3.3,  8.8, 12.1, 2.5, 15.2, 9.0)

tnh4expout<-c(0,
              1,
              2,
              3,
              4,
              5,
              6,
              7,
              8,
              9,
              10,
              11,
              14,
              15)
nh4expoutFull <-c(19.7,
                  5.5,
                  4.4,
                  2.0,
                  4.3,
                  4.3,
                  2.6,
                  1.9,
                  5.8,
                  3.25,
                  10.4,
                  2.1,
                  15.2,
                  9.0)
nh4expminoutFull<-c(18.8,
                    3.2,
                    2.6,
                    1.4,
                    2.3,
                    2.7,
                    2.3,
                    1.7,
                    4.8,
                    2.6,
                    10.1,
                    1.5,
                    15.2,
                    7.1)
nh4expmaxoutFull<-c(20.5,
                    7.7,
                    5.7,
                    8.0,
                    6.6,
                    4.5,
                    3.3,
                    4.5,
                    8.8,
                    3.9,
                    12.1,
                    2.5,
                    15.2,
                    9)



DINexpout<-dinexpout + nh4expout

# measublue PO4 (uM)
po4expout<-c(29.1,15.2,13.0,13.2,17.5,14.5,14.5,13.1,15.3)
po4expminout<-c(16.3,15.0,12.4,12.9,13.8,13.4,14.4,12.8,12.2)
po4expmaxout<-c(52.4,16.3,13.3,14.6,24.4,14.8,14.5,13.4,21.3)

#measublue Si(OH)4 in uM
siexpout<-c(23.5,15.1,13.3,3.8,0.8,0.7,0.7,0.7,0.6)
siexpminout<-c(12.6,11.1,11.5,2.9,0.1,0.2,0.7,0.6,0.3)
siexpmaxout<-c(29.5,20.3,13.8,5.6,1.1,1.0,1.5,0.9,0.9)


#############################################################
###############################################################

chloutlier<-chlexpmin[8]
chlexpmin[8]<-chlexp[8]
chlexp[8]<-mean(chlexp[8], chlexpmax[8])


DINexp<-dinexp+nh4exp
DINexpmax<-dinexpmax+nh4expmax
DINexpmin<-dinexpmin+nh4expmin
DINexpout<-dinexpout+nh4expout
DINexpmaxout<-dinexpmaxout+nh4expmaxout
DINexpminout<-dinexpminout+nh4expminout
tDINexp <- tdinexp
tDINexpout <- tdinexpout
# DIN cells plot


#########################################################3
##########################################################
#########################################################


pdf(file="FigS2_modpap_rev3.pdf", width=12, height=6)
#tiff(filename="Fig4_modpaper_rev2.tiff", width=600, compression="none")

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(1,2))

plot(texp, cexpmax, type='n', 
     ylab = expression(paste("POC [",mu,"g mL"^"-1","]")), 
     xlab = "time [days]",ylim=c(0,12), xlim=c(0,16),
     axes=F)
axis(2, las=1, at=c(0,2,4,6,8,10,12))
axis(1, las=2, at=c(0,2,4,6,8,10,12,14,16))
lines(texp, cexpmin, col="grey")
lines(texp, cexpmax, col="grey")
polygon(c(texp, rev(texp)), c(cexpmax, rev(cexpmin)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(texp, cexp, type="p", bg="red", pch=21)
points(texpout, cexpmaxout, type='n')
lines(texpout, cexpminout, col="grey")
lines(texpout, cexpmaxout, col="grey")
polygon(c(texpout, rev(texpout)), c(cexpmaxout, rev(cexpminout)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(texpout, cexpout, type="p", bg="blue", pch=21)
title(main="(a)",adj=0, cex=0.75)
abline(v=8, col="grey")
lines(tspan, out_axenic[,2], col="blue")
lines(tspan, out_bact[,2], col="red")
lines(tspan, out_axenic_noEX[,2], col="blue", lty="dotted", lwd=3)
lines(tspan, out_bact_noEX[,2], col="red",lty="dotted", lwd=3)
legend("topleft", legend=c("BAC+ exp","BAC- exp", "EXT BAC+","EXT BAC-", 
                           "EXT +5Xf BAC+","EXT +5Xf BAC-"), 
       col=c("black","black","red","blue","red","blue"), 
       pt.bg=c("red","blue",NA,NA,NA,NA),
       lty=c(0,0,1,1,3,3), pch=c(21,21,NA,NA,NA,NA),
       box.lty=0, bg="transparent", cex=0.7)

# PON plot
plot(texp, nexpmax, type='n', 
     ylab = expression(paste("PON [",mu,"g mL"^"-1","]")), 
     xlab = "time [days]", xlim=c(0,16),ylim=c(0,3),
     axes=F)
axis(2, las=1, at=c(0,0.5,1,1.5,2,2.5,3))
axis(1, las=2, at=c(0,2,4,6,8,10,12,14,16))
lines(texp, nexpmin, col="grey")
lines(texp, nexpmax, col="grey")
polygon(c(texp, rev(texp)), c(nexpmax, rev(nexpmin)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(texp, nexp, type="p", bg="red", pch=21)
points(texpout, nexpmaxout, type='n')
lines(texpout, nexpminout, col="grey")
lines(texpout, nexpmaxout, col="grey")
polygon(c(texpout, rev(texpout)), c(nexpmaxout, rev(nexpminout)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(texpout, nexpout, type="p", bg="blue", pch=21)
title(main="(b)",adj=0, cex=0.75)
abline(v=8, col="grey")
lines(tspan, out_axenic[,4], col="blue")
lines(tspan, out_bact[,4], col="red")
lines(tspan, out_axenic_noEX[,4]  , col="blue", lty="dotted", lwd=3)
lines(tspan, out_bact_noEX[,4] , col="red",lty="dotted", lwd=3)


dev.off()


#####Model fit to experiment without DOM excretion
cmod <- out_bact_noEX[c(1,101,201,301, 401, 501, 601, 701, 801, 901, 1001,  1301, 1401, 1501),2]
costfc <- sqrt(sum((cexp - cmod)^2 ))
costfc <- sqrt(sum(scale(cexp - cmod)^2 ))
# Chl
chlmod <- out_bact_noEX[c(1,101,201,301, 401, 501, 601, 701, 801, 901, 1001,  1301, 1401, 1501),3]
costfchl <- sqrt( sum((chlexp - chlmod)^2 ) )
costfchl <- sqrt( sum(scale(chlexp - chlmod)^2 ) )
# N
nmod <- out_bact_noEX[c(1,101,201,301, 401, 501, 601, 701, 801, 901, 1001,  1301, 1401, 1501),4]
  
costfn <- sqrt( sum( (nexp - nmod)^2 ) )
costfn <- sqrt( sum(scale(nexp - nmod)^2 ) )
# DIN
# t= 2  4  6  8  9 11 14 15
dinmod <- out_bact_noEX[c(101,301, 501,  701, 801,  1001,  1301, 1401),5]+
  out_bact_noEX[c(101,301, 501,  701, 801,  1001,  1301, 1401),12]
costfdin <- sqrt( sum( (DINexp[-1] - dinmod)^2 ) )
costfdin <- sqrt( sum(scale(DINexp[-1] - dinmod)^2 ) )
#Total
costTot <- sqrt( sum( (nmod-nexp)^2/var(nexp) ) +
                   sum( (chlmod - chlexp)^2 / var(chlexp)) +
                   sum( (cmod - cexp)^2 / var(cexp))+
                   sum( (dinmod -  DINexp[-1])^2 / var(DINexp[-1]))) 
costTot


