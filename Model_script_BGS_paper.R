##################################################################################
### Script for running, plotting, and fitting the EXT model to BAC- and BAC+ data
### Tobias Reiner Vonnahme, UiT the Arctic University of Norway #################
### v1.0 09.11.2020 #############################################################
### Replication script to: https://doi.org/10.5194/bg-2020-314 ##################
#################################################################################
### Geider98v3 needed ###########################################################



# Input: t, x (Cf, Cr, Chl, N, r, DOC, resdoc, pcho, DON, TEPC, bactc, bactn, DIN)

#parameters START

#parameters and variables
#___________________________________________________________
# C   - Carbon content (starting condition, modelled)
# Chl - Chlorophyll content (starting condition, modelled)
# N   - Nitrogen content (starting condition, modelled)
# pc  - C-specific rate of photosynthesis (calculated)
# Pcmax - maximum value of pc at given temperature (calculated)
# Pcref - Value of PCmax at a given reference temperature (calculated)
Pcref <-   0.8#
# cost - Cost of biosynthesis (parameter)
cost <- 1 # 
# Vcn  - max N uptake (calculated)
# Rc   - carbon based maintenance metabolic rate constant (parameter)
Rc <- 0.01 # 
xf <- 0.06
# Rchl and RN - Chl and N based maintenance metabolic rates(calculated)
# pchl - carbon specific light limited photosynt. rate (calculated)
# thetac - Chl : C ratio in the cell (calculated)
# thetan - Chl : N ratio in the cell (Calculated)
# thetanmax - max thetan (parameter)
thetanmax <- 1.7 # 
# Q    - N : C ratio in the cell (calculated)
# Qmin - minimum value of Q (parameter)
Qmin <- 0.05 # 
# Qmax - maximum value of Q (parameter)
Qmax <- 0.3 #
# Alfchl - Chl specific linitial C assimilation rate (parameter)
Alfchl <- 0.076 # 
# I   - incident scalar irradiance (variable, model input)
I <- 100 # (measured)
# my  - specific growth rate (calculated)
# DIN - inorganic n in medium (variable, model input) in uM
#DIN <- 35
# n   - shape factor describing dependence of Vcmax on Q
n <- 3.45 # 
# kn  - half saturation constant for nitrogen uptake (parameter)
kn <- 2 #
kn2 <- 6.97#
nh4thres <- 1.19#
SiPS <- 0.20##Photosynthesis reduction after Si limitation

Vmax <- 0.1 #
ks <-7.6 #
smin <- 1.82 # 

rem <-  10 # remineralisation rate of excreted don (correction for bact C)
remd <-  4.86#4.83#4.9 # remineralisation rate of background don (correction for bact C)
mubact <-0.04 # bacterial max growth rate
bactmax <-15/1000 # bacterial capacity
paramMod<-c(cost=cost, Rc=Rc, thetanmax=thetanmax, Qmin=Qmin, Qmax=Qmax, 
            Alfchl=Alfchl, I=I, n=n, kn=kn, Pcref=Pcref, xf=xf,
            ks=ks, Vmax=Vmax, smin=smin, rem=rem, remd=remd,
            mubact=mubact, bactmax=bactmax,
            kn2=kn2, nh4thres=nh4thres, SiPS=SiPS) 
pars<-paramMod
# state variable END

tspan <- seq(1, 16, 0.01)

require(deSolve)
source('geider98v3.R')

# axenic
# Input: t, x (Cf, Cr, Chl, N, r, DOC, resdoc, pcho, DON, TEPC, bactc, bactn, DIN)
# state variable START
C0 <- 1203  / 1000
Chl0 <-62.3  / 1000
N0 <-89.2  / 1000
DIN0 <- 89#
DIN0<-49.1#
nh40 <- 4.4#
cell0 <- 7980 /1000
Si0 <- 23
PSi0 <- 7980  * 104 *10e-07 *3.9  # 
don0 <- 2800/16  
bact0 <- 0#

x0 <- c(C0, Chl0, N0, DIN0, PSi0, Si0, bact0, 0, 0, don0, nh40) # 
x <- x0

out_axenic <- ode(x0, times = tspan, func = geider98v2, parms = paramMod, method= "rk")


# with bact

# state variable START
C0 <- 1203  / 1000
Chl0 <-62.3  / 1000
N0 <-89.2  / 1000
DIN0 <- 89#
DIN0<-56.4 #
nh40 <- 20.5 #
cell0 <- 7980 /1000
Si0 <- 23
PSi0 <- 7980  * 104 *10e-07 *3.9  #
don0 <- 2800/16  
bact0 <- 0.18/1000

x0 <- c(C0, Chl0, N0, DIN0, PSi0, Si0, bact0, 0, 0, don0, nh40) # 
x <- x0

out_bact <- ode(x0, times = tspan, func = geider98v2, parms = paramMod, method= "rk")

################################

Rc <- Rc + xf #
xf <- 0.0#

paramMod<-c(cost=cost, Rc=Rc, thetanmax=thetanmax, Qmin=Qmin, Qmax=Qmax, 
            Alfchl=Alfchl, I=I, n=n, kn=kn, Pcref=Pcref, xf=xf,
            ks=ks, Vmax=Vmax, smin=smin, rem=rem, remd=remd,
            mubact=mubact, bactmax=bactmax,
            kn2=kn2, nh4thres=nh4thres, SiPS= SiPS) 
pars<-paramMod

# axenic
# Input: t, x (Cf, Cr, Chl, N, r, DOC, resdoc, pcho, DON, TEPC, bactc, bactn, DIN)
# state variable START
C0 <- 1203  / 1000
Chl0 <-62.3  / 1000
N0 <-89.2  / 1000
DIN0 <- 89#
DIN0<-49.1#
nh40 <- 4.4#
cell0 <- 7980 /1000
Si0 <- 23
PSi0 <- 7980  * 104 *10e-07 *3.9  #
don0 <- 2800/16  
bact0 <- 0#

x0 <- c(C0, Chl0, N0, DIN0, PSi0, Si0, bact0, 0, 0, don0, nh40) # 
x <- x0

out_axenic_noEX <- ode(x0, times = tspan, func = geider98v2, parms = paramMod, method= "rk")

# with bact

# state variable START
C0 <- 1203  / 1000
Chl0 <-62.3  / 1000
N0 <-89.2  / 1000
DIN0 <- 89#
DIN0<-56.4 #
nh40 <- 20.5 #
cell0 <- 7980 /1000
Si0 <- 23
PSi0 <- 7980  * 104 *10e-07 *3.9  #
don0 <- 2800/16  
bact0 <- 0.18/1000

x0 <- c(C0, Chl0, N0, DIN0, PSi0, Si0, bact0, 0, 0, don0, nh40) # 
x <- x0

out_bact_noEX <- ode(x0, times = tspan, func = geider98v2, parms = paramMod, method= "rk")



##########################




# measured values in the experiment

### Model end
# measured values in the experiment
#time in days


#### HERE WITH BACTERIA!!!!
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
### withOUT Bacteria

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

DINexpout<-dinexpout + nh4expout

# measublue PO4 (uM)
po4expout<-c(29.1,15.2,13.0,13.2,17.5,14.5,14.5,13.1,15.3)
po4expminout<-c(16.3,15.0,12.4,12.9,13.8,13.4,14.4,12.8,12.2)
po4expmaxout<-c(52.4,16.3,13.3,14.6,24.4,14.8,14.5,13.4,21.3)

#measublue Si(OH)4 in uM
siexpout<-c(23.5,15.1,13.3,3.8,0.8,0.7,0.7,0.7,0.6)
siexpminout<-c(12.6,11.1,11.5,2.9,0.1,0.2,0.7,0.6,0.3)
siexpmaxout<-c(29.5,20.3,13.8,5.6,1.1,1.0,1.5,0.9,0.9)

##############################################################
# POC plot

#pdf(file="Fig5_modpap_rev2.pdf")#, width=600, height=750)
tiff(file="Fig5_modpap_rev2.tiff", width=600, height=750, compression="none")

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(3,2))

#par(mfrow=c(1,2))

plot(texp, cexpmax, type='n', 
     ylab = expression(paste("Carbon [",mu,"g mL"^"-1","]")), 
     xlab = "time [days ]",ylim=c(0,10), xlim=c(0,15),
     axes=F)
axis(2, las=1, at=c(0,2,4,6,8,10))
axis(1, las=1, at=c(0,2,4,6,8,10,12,14))
lines(texp, cexpmin, col="grey")
lines(texp, cexpmax, col="grey")
polygon(c(texp, rev(texp)), c(cexpmax, rev(cexpmin)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(texp, cexp, type="p", col="red")
points(texpout, cexpmaxout, type='n')
lines(texpout, cexpminout, col="grey")
lines(texpout, cexpmaxout, col="grey")
polygon(c(texpout, rev(texpout)), c(cexpmaxout, rev(cexpminout)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(texpout, cexpout, type="p", col="blue")
title(main="a)",adj=0, cex=0.75)
abline(v=8, col="grey")
lines(tspan, out_axenic[,2], col="blue")
lines(tspan, out_bact[,2], col="red")
lines(tspan, out_axenic_noEX[,2], col="blue", lty="dotted")
lines(tspan, out_bact_noEX[,2], col="red",lty="dotted")
legend("topleft", legend=c("BAC+ exp","BAC- exp", "EXT BAC+","EXT BAC-", "EXT -excr BAC+","EXT -excr BAC-"), 
       col=c("red","blue","red","blue","red","blue"), 
       lty=c(0,0,1,1,3,3), pch=c(1,1,NA,NA,NA,NA),
       box.lty=0, bg="transparent", cex=0.7)

# PON plot
plot(texp, nexpmax, type='n', 
     ylab = expression(paste("PON [",mu,"g mL"^"-1","]")), 
     xlab = "time [days ]", xlim=c(0,15),ylim=c(0,2),
     axes=F)
axis(2, las=1, at=c(0,0.5,1,1.5,2))
axis(1, las=1, at=c(0,2,4,6,8,10,12,14))
lines(texp, nexpmin, col="grey")
lines(texp, nexpmax, col="grey")
polygon(c(texp, rev(texp)), c(nexpmax, rev(nexpmin)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(texp, nexp, type="p", col="red")
points(texpout, nexpmaxout, type='n')
lines(texpout, nexpminout, col="grey")
lines(texpout, nexpmaxout, col="grey")
polygon(c(texpout, rev(texpout)), c(nexpmaxout, rev(nexpminout)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(texpout, nexpout, type="p", col="blue")
title(main="(b)",adj=0, cex=0.75)
abline(v=8, col="grey")
lines(tspan, out_axenic[,4], col="blue")
lines(tspan, out_bact[,4], col="red")
lines(tspan, out_axenic_noEX[,4], col="blue", lty="dotted")
lines(tspan, out_bact_noEX[,4], col="red",lty="dotted")

# Chl plot
plot(texp, chlexpmax, type='n', 
     ylab = expression(paste("Chlorophyll a [",mu,"g mL"^"-1","]")), 
     xlab = "time [days ]", xlim=c(0,15),ylim=c(0,0.7),
     axes=F)
axis(2, las=1, at=c(0,0.2,0.4,0.6))
axis(1, las=1, at=c(0,2,4,6,8,10,12,14))
lines(texp, chlexpmin, col="grey")
lines(texp, chlexpmax, col="grey")
polygon(c(texp, rev(texp)), c(chlexpmax, rev(chlexpmin)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(texp, chlexp, type="p", col="red")
points(texpout, chlexpmaxout, type='n')
lines(texpout, chlexpminout, col="grey")
lines(texpout, chlexpmaxout, col="grey")
polygon(c(texpout, rev(texpout)), c(chlexpmaxout, rev(chlexpminout)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(texpout, chlexpout, type="p", col="blue")
title(main="(c)",adj=0, cex=0.75)
abline(v=8, col="grey")
lines(tspan, out_axenic[,3], col="blue")
lines(tspan, out_bact[,3], col="red")
lines(tspan, out_axenic_noEX[,3], col="blue", lty="dotted")
lines(tspan, out_bact_noEX[,3], col="red",lty="dotted")



# POC:PON plot

plot(texp, cexpmax/nexpmin, type='n', 
     ylab = expression(paste("C : N ratio [",mu,"g C : ",mu,"g N]")), 
     xlab = "time [days ]",ylim=c(0,25), xlim=c(0,16),
     axes=F)
axis(2, las=1, at=c(0,5,10,15,20))
axis(1, las=1, at=c(0,2,4,6,8,10,12,14))
lines(texp, cexpmin/nexpmax, col="grey")
lines(texp, cexpmax/nexpmin, col="grey")
polygon(c(texp, rev(texp)), c(cexpmax/nexpmin, rev(cexpmin/nexpmax)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(texp, cexp/nexp, type="p", col="red")
points(texpout, cexpmaxout/nexpminout, type='n')
lines(texpout, cexpminout/nexpmaxout, col="grey")
lines(texpout, cexpmaxout/nexpminout, col="grey")
polygon(c(texpout, rev(texpout)), c(cexpmaxout/nexpminout, rev(cexpminout/nexpmaxout)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(texpout, cexpout/nexpout, type="p", col="blue")
title(main="(d)",adj=0, cex=0.75)
abline(v=8, col="grey")
lines(tspan, out_axenic[,2]/out_axenic[,4], col="blue")
lines(tspan, out_bact[,2]/out_bact[,4], col="red")
lines(tspan, out_axenic_noEX[,2]/out_axenic_noEX[,4], col="blue", lty="dotted")
lines(tspan, out_bact_noEX[,2]/out_bact_noEX[,4], col="red", lty="dotted")


# C:Chl plot

plot(texp, cexpmax/chlexpmin, type='n', 
     ylab = expression(paste("C : Chl ratio [",mu,"g C : ",mu,"g Chl]")), 
     xlab = "time [days ]",ylim=c(0,40), xlim=c(0,16),
     axes=F)
axis(2, las=1, at=c(0,10,20,30,40))
axis(1, las=1, at=c(0,2,4,6,8,10,12,14))
lines(texp, cexpmin/chlexpmax, col="grey")
lines(texp, cexpmax/chlexpmin, col="grey")
polygon(c(texp, rev(texp)), c(cexpmax/chlexpmin, rev(cexpmin/chlexpmax)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(texp, cexp/chlexp, type="p", col="red")
points(texpout, cexpmaxout/chlexpminout, type='n')
lines(texpout, cexpminout/chlexpmaxout, col="grey")
lines(texpout, cexpmaxout/chlexpminout, col="grey")
polygon(c(texpout, rev(texpout)), c(cexpmaxout/chlexpminout, rev(cexpminout/chlexpmaxout)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(texpout, cexpout/chlexpout, type="p", col="blue")
title(main="(e)",adj=0, cex=0.75)
abline(v=8, col="grey")
lines(tspan, out_axenic[,2]/out_axenic[,3], col="blue")
lines(tspan, out_bact[,2]/out_bact[,3], col="red")
lines(tspan, out_axenic_noEX[,2]/out_axenic_noEX[,3], col="blue", lty="dotted")
lines(tspan, out_bact_noEX[,2]/out_bact_noEX[,3], col="red", lty="dotted")


# N/Chl plot

plot(texp, nexpmax/chlexpmin, type='n', 
     ylab = expression(paste("N : Chl ratio [",mu,"g N : ",mu,"g Chl]")), 
     xlab = "time [days ]",ylim=c(0,6.2), xlim=c(0,16),
     axes=F)
axis(2, las=1, at=c(0,1,2,3,4,5))
axis(1, las=1, at=c(0,2,4,6,8,10,12,14))
lines(texp, nexpmin/chlexpmax, col="grey")
lines(texp, nexpmax/chlexpmin, col="grey")
polygon(c(texp, rev(texp)), c(nexpmax/chlexpmin, rev(nexpmin/chlexpmax)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(texp, nexp/chlexp, type="p", col="red")
points(texpout, nexpmaxout/chlexpminout, type='n')
lines(texpout, nexpminout/chlexpmaxout, col="grey")
lines(texpout, nexpmaxout/chlexpminout, col="grey")
polygon(c(texpout, rev(texpout)), c(nexpmaxout/chlexpminout, rev(nexpminout/chlexpmaxout)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(texpout, nexpout/chlexpout, type="p", col="blue")
title(main="(f)",adj=0, cex=0.75)
abline(v=8, col="grey")
lines(tspan, out_axenic[,4]/out_axenic[,3], col="blue")
lines(tspan, out_bact[,4]/out_bact[,3], col="red")
lines(tspan, out_axenic_noEX[,4]/out_axenic_noEX[,3], col="blue", lty="dotted")
lines(tspan, out_bact_noEX[,4]/out_bact_noEX[,3], col="red", lty="dotted")


dev.off()


##########################################################
#pdf(file="Fig6_modpap_rev2.pdf")#, width=600)
tiff(file="Fig6_modpap_rev2.tiff", width=600, compression="none")
par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))

# DIN cells plot

DINexp<-dinexp+nh4exp
DINexpmax<-dinexpmax+nh4expmax
DINexpmin<-dinexpmin+nh4expmin
DINexpout<-dinexpout+nh4expout
DINexpmaxout<-dinexpmaxout+nh4expmaxout
DINexpminout<-dinexpminout+nh4expminout

plot(tdinexp, DINexpmax, type='n', 
     ylab = expression(paste("DIN (NO "[X],"  + NH "[4]," [",mu,"mol L"^"-1","]")), 
     xlab = "time [days]", ylim=c(0,90), xlim=c(0,16),
     axes=F)
axis(2, las=1, at=c(0,20,40,60,80))
axis(1, las=1, at=c(0,2,4,6,8,10,12,14))
title(main="(a)",adj=0, cex=0.75)
lines(tdinexp, DINexpmin, col="grey")
lines(tdinexp, DINexpmax, col="grey")
polygon(c(tdinexp, rev(tdinexp)), c(DINexpmax, rev(DINexpmin)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(tdinexp, DINexp, type="p", col="red")
points(tdinexpout, DINexpmaxout, type='n')
lines(tdinexpout, DINexpminout, col="grey")
lines(tdinexpout, DINexpmaxout, col="grey")
polygon(c(tdinexpout, rev(tdinexpout)), c(DINexpmaxout, rev(DINexpminout)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(tdinexpout, DINexpout, type="p", col="blue")
abline(v=8, col="grey")
#abline(h=0.5, lty="dotted")
lines(tspan, out_axenic[,5]+out_axenic[,12], col="blue")
lines(tspan, out_bact[,5]+out_bact[,12], col="red")
lines(tspan, out_axenic_noEX[,5]+out_axenic_noEX[,12], col="blue", lty="dotted")
lines(tspan, out_bact_noEX[,5]+out_bact_noEX[,12], col="red", lty="dotted")


# NOx cells plot
plot(tdinexp, dinexpmax, type='n', 
     ylab = expression(paste("NO "[X]," [",mu,"mol L"^"-1","]")), 
     xlab = "time [days]", ylim=c(0,90), xlim=c(0,16),
     axes=F)
axis(2, las=1, at=c(0,20,40,60,80))
axis(1, las=1, at=c(0,2,4,6,8,10,12,14))
title(main="(b)",adj=0, cex=0.75)
lines(tdinexp, dinexpmin, col="grey")
lines(tdinexp, dinexpmax, col="grey")
polygon(c(tdinexp, rev(tdinexp)), c(dinexpmax, rev(dinexpmin)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(tdinexp, dinexp, type="p", col="red")
points(tdinexpout, dinexpmaxout, type='n')
lines(tdinexpout, dinexpminout, col="grey")
lines(tdinexpout, dinexpmaxout, col="grey")
polygon(c(tdinexpout, rev(tdinexpout)), c(dinexpmaxout, rev(dinexpminout)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(tdinexpout, dinexpout, type="p", col="blue")
abline(v=8, col="grey")
#abline(h=0.5, lty="dotted")
lines(tspan, out_axenic[,5], col="blue")
lines(tspan, out_bact[,5], col="red")
lines(tspan, out_axenic_noEX[,5], col="blue",lty="dotted")
lines(tspan, out_bact_noEX[,5], col="red",lty="dotted")


# NH4  plot
plot(tdinexp, nh4expmax, type='n', 
     ylab = expression(paste("NH "[4]," [",mu,"mol L"^"-1","]")), 
     xlab = "time [days ]", ylim=c(0,21), xlim=c(0,16), 
     axes=F)
axis(2, las=1, at=c(0,10,20))
axis(1, las=1, at=c(0,2,4,6,8,10,12,14))
title(main="(c)",adj=0, cex=0.75)
lines(tdinexp, nh4expmin, col="grey")
lines(tdinexp, nh4expmax, col="grey")
polygon(c(tdinexp, rev(tdinexp)), c(nh4expmax, rev(nh4expmin)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(tdinexp, nh4exp, type="p", col="red")
points(tdinexpout, nh4expmaxout, type='n')
lines(tdinexpout, nh4expminout, col="grey")
lines(tdinexpout, nh4expmaxout, col="grey")
polygon(c(tdinexpout, rev(tdinexpout)), c(nh4expmaxout, rev(nh4expminout)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(tdinexpout, nh4expout, type="p", col="blue")
abline(v=8, col="grey")
#abline(h=8, lty="dotted")
lines(tspan, out_axenic[,12], col="blue")
lines(tspan, out_bact[,12], col="red")
lines(tspan, out_axenic_noEX[,12], col="blue",lty="dotted")
lines(tspan, out_bact_noEX[,12], col="red",lty="dotted")



# Si plot
plot(tdinexp, siexpmax, type='n', 
     ylab = expression(paste("Si(OH) "[4]," [",mu,"mol L"^"-1","]")), 
     xlab = "time [days ]", ylim=c(0,29), xlim=c(0,16),
     axes=F)
axis(2, las=1, at=c(0,10,20))
axis(1, las=1, at=c(0,2,4,6,8,10,12,14))
title(main="(d)",adj=0, cex=0.75)
lines(tdinexp, siexpmin, col="grey")
lines(tdinexp, siexpmax, col="grey")
polygon(c(tdinexp, rev(tdinexp)), c(siexpmax, rev(siexpmin)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(tdinexp, siexp, type="p", col="red")
points(tdinexpout, siexpmaxout, type='n')
lines(tdinexpout, siexpminout, col="grey")
lines(tdinexpout, siexpmaxout, col="grey")
polygon(c(tdinexpout, rev(tdinexpout)), c(siexpmaxout, rev(siexpminout)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(tdinexpout, siexpout, type="p", col="blue")
abline(v=8, col="grey")
#abline(h=4.6, lty="dotted")
lines(tspan, out_axenic[,7], col="blue")
lines(tspan, out_bact[,7], col="red")
lines(tspan, out_axenic_noEX[,7], col="blue",lty="dotted")
lines(tspan, out_bact_noEX[,7], col="red",lty="dotted")
legend("topright", legend=c("BAC+ exp","BAC- exp", "EXT BAC+","EXT BAC-", "EXT -excr BAC+","EXT -excr BAC-"), 
       col=c("red","blue","red","blue","red","blue"), 
       lty=c(0,0,1,1,3,3), pch=c(1,1,NA,NA,NA,NA),
       box.lty=0, bg="transparent", cex=0.7)

#############################################

dev.off()


#############################################################################

#########
#minimizing (least sum of squares) (MANUAL)
# C
# texp: 1  2  3  4  5  6  7  8  9 10 11 14 15 16
# x0 <- c(C0, Chl0, N0, DIN0,
cmod <- out_axenic[c(1,101,201,301, 401, 501, 601, 701, 801, 901, 1001,  1301, 1401, 1501),2]
# Chl
chlmod <- out_axenic[c(1,101,201,301, 401, 501, 601, 701, 801, 901, 1001,  1301, 1401, 1501),3]
# N
nmod <- out_axenic[c(1,101,201,301, 401, 501, 601, 701, 801, 901, 1001,  1301, 1401, 1501),4]

# DIN
# t= 2  4  6  8  9 11 14 15
dinmod <- out_axenic[c(101,301, 501,  701, 801,  1001,  1301, 1401),5]+
  out_axenic[c(101,301, 501,  701, 801,  1001,  1301, 1401),12]

#Total
costTot <- sqrt( sum( (nmod-nexpout)^2/var(nexpout) ) +
                   sum( (chlmod - chlexpout)^2 / var(chlexpout)) +
                   sum( (cmod - cexpout)^2 / var(cexpout))+
                   sum( (dinmod -  DINexpout[-1])^2 / var(DINexpout[-1]))) 
costTot

##########################################################################

cmod <- out_bact[c(1,101,201,301, 401, 501, 601, 701, 801, 901, 1001,  1301, 1401, 1501),2]
costfc <- sqrt(sum((cexp - cmod)^2 ))
costfc <- sqrt(sum(scale(cexp - cmod)^2 ))
# Chl
chlmod <- out_bact[c(1,101,201,301, 401, 501, 601, 701, 801, 901, 1001,  1301, 1401, 1501),3]
costfchl <- sqrt( sum((chlexp - chlmod)^2 ) )
costfchl <- sqrt( sum(scale(chlexp - chlmod)^2 ) )
# N
nmod <- out_bact[c(1,101,201,301, 401, 501, 601, 701, 801, 901, 1001,  1301, 1401, 1501),4]
costfn <- sqrt( sum( (nexp - nmod)^2 ) )
costfn <- sqrt( sum(scale(nexp - nmod)^2 ) )
# DIN
# t= 2  4  6  8  9 11 14 15
dinmod <- out_bact[c(101,301, 501,  701, 801,  1001,  1301, 1401),5]+
  out_bact[c(101,301, 501,  701, 801,  1001,  1301, 1401),12]
costfdin <- sqrt( sum( (DINexp[-1] - dinmod)^2 ) )
costfdin <- sqrt( sum(scale(DINexp[-1] - dinmod)^2 ) )
#Total
costTot <- sqrt( sum( (nmod-nexp)^2/var(nexp) ) +
                   sum( (chlmod - chlexp)^2 / var(chlexp)) +
                   sum( (cmod - cexp)^2 / var(cexp))+
                   sum( (dinmod -  DINexp[-1])^2 / var(DINexp[-1]))) 
costTot




#####################################################
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











###################################################
# TUNING

##############################################################

plot(texp, cexp, col="black")
points(texp, nexp*10, col="red")
points(texp, chlexp*10, col="green")
points(tdinexp, DINexp/10, col="blue")

observed <- as.data.frame(cbind (c(rep("C",length(texp)),rep("Chl",length(texp)),rep("N",length(texp)),rep("DIN",length(tdinexp))),
                                 as.numeric(c(texp,texp,texp,tdinexp)),
                                 as.numeric(c(cexp,chlexp*10,nexp*10,(DINexp)/10))),stringsAsFactors = F)
colnames(observed)<-c("name","time","val")
#observed$err<- c(cerr, nerr, chlerr)  
observed$time<-as.numeric(observed$time)
observed$val<-as.numeric(observed$val)


####Cost function from FME

require(FME)
CostFUN <- function (pars) {
  out <- ode(x0, times = tspan, func = geider98v2,#_EPS_bact_131118, 
             parms = pars, method= "rk")
  modelled<-as.data.frame(cbind(out[,1],out[,2],out[,3]*10, out[,4]*10, out[,5]/10 +out[,12]/10))
  colnames(modelled)<-c("time","C","Chl","N","DIN")
  costFUN<-modCost(model = modelled, obs = observed, x="time", y = "val")
}
#x0 <- c(C0, Chl0, N0, DIN0, PSi0, Si0, bact0, 0, 0, don0, nh40) # 

CostFUN(pars)$model
par(mfrow=c(1,1))
plot(CostFUN(pars))#, ylim=c(-0.5,0.5))
abline(h=0)

Pars<-pars[c(11:21)]
SFUN<- sensFun(CostFUN, Pars)
summary(SFUN)
#plot(SFUN, which = c("C","N","Chl"), lwd=2)
#pairs(SFUN, which = c("C","N","Chl"), col=c("black","red","green"))

ident<-collin(SFUN)
ident[ident$collinearity>20 & ident$N == 2,]
head(ident, n = 20)
collin(SFUN, parset = c("n","kn"))

pars
Pars<-pars[c(15,16:19,21)]
SFUN<- sensFun(CostFUN, Pars)
ident<-collin(SFUN)
ident
ident[ident$collinearity>20 & ident$N == 2,]
#sens ks, Vmax, smin, nh4thresh is 0!?


lpars <- log(pars)
CostFUN2 <- function(lpars){
  CostFUN(c(exp(lpars)))
} #Add fixed (non-tuned) variables (otherwise it might not work unless they are loaded before globally)
Pars <- pars[c(15,16,19,20, 21)]
lowbound <- c   ( 10,0.1,   0.5,      0.1,  0.01)
highbound <- c  ( 100, 10,    9.3,      10,   0.5)

Fit <- modFit(f = CostFUN2, p = log(Pars), method="Pseudo",
              lower=log(lowbound), upper=log(highbound))
Pars2 <- exp(coef(Fit))
exp(coef(Fit))

Fit2 <- modFit(f = CostFUN2, p = log(Pars2), method="Nelder-Mead",
              lower=log(lowbound), upper=log(highbound))
Pars3 <- exp(coef(Fit2))
exp(coef(Fit2))