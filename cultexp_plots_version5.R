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
dinexp<-c(70,69.2,49.3,37.0,37.6,34.1,16.5,4.7,0.0)
dinexpmin<-c(38,30.1,46.0,29.9,34.5,33.4,13.9,0.0,0.0)
dinexpmax<-c(90,69.9,50.3,48.1,41.8,34.8,19.1,17.5,0.1)

nh4exp<-c(20,20.5,11.2,7.7,16.3,14.5,12.2,12.8,13.3)
nh4expmin<-c(19,11.8,8.5,7.3,15.8,13.8,12.2,12.7,13.1)
nh4expmax<-c(21,20.5,11.2,8.1,16.3,14.5,12.2,12.9,13.5)

po4exp<-c(29,18.2,14.2,16.3,27.9,38.4,33.9,47.0,11.3)
po4expmin<-c(16,16.8,14.0,15.1,16.6,36.7,31.7,2.1,4.1)
po4expmax<-c(52,18.5,14.7,16.6,43.7,40.2,36.1,53.3,18.5)

siexp<-c(24,20.8,10.9,4.5,2.2,0.5,0.6,1.8,0.5)
siexpmin<-c(13,9.2,10.0,3.8,1.9,0.3,0.5,0.3,0.4)
siexpmax<-c(30,21.4,11.6,5.0,4.0,0.7,0.8,2.2,0.7)

algexp<-c(6.9E+03,1.2E+04,1.5E+04,3.2E+04,3.8E+04,6.6E+04,8.7E+04,1.3E+05,
          1.3E+05,7.6E+04,9.2E+04,7.2E+04,8.9E+04,7.2E+04)
algexpmin<-c(6.0E+03,1.1E+04,1.5E+04,2.8E+04,3.6E+04,6.4E+04,8.4E+04,
             1.2E+05,1.2E+05,6.1E+04,8.0E+04,5.7E+04,8.3E+04,7.0E+04)
algexpmax<-c(7.3E+03,1.3E+04,1.9E+04,3.4E+04,4.3E+04,6.7E+04,9.1E+04,
             1.4E+05,1.3E+05,
             9.2E+04 + 0.3*9.2E+04,
             1.0E+05+ 0.3*1.0E+05,
             8.7E+04+ 0.3*8.7E+04,
             9.7E+04+ 0.3*9.7E+04,
             8.2E+04+ 0.3*8.2E+04)

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
bactexpout<-c(2.7E+03,1.2E+03,1.5E+03,6.7E+03,4.2E+03  ,4.7E+03,1.5E+04,3.8E+03,
              2.5E+04,1.4E+05,4.6E+04,3.0E+06)
bactexpminout<-c(1.3E+03,1.0E+03,1.5E+03,3.8E+03, 2.3E+03 ,2.7E+03,7.6E+03,3.5E+03,
                 1.6E+04,4.0E+04,2.9E+04,3.0E+06)
bactexpmaxout<-c(4.1E+03,2.9E+03,1.5E+03,4.2E+04,1.0E+04,   6.7E+03,2.1E+04,4.3E+03,
                 5.3E+04,2.3E+05,6.4E+04,3.0E+06)

# measured NO3 and NO2 (uM) 
tdinexpout<- c(0,2,4,6,8,9,11,14,15)
dinexpout<-c(70,49.1,46.8,26.9,2.5,0.5,0.3,0.2,0.7)
dinexpminout<-c(38,33.4,40.1,18.1,1.9,0.4,0.1,0.1,0.0)
dinexpmaxout<-c(90,63.7,47.7,28.5,3.4,0.7,1.0,0.3,1.9)
dinerrout<- abs(scale(dinexpmax - dinexpmin))

#measured nh4 (uM) (max >100 are removed, too high filtration pressure)
nh4expout<-c(20,4.4,4.3,2.6,5.8,3.25,2.1,15.2,8.0)
nh4expminout<-c(19,2.6,2.3,2.3,4.8,2.6,1.5,15.2,7.1)
nh4expmaxout<-c(21,7.7, 8.0, 3.3,  8.8, 12.1, 2.5, 15.2, 9.0)

# measured PO4 (uM)
po4expout<-c(29,15.2,13.0,13.2,17.5,14.5,14.5,13.1,15.3)
po4expminout<-c(16,15.0,12.4,12.9,13.8,13.4,14.4,12.8,12.2)
po4expmaxout<-c(32,16.3,13.3,14.6,24.4,14.8,14.5,13.4,21.3)

#measured Si(OH)4 in uM
siexpout<-c(24,15.1,13.3,3.8,0.8,0.7,0.7,0.7,0.6)
siexpminout<-c(13,11.1,11.5,2.9,0.1,0.2,0.7,0.6,0.3)
siexpmaxout<-c(30,20.3,13.8,5.6,1.1,1.0,1.5,0.9,0.9)




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



##################################################











pdf(file="Fig3_modpaper_rev3.pdf",width=12, height=6)

#tiff(filename="Fig3_modpaper_rev2.tiff", width=600, height=300,compression="none")
par(mar=c(5.1,5.1,4.1,2.1))

par(mfrow=c(1,2))
# algae cells plot
max(c(algexpmax,algexpmaxout))
options(scipen=1)
plot(texp, algexpmax/10e03, type='n', 
     ylab = expression(paste("Algae cells [10"^"4","cells mL"^"-1","]")), 
     xlab = "time [days]", ylim=c(0,15), xlim=c(0,16), axes=F)
axis(2, las=1, at=c(0,5,10,15))
axis(1, las=2, at=c(0,2,4,6,8,10,12,14,16))
title(main="(a)",adj=0, cex=0.75)
lines(texp, algexpmin/10e03, col="grey")
lines(texp, algexpmax/10e03, col="grey")
polygon(c(texp, rev(texp)), c(algexpmax/10e03, rev(algexpmin/10e03)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(texp, algexp/10e03, type="p", bg="red", pch=21, cex=1.2)
points(texpout, algexpmaxout/10e03, type='n')
lines(texpout, algexpminout/10e03, col="grey")
lines(texpout, algexpmaxout/10e03, col="grey")
polygon(c(texpout, rev(texpout)), c(algexpmaxout/10e03, rev(algexpminout/10e03)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(texpout, algexpout/10e03, type="p", bg="blue", pch=21, cex=1.2)

# bacteria cells plot
max(c(bactexpmax,bactexpmaxout))
plot(tbactexp, bactexpmax/10e06, type='n', 
     ylab = expression(paste("Bacteria [10"^"7","cells mL"^"-1","]")), 
     xlab = "time [days]", ylim=c(0,7),xlim=c(0,14), axes=F)
axis(2, las=1, at=c(0,1,2,3,4,5,6,7))
axis(1, las=2, at=c(0,2,4,6,8,10,12,14))
title(main="(b)",adj=0, cex=0.75)
lines(tbactexp[-12], bactexpmin[-12]/10e06, col="grey")
lines(tbactexp[-12], bactexpmax[-12]/10e06, col="grey")
polygon(c(tbactexp[-12], rev(tbactexp[-12])), c(bactexpmax[-12]/10e06, rev(bactexpmin[-12]/10e06)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(tbactexp[-12], bactexp[-12]/10e06, type="p", bg="red", pch=21, cex=1.2)
points(tbactexpout[-12], bactexpmaxout[-12]/10e06, type='n')
lines(tbactexpout[-12], bactexpminout[-12]/10e06, col="grey")
lines(tbactexpout[-12], bactexpmaxout[-12]/10e06, col="grey")
polygon(c(tbactexpout[-12], rev(tbactexpout[-12])), c(bactexpmaxout[-12]/10e06, rev(bactexpminout[-12]/10e06)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(tbactexpout[-12], bactexpout[-12]/10e06, type="p", bg="blue", pch=21, cex=1.2)
points(tbactexp[12],bactexp[12]/10e06, pch=8, col="red")
points(tbactexp[12],bactexpout[12]/10e06, pch=8, col="blue")
abline(h=0)

legend("topleft", legend=c("BAC+ exp","BAC- exp"), 
       pt.bg=c("red","blue"), 
       lty=c(0,0), pch=c(21,21),
       box.lty=0, bg="transparent")

dev.off()
















##################################################################################
pdf(file="Fig2_modpaper_rev3.pdf")#, width=600)
#tiff(filename="Fig2_modpaper_rev2.tiff", width=600, compression="none")
par(mar=c(5.1,5.1,4.1,2.1))


par(mfrow=c(2,2))

# NOx cells plot
plot(tdinexp, dinexpmax, type='n', 
     ylab = expression(paste("NO "[X]," [",mu,"mol L"^"-1","]")), 
     xlab = "time [days]", ylim=c(0,100), xlim=c(0,16),
     axes=F)
axis(2, las=1, at=c(0,20,40,60,80,100))
axis(1, las=2, at=c(0,2,4,6,8,10,12,14,16))
title(main="(a)",adj=0, cex=0.75)
lines(tdinexp, dinexpmin, col="grey")
lines(tdinexp, dinexpmax, col="grey")
polygon(c(tdinexp, rev(tdinexp)), c(dinexpmax, rev(dinexpmin)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(tdinexp, dinexp, type="p", bg="red", pch=21)
points(tdinexpout, dinexpmaxout, type='n')
lines(tdinexpout, dinexpminout, col="grey")
lines(tdinexpout, dinexpmaxout, col="grey")
polygon(c(tdinexpout, rev(tdinexpout)), c(dinexpmaxout, rev(dinexpminout)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(tdinexpout, dinexpout, type="p", bg="blue", pch=21, cex=1.2)
abline(v=8, col="grey")
#abline(h=0.5, lty="dotted")
#lines(tspan, out_axenic[,5], col="blue")
#lines(tspan, out_bact[,5], col="red")
#lines(tspan, out_axenic_noEX[,5], col="blue",lty="dotted")
#lines(tspan, out_bact_noEX[,5], col="red",lty="dotted")


# NH4  plot
plot(tnh4exp, nh4expmaxFull, type='n', 
     ylab = expression(paste("NH "[4]," [",mu,"mol L"^"-1","]")), 
     xlab = "time [days]", ylim=c(0,25), xlim=c(0,16), 
     axes=F)
axis(2, las=1, at=c(0,5,10,15,20,25))
axis(1, las=2, at=c(0,2,4,6,8,10,12,14,16))
title(main="(b)",adj=0, cex=0.75)
lines(tnh4exp, nh4expminFull, col="grey")
lines(tnh4exp, nh4expmaxFull, col="grey")
polygon(c(tnh4exp, rev(tnh4exp)), c(nh4expmaxFull, rev(nh4expminFull)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(tnh4exp, nh4expFull, type="p", bg="red", pch=21)
points(tnh4expout, nh4expmaxoutFull, type='n')
lines(tnh4expout, nh4expminoutFull, col="grey")
lines(tnh4expout, nh4expmaxoutFull, col="grey")
polygon(c(tnh4expout, rev(tnh4expout)), c(nh4expmaxoutFull, rev(nh4expminoutFull)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(tnh4expout, nh4expoutFull, type="p", bg="blue", pch=21)
abline(v=8, col="grey")
#abline(h=8, lty="dotted")
#lines(tspan, out_axenic[,12], col="blue")
#lines(tspan, out_bact[,12], col="red")
#lines(tspan, out_axenic_noEX[,12], col="blue",lty="dotted")
#lines(tspan, out_bact_noEX[,12], col="red",lty="dotted")

#po4outlier<-po4expmin[8]
#po4expmin[8]<-(po4exp[8])
#po4exp[8]<-mean(po4exp[8],po4expmax[8])


# PO4 plot
max(c(po4expmax,po4expmaxout))
plot(tdinexp, po4expmax, type='n', 
     ylab = expression(paste("PO "[4]," [",mu,"mol L"^"-1","]")), 
     xlab = "time [days]", xlim=c(0,16),ylim=c(0,54), axes=F)
axis(2, las=1, at=c(0,10,20,30,40,50))
axis(1, las=2, at=c(0,2,4,6,8,10,12,14,16))
title(main="(c)",adj=0, cex=0.75)
lines(tdinexp, po4expmin, col="grey")
lines(tdinexp, po4expmax, col="grey")
polygon(c(tdinexp, rev(tdinexp)), c(po4expmax, rev(po4expmin)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(tdinexp, po4exp, type="p", bg="red", pch=21)
points(tdinexpout, po4expmaxout, type='n')
lines(tdinexpout, po4expminout, col="grey")
lines(tdinexpout, po4expmaxout, col="grey")
polygon(c(tdinexpout, rev(tdinexpout)), c(po4expmaxout, rev(po4expminout)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(tdinexpout, po4expout, type="p", bg="blue", pch=21)
points(14,po4outlier, pch=8, col="red")
abline(v=8, col="grey")
#abline(h=8, lty="dotted")

# Si plot
plot(tdinexp, siexpmax, type='n', 
     ylab = expression(paste("Si(OH) "[4]," [",mu,"mol L"^"-1","]")),
     xlab = "time [days]", ylim=c(0,30), xlim=c(0,16),
     axes=F)
axis(2, las=1, at=c(0,10,20,30))
axis(1, las=2, at=c(0,2,4,6,8,10,12,14,16))
title(main="(d)",adj=0, cex=0.75)
lines(tdinexp, siexpmin, col="grey")
lines(tdinexp, siexpmax, col="grey")
polygon(c(tdinexp, rev(tdinexp)), c(siexpmax, rev(siexpmin)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(tdinexp, siexp, type="p", bg="red", pch=21)
points(tdinexpout, siexpmaxout, type='n')
lines(tdinexpout, siexpminout, col="grey")
lines(tdinexpout, siexpmaxout, col="grey")
polygon(c(tdinexpout, rev(tdinexpout)), c(siexpmaxout, rev(siexpminout)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(tdinexpout, siexpout, type="p", bg="blue", pch=21)
abline(v=8, col="grey")
#abline(h=4.6, lty="dotted")
#lines(tspan, out_axenic[,7], col="blue")
#lines(tspan, out_bact[,7], col="red")
#lines(tspan, out_axenic_noEX[,7], col="blue",lty="dotted")
#lines(tspan, out_bact_noEX[,7], col="red",lty="dotted")

legend("topright", legend=c("BAC+ exp","BAC- exp"), 
       pt.bg=c("red","blue"), 
       lty=c(0,0), pch=c(21,21),
       box.lty=0, bg="transparent")

dev.off()




#######################################################################

























##############################################################
# POC plot

pdf(file="Fig4_modpap_rev3.pdf")# width=600)
#tiff(filename="Fig4_modpaper_rev2.tiff", width=600, compression="none")

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))

plot(texp, cexpmax, type='n', 
     ylab = expression(paste("POC [",mu,"g mL"^"-1","]")), 
     xlab = "time [days]",ylim=c(0,10), xlim=c(0,16),
     axes=F)
axis(2, las=1, at=c(0,2,4,6,8,10))
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
#lines(tspan, out_axenic[,2], col="blue")
#lines(tspan, out_bact[,2], col="red")

legend("topleft", legend=c("BAC+ exp","BAC- exp"), 
       pt.bg=c("red","blue"), 
       lty=c(0,0), pch=c(21,21),
       box.lty=0, bg="transparent")

# PON plot
plot(texp, nexpmax, type='n', 
     ylab = expression(paste("PON [",mu,"g mL"^"-1","]")), 
     xlab = "time [days]", xlim=c(0,16),ylim=c(0,2),
     axes=F)
axis(2, las=1, at=c(0,0.5,1,1.5,2))
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
#lines(tspan, out_axenic[,4], col="blue")
#lines(tspan, out_bact[,4], col="red")



# POC:PON plot

plot(texp, cexpmax/nexpmin, type='n', 
     ylab = expression(paste("C : N ratio [",mu,"g C : ",mu,"g N]")), 
     xlab = "time [days]",ylim=c(0,25), xlim=c(0,16),
     axes=F)
axis(2, las=1, at=c(0,5,10,15,20,25))
axis(1, las=2, at=c(0,2,4,6,8,10,12,14,16))
lines(texp, cexpmin/nexpmax, col="grey")
lines(texp, cexpmax/nexpmin, col="grey")
polygon(c(texp, rev(texp)), c(cexpmax/nexpmin, rev(cexpmin/nexpmax)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(texp, cexp/nexp, type="p", bg="red", pch=21)
points(texpout, cexpmaxout/nexpminout, type='n')
lines(texpout, cexpminout/nexpmaxout, col="grey")
lines(texpout, cexpmaxout/nexpminout, col="grey")
polygon(c(texpout, rev(texpout)), c(cexpmaxout/nexpminout, rev(cexpminout/nexpmaxout)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(texpout, cexpout/nexpout, type="p", bg="blue", pch=21)
title(main="(c)",adj=0, cex=0.75)
abline(v=8, col="grey")
#lines(tspan, out_axenic[,2]/out_axenic[,4], col="blue")
#lines(tspan, out_bact[,2]/out_bact[,4], col="red")

#chloutlier<-chlexpmin[8]
#chlexpmin[8]<-chlexp[8]
#chlexp[8]<-mean(chlexp[8], chlexpmax[8])

# Chl plot
plot(texp, chlexpmax, type='n', 
     ylab = expression(paste("Chlorophyll a [",mu,"g mL"^"-1","]")), 
     xlab = "time [days]", xlim=c(0,16),ylim=c(0,0.8),
     axes=F)
axis(2, las=1, at=c(0,0.2,0.4,0.6,0.8))
axis(1, las=2, at=c(0,2,4,6,8,10,12,14,16))
lines(texp, chlexpmin, col="grey")
lines(texp, chlexpmax, col="grey")
polygon(c(texp, rev(texp)), c(chlexpmax, rev(chlexpmin)),
        col=rgb(1, 0, 0, 0.2), border=NA)
lines(texp, chlexp, type="p", bg="red", pch=21)
points(texpout, chlexpmaxout, type='n')
lines(texpout, chlexpminout, col="grey")
lines(texpout, chlexpmaxout, col="grey")
polygon(c(texpout, rev(texpout)), c(chlexpmaxout, rev(chlexpminout)),
        col=rgb(0, 0, 1, 0.2), border=NA)
lines(texpout, chlexpout, type="p", bg="blue", pch=21)
title(main="(d)",adj=0, cex=0.75)
points(8,chloutlier, pch=8, col="red")
abline(v=8, col="grey")
#lines(tspan, out_axenic[,3], col="blue")
#lines(tspan, out_bact[,3], col="red")

dev.off()
