##############################################################
# POC plot

pdf(file="Fig4_modpap_rev2.pdf")

par(mar=c(5.1,5.1,4.1,2.1))
par(mfrow=c(2,2))

plot(texp, cexpmax, type='n', 
     ylab = expression(paste("POC [",mu,"g mL"^"-1","]")), 
     xlab = "time [days ]",ylim=c(0,10), xlim=c(0,16),
     axes=F)
axis(2, las=1, at=c(0,2,4,6,8,10))
axis(1, las=2, at=c(0,2,4,6,8,10,12,14))
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
#lines(tspan, out_axenic[,2], col="blue")
#lines(tspan, out_bact[,2], col="red")

legend("topleft", legend=c("BAC+ exp","BAC- exp"), 
       col=c("red","blue"), 
       lty=c(0,0), pch=c(1,1),
       box.lty=0, bg="transparent")

# PON plot
plot(texp, nexpmax, type='n', 
     ylab = expression(paste("PON [",mu,"g mL"^"-1","]")), 
     xlab = "time [days ]", xlim=c(0,16),ylim=c(0,2),
     axes=F)
axis(2, las=1, at=c(0,0.5,1,1.5,2))
axis(1, las=2, at=c(0,2,4,6,8,10,12,14))
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
title(main="b)",adj=0, cex=0.75)
abline(v=8, col="grey")
#lines(tspan, out_axenic[,4], col="blue")
#lines(tspan, out_bact[,4], col="red")



# POC:PON plot

plot(texp, cexpmax/nexpmin, type='n', 
     ylab = expression(paste("C : N ratio [",mu,"g C : ",mu,"g N]")), 
     xlab = "time [days ]",ylim=c(0,25), xlim=c(0,16),
     axes=F)
axis(2, las=1, at=c(0,5,10,15,20))
axis(1, las=2, at=c(0,2,4,6,8,10,12,14))
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
title(main="c)",adj=0, cex=0.75)
abline(v=8, col="grey")
#lines(tspan, out_axenic[,2]/out_axenic[,4], col="blue")
#lines(tspan, out_bact[,2]/out_bact[,4], col="red")


# Chl plot
plot(texp, chlexpmax, type='n', 
     ylab = expression(paste("Chlorophyll a [",mu,"g mL"^"-1","]")), 
     xlab = "time [days ]", xlim=c(0,16),ylim=c(0,0.7),
     axes=F)
axis(2, las=1, at=c(0,0.2,0.4,0.6))
axis(1, las=2, at=c(0,2,4,6,8,10,12,14))
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
title(main="d)",adj=0, cex=0.75)
abline(v=8, col="grey")
#lines(tspan, out_axenic[,3], col="blue")
#lines(tspan, out_bact[,3], col="red")

dev.off()
