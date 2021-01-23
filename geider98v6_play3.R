########################################################
########################################################
### Modelling script/function (extended model)#########
### v2.0 by Tobias Reiner Vonnahme #####################
########################################################

# Function for modelling Photosynthesis, Chlorophyll synthesis, and Nitrogen uptake
# in a single diatom species under varying DIN and Si concentrations
# and remineralisation of refractory excreted nitrogen and more recalcitrant background DON

# Based on:
# Geider 1998: Baseline model (C:N.Chl pool based)
# Spilling et al. 2009: Silicate uptake after Michaelis-Mnton kinetics
# Werner 1987: Response of Photosynthesis and Chl synthesis under Si limitation
# Remineralisation of excreted and background DON when C:N <10 (Tezuka 1990)
# Remineralisation as function of bacteria (bact growth as logistic function)    

# input x <- [C, Chl, N, DIN, pSi, dSi, 0, 0]
# output - rates of change in C Chl and N and DIN pSi, dSi, VSi and Vcn

#parameters and variables
#___________________________________________________________
# C   - Carbon content (starting condition, modelled)
# Chl - Chlorophyll content (starting condition, modelled)
# N   - Nitrogen content (starting condition, modelled)
# pc  - C-specific rate of photosynthesis (calculated)
# Pcmax - maximum value of pc at given temperature (calculated)
# Pcref - Value of PCmax at a given reference temperature (calculated)
# cost - Cost of biosynthesis (parameter)
# Vcn  - max N uptake (calculated)
# Rc   - carbon based maintenance metabolic rate constant (parameter)
#xf  -  excretion rate of C and N)
# Rchl and RN - Chl and N based maintenance metabolic rates(calculated)
# pchl - carbon specific light limited photosynt. rate (calculated)
# thetac - Chl : C ratio in the cell (calculated)
# thetan - Chl : N ratio in the cell (Calculated)
# thetanmax - max thetan (parameter)
#thetanmax <- 1.7 # (measured) #g chl a / g N
# Q    - N : C ratio in the cell (calculated)
# Qmin - minimum value of Q (parameter)
# Qmax - maximum value of Q (parameter)
# Alfchl - Chl specific linitial C assimilation rate (parameter)
# I   - incident scalar irradiance (variable, model input)
# my  - specific growth rate (calculated)
# DIN - inorganic n in medium (variable, model input) in uM
# n   - shape factor describing dependence of Vcmax on Q
# kn  - half saturation constant for nitrogen uptake (parameter)
#Vmax - max Si uptake rate
#ks - half saturation constant Si uptake
#smin - min dis Si allowing Si assimilation
#rem - remineralisation rate of excreted don
#remd - remineralisationr ate of dissolved organic nitrogen backgroun

geider98v2 <- function(t,x,paramMod){
  with(as.list(c(x, paramMod)),{
    
    #Input variables
    Rc <- Rc - xf
    
    C <- x[1]
    Chl <- x[2]
    N <- x[3]
    DIN <-x[4]
    S <- x[5]
    ds <- x[6]
    thetac <- Chl / C
    Q <- N / C
    Qs <- S / C
    thetan <- Chl / N
    Rchl <- Rc
    Rn <- Rc
    bact <- x[7]
    donr <- x[10]
    nh4 <- x[11]
    donl <- x[12]
   
    
    #Bacteria growth (logistic growth) -> May be changed to Michaelis_Menton kinetics 
    #dbactdt <- bact*(mubact)*(bactmax-bact)
    dbactdt <- bact*(mubact)*(1-bact/bactmax)
    
    #Photosynthesis (Geider et al 1998)
    if (Q < Qmin) {Q <- Qmin} else if (Q > Qmax) {Q <- Qmax} else Q <- Q 
    Pcmax <- Pcref * ((Q - Qmin)/(Qmax - Qmin))
    Ik <- Pcmax / (Alfchl * thetac)
    Pc <- Pcmax * (1 - exp(-I/Ik))
    
    #max. N uptake (Geider et al 1998)
    Vcref <- Pcref * Qmax
    Vcn1 <- Vcref * ((Qmax - Q)/(Qmax - Qmin))^n * DIN / (DIN + kn)
    Vcn2 <- (0.05 * 0.2 * Q) * 0.42/0.005 * nh4 / (nh4 + kn2) # SHANIM Flynn et al 1997
    if (is.na(nh4)){nh4<-0}
    if (nh4 > nh4thres) {
      Vcn <- Vcn2 + 0.2 * Vcn1
    } else {
      Vcn <- Vcn1 + Vcn2
    }
    
    # max. S uptake
    # Spilling et al. 2009 (Si uptake in diatoms after Michaelis menton)
    V <- Vmax * (ds - smin)/(ks-smin) 
    Vcs <- V * ds * C
    
    #Carbon Synthesis (Geider et al. 1998)
    # apparent Photosynthesis is reduced under Si limitation (Werner 1987)
    # Carbon excretion (Function of Carbon content)
    
    ds[!is.finite(ds)] <- 0
    if (is.na(ds) | is.null(ds) | ds < smin *2){  
      dcdt <- SiPS*((Pc - cost * Vcn -Rc - xf) * C )
      my <- SiPS*((Pc - cost * Vcn -Rc))
      rem <- rem2
    }else {
      dcdt <- (Pc - cost * Vcn -Rc - xf) * C 
      my <- (Pc - cost * Vcn -Rc -xf)
      xf <- xf
    }
    
    #Chl Synthesis (Geider 1998)
    # Chl synthesis is inhibited when Si is limitng after 5 h (Werner 1987)
    pchl <- thetanmax * (Pc / (Alfchl * thetac * I))
    if (is.na(ds) | is.null(ds) | ds < smin*2){
      dchldt <- -Rchl * Chl
    }else {
      dchldt <- ((pchl * Vcn)/thetac - Rchl) * Chl 
    }
    
    #Nitrogen uptake (Geider 1998)
    # PLUS Nitrogen excretion
    dndt <- (Vcn / Q - Rn - xf) * N 
    
    # Si uptake
    dsdt <- Vcs
    ddsdt <- - Vcs 
    if (is.na(dcdt)){ dcdt <-0 }else {dcdt <- dcdt }
    if (is.na(dchldt)){ dchldt <-0 }else {dchldt <- dchldt }
    if (is.na(dndt)){ dndt <-0 }else {dndt <- dndt}
    
    # DIN consumption (Geider et al. 1998) and production by remineralisation
    # Remineralisation of refractory freshly excreted DON (xf rem N) and old DON (don remd)
    # Remineralisation only if C:N <10 (Tezuka 1990) and after lag pahse for this experiment (4 days)
    
    #bacteria N demand
    Bacn <- dbactdt/6.625
    
    if (Q > 0.1 & t > 4){ # 
      dnh4dt <- (-(Vcn2 / Q)* N +
                   bact * donr * remd +
                   bact * donl * rem)/14 * 1e3 # minus bacterial N uptake after Redfield
      ddonrdt <- (- bact * donr * remd)/14 * 1e3
      ddonldt <- (xf * N - bact * donl * rem - Bacn) /14 * 1e3
    }else if (nh4 < nh4thres){
      dnh4dt <- (-(Vcn2 / Q)* N)/14 * 1e3 # minus bacterial N uptake after Redfield
      ddonldt <- (xf * N)/14 * 1e3
      ddonrdt <- 0
      } else {
      dnh4dt <- (-(Vcn2 / Q)* N -
                   Bacn)/14 * 1e3 # minus bacterial N uptake after Redfield
      ddonldt <- (xf * N)/14 * 1e3
      ddonrdt <- 0
    }
    #IDEA : when rem happening: BacN from DON
    # when NH4 over threhs: Baxcn From NH4
    # when NH4 belwo thresh: Bacn from NO3
    
    if (nh4 < nh4thres){
      ddindt <- (-(Vcn1 / Q)*N - 
                   Bacn)/14 * 1e3 # minus bacterial N uptake after Redfield
    } else {
      ddindt <- (- 0.2*((Vcn1 / Q)*N))/14 * 1e3 # no bacteria uptake since NH4 is abundant
    }
    
    dcdt[!is.finite(dcdt)] <- 0
    dndt[!is.finite(dndt)] <- 0
    if (is.na(ddindt)){ ddindt <-0 }else {ddindt <- ddindt}
    if (ddindt > DIN){ddindt > DIN} else {ddondt <- ddindt}
    
    #depsndt <-(donr + donl) * xeps
    #depscdt <- depsndt/14*12*6.625
    #ddonrdt <- ddonrdt - donr/donl * xeps
    #ddonldt <- ddonldt - (1-(donr/donl)) * xeps
    
    
    list(c(dcdt,dchldt, dndt, ddindt, dsdt, ddsdt, dbactdt, V, Vcn, ddonrdt, dnh4dt, ddonldt))
    
  }) 
}

return