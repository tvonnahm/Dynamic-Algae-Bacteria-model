geider98old <- function(t,x,paramMod){
  with(as.list(c(x, paramMod)),{

# Function for modelling Photosynthesis, Chlorophyll synthesis, and Nitrogen uptake
# in a single diatom species under varying DIN and Si concentrations
# and remineralisation of refractory excreted nitrogen and more recalcitrant background DON

# Based on:
# Geider 1998: Baseline model (C:N.Chl pool based)
# Spilling et al. 2009: Silicate uptake after Michaelis-menton kinetics
# Werner 1987: Response of Photosynthesis and Chl synthesis under Si limitation
# Remineralisation of excreted and background DON when C:N <10 (Tezuka 1990)
# Remineralisation as function of bacteria (bact growth as logistic function)    
    
# Additional ideas:
# Model NH4 and NO3 uptake with different kinetics
# Model Algae cell numbers (if Si not limiting & C:N ratio below threshhole -> C*CF = cells)

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
    #Pcref <-  0.9#0.9#0.9 # (TUNE) 1.5/3.5
    # cost - Cost of biosynthesis (parameter)
    #cost <- 1 # (TUNE) #gC / gN 2
    # Vcn  - max N uptake (calculated)
    # Rc   - carbon based maintenance metabolic rate constant (parameter)
    #Rc <- 0.02 #(TUNE) 0.025
    #xf <- 0.04 (excretion rate of C and N)
    # Rchl and RN - Chl and N based maintenance metabolic rates(calculated)
    # pchl - carbon specific light limited photosynt. rate (calculated)
    # thetac - Chl : C ratio in the cell (calculated)
    # thetan - Chl : N ratio in the cell (Calculated)
    # thetanmax - max thetan (parameter)
    #thetanmax <- 1.7 # (measured) #g chl a / g N
    # Q    - N : C ratio in the cell (calculated)
    # Qmin - minimum value of Q (parameter)
    #Qmin <- 0.05 # 0.06 (0.06 measured) g n / g C 0.04
    # Qmax - maximum value of Q (parameter)
    #Qmax <- 0.35 #(0.025 measured) 0.17/024
    # Alfchl - Chl specific linitial C assimilation rate (parameter)
    #Alfchl <- 0.1 # 0.075 (TUNE) # for Pavlova lutheri 1
    # I   - incident scalar irradiance (variable, model input)
    #I <- 100 # (measured)
    # my  - specific growth rate (calculated)
    # DIN - inorganic n in medium (variable, model input) in uM
    #DIN <- 35
    # n   - shape factor describing dependence of Vcmax on Q
    #n <- 3.9#4 #(TUNE) 1
    # kn  - half saturation constant for nitrogen uptake (parameter)
    #kn <- 1 #(TUNE) 2/3
    #Qsmax <-0.13
    #Qsmin <-0.0037
    #ns <-1
    #Rs <- Rc
    #Vmax <- 0.33 #0.32-0.9 # max Si uptake rate
    #ks <-10 #0.5-7.6 # half saturation constant Si uptake
    #smin <- 1.6 #1.5-6 #min dis Si allowing Si assimilation
    #rem <- 0#0.8 # remineralisation rate of excreted don
    #remd <- 0#0.7#remineralisationr ate of dissolved organic nitrogen backgroun

#Input variables
C <- x[1]
Chl <- x[2]
N <- x[3]
DIN <-x[4]
thetac <- Chl / C
Q <- N / C
thetan <- Chl / N
Rchl <- Rc
Rn <- Rc



#Photosynthesis (Geider 1998)
Pcmax <- Pcref * ((Q - Qmin)/(Qmax - Qmin))
Ik <- Pcmax / (Alfchl * thetac)
Pc <- Pcmax * (1 - exp(-I/Ik))

#max. N uptake (Geider 1998)
Vcref <- Pcref * Qmax
Vcn <- Vcref * ((Qmax - Q)/(Qmax - Qmin))^n * DIN / (DIN + kn)


#Carbon Synthesis (Geider 1998)
# apparent Photosynthesis is reduced by 80% under Si limitation (Werner 1987)
# Carbon excretion (Function of Carbon content)

  dcdt <- (Pc - cost * Vcn -Rc - xf) * C 
  my <- (Pc - cost * Vcn -Rc)

#Chl Synthesis (Geider 1998)
# Chl synthesis is inhibited when Si is limitng after 5 h (Werner 1987)
  pchl <- thetanmax * (Pc / (Alfchl * thetac * I))
 dchldt <- ((pchl * Vcn)/thetac - Rchl) * Chl 

#Nitrogen uptake (Geider 1998)
# PLUS Nitrgen excretion
dndt <- (Vcn / Q - Rn - xf) * N 


#from eq 1 and 2
#dQdt <- Vcn - my * Q
#dthetaCdt <- Vcn * pchl - thetac * my

if (is.na(dcdt)){ dcdt <-0 }else {dcdt <- dcdt }
if (is.na(dchldt)){ dchldt <-0 }else {dchldt <- dchldt }
if (is.na(dndt)){ dndt <-0 }else {dndt <- dndt}


# DIN consumption (Geider 1998) and production by remineralisation


ddindt <- (- ((Vcn / Q)*N))/14 * 1e3 # minus bacterial N uptake after Redfield




if (is.na(ddindt)){ ddindt <-0 }else {ddindt <- ddindt}

list(c(dcdt,dchldt, dndt, ddindt))

    
  }) 
}

return