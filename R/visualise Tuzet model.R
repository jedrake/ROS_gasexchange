
# Following code allows to visualise the solution to the Tuzet model
# It is the intersection of two relationships between gs and psil
# - supply curve E = kh (psil - psis) = gs D
# - demand curve gs = g1 A f(psil) / Ca  

# The f(psi-L) function used in the Tuzet-Leuning model
fpsilfun <- function(psil, psiv, sf){
  (1 + exp(sf*psiv) )/ (1 + exp(sf*(psiv-psil)))
}


# range of psil values
psil <- seq(-5,-0.1,by=0.1)

# supply side stomatal conductance: 
# E = gs*D = kh*(psil-psis) so gs = kh*(psil-psis)/D
# note, this assumes that kh is constant, not dependent on psil
# could also try kh = f(psil) or Sperry E = integral of kh from psis to psil
gsE <- function(psil, D = 1.5, kh = 0.3, psis = -0.1) {
  return (kh*(psis-psil)/D)
}

# demand side stomatal conductance: solving the following three equations
# A = Vcmax *(Ci-Gamma)/(Ci+Km) - Rday
# Ci = Ca - 1.6 A / gs
# gs = g0 + g1*f(psi-L)*A/Ca
# yields quadratic
# Code from Maestra for this quadratic
# GSDIVA = G1 / (CS - GAMMA) / (1 + VPD/D0L) * FSOIL
# A = G0 + GSDIVA * (VCMAX - RD)
# B = (1. - CS*GSDIVA) * (VCMAX - RD) + G0 * (KM - CS)
# +      - GSDIVA * (VCMAX*GAMMASTAR + KM*RD)
# C = -(1. - CS*GSDIVA) * (VCMAX*GAMMASTAR + KM*RD) - G0*KM*CS
# CIC = QUADP(A,B,C,IQERROR)
# Note: g0 > 0 for stability, can take 0.001 but must be > 0
gsA <- function(psil, Vcmax=90, Jmax = 140, g0 = 0.001, g1 = 4.0, Ca = 400,
                sf = 25, psiv = -1.5, gamma = 42, Km = 700, Rd = 0) {
  f <- fpsilfun(psil,psiv,sf)
  gsdiva <- g1 * f / Ca / 1.6
  g0c <- g0/1.6
  A <- g0c + gsdiva*(Vcmax-Rd)
  B <- (1 - Ca*gsdiva)*(Vcmax-Rd) + g0c*(Km-Ca) - gsdiva*(Vcmax*gamma + Km*Rd)
  C <- -(1-Ca*gsdiva) * (Vcmax*gamma + Km*Rd) - g0c*Km*Ca
  Ci <- (-B + sqrt(B^2-4*A*C))/(2*A)
  Photo <- Vcmax*(Ci-gamma)/(Ci+Km) - Rd
  Cond <- 1.6*Photo/(Ca-Ci)
  ret <- data.frame(Photo,Cond,Ci)
  return(ret)
}

# The simplified demand curve, with no g0 or Rd ie Ci fixed
gsAsimp <- function(psil, Vcmax=90, Jmax = 140, g1 = 4.0, Ca = 400,
                  sf = 25, psiv = -1.5, gamma = 42, Km = 700) {
  f <- fpsilfun(psil,psiv,sf)
  Ci <- Ca * (1-1.6/g1/f)
  Photo <- Vcmax*(Ci-gamma)/(Ci+Km)
  Cond <- g1*f*Photo/Ca
#  ret <- (g1*f*(Ca - gamma)-1.6*Ca) / (g1*f*(Ca + Km)-1.6*Ca)
  ret <- data.frame(Photo,Cond,Ci)
  return(ret)
}

# range of psil values
psil <- seq(-5,-0.1,by=0.1)

# visualise demand curves - see the difference g0 makes
plot(psil,gsA(psil)$Cond,type='l')
points(psil,gsAsimp(psil)$Cond,type='l',col="red")

# Visualise full model - solution is given by intersection of curves
plot(psil,gsE(psil,D=1),type='l',ylim=c(0,0.5))  # supply curve
points(psil,gsE(psil,D=2),type='l',col="red")   # increases D
points(psil,gsE(psil,D=2,psis=-0.8),type='l',col="brown")  # reduces SWC
points(psil,gsA(psil,sf=9.6, psiv=-1)$Cond,type='l',col="green")  # demand curve, Pira
points(psil,gsA(psil,sf=2.5, psiv=-1)$Cond,type='l',col="blue")   # varies sf (Eusi)
points(psil,gsA(psil,sf=1.6, psiv=-2.5)$Cond,type='l',col="purple")  # varies psiv

# Find max possible values for gs, A and Ci given G1, Vcmax
g1fixed <- 10
Vfixed <- 80
maxCi <- 400*(1-1.6/g1fixed)
maxA <- Vfixed*(maxCi-42)/(maxCi+700)
maxgs <- g1fixed*maxA/400
