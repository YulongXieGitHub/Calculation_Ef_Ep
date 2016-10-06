# Ef_Ep.R
#
# Based on the formula in Fei Gao's manuscript "GaoFei_manuscripts.pdf"
#
# 1. unit conversion
J_2_eV <- 6.24150934e18			# Joel (kg*m2/s2) to eV

# 2. 
C.Avogadro <- 6.0221413e23		# Avogadro number
m.elec     <- 9.1093829e-31		# electron mass in kg
e.elec     <- 1.60217657e-19		# electron charge in Coulombs
v.light    <- 2.997925e8		# light speed in m/s
h          <- 6.62606957e-34		# Planck constant in m2 kg / s
h_hat      <- h / (2 * pi)		# Planck constant in m2 kg / s divided by 2 pi
h_hat.eVs  <- h_hat * J_2_eV		# Planck constant in eVs divided by 2 pi
mc2.J      <- m.elec * v.light^2	# mc2 in J (kg*m2/s2)
mc2.eV     <- mc2.J * J_2_eV		# mc2 in eV


# @mat_array	= ('SrI2','CaF2', 'BaF2','NaI', 'LaBr3','CsI', 'CZT',           'SiGe');	# the Material names
# @bandGap_array= ( 4.5,   11.8,   11.0,   5.9,   5.6,    6.1,   1.57,           0.945);	# band gap
# @fermi_array	= ( 9.27,  14.842, 15.17,  8.762,12.55,   6.97, 17.88534,       11.75177);	# Fermi energy
# @plasma_array	= (13.293, 18.9202,19.24, 12.742,16.68,  10.73, 21.76116,       15.88132);	# plasma energy
# @Zi_array	= (16,     16,     16,     8,    24,      8,    10,              8);		# Z_i (number of electron in the outmost shells of the molecule
# @cut_array	= ( 27.4,  78.5,   70.5,  33.9,  35.0,   36.0,  30.6,    21.6);			# this cut energy should be slightly larger than the maximum energies found from plasmon_findmin.pl (if much larger, should cause no problem with only impact on energy grid)
# @a0		= (  0.001, 0.001,  0.001, 0.001, 0.001,  0.001, 0.001,   0.001);		# parameter a0 (previously used 0.001 for other materials)
# @ecut		= (  5.0,   5.0,    5.0,   5.0,   5.0,    5.0,   3.5,     2.4 );		# ecut (about the magnitude of band gap for CsI, SiGe)





# band gap (eV)
      Eg  <- c(-999,  -999,  -999,  -999,  5.6,    6.1,  11.,   11.8,  5.9, -999,   1.57, 4.5,    4.56, 5.5 )
names(Eg) <- c("SiO2","Ge",  "Si" ,"GaAs","LaBr3","CsI","BaF2","CaF2","NaI","SiGe","CZT","SrI2", "YAG","YAP")

# material density in g/cm3
      d  <- c(2.32,   5.323, 2.33, 5.31,  5.29,   4.51, 4.89,  3.186, 3.667,3.827, 5.81, 4.549,  4.56, 5.5)
names(d) <- c("SiO2","Ge",  "Si" ,"GaAs","LaBr3","CsI","BaF2","CaF2","NaI","SiGe","CZT","SrI2", "YAG","YAP")

# number of valence electron: YAG, molecular formula: Y3Al5O12.  Y (Z=39, 3 valence electron, Al (Z=13, 3 valence electron, O: (Z=16, 6 valence electron: no.val.elec = 3*3 + 3 * 5 + 6 * 12 = 9+15+72 = 96
#                             YAP, molecular formula: YAlO3.                                                                                              no.val.elec = 3 + 3 + 18 = 24
      no.val  <- c( 16,    4,    4,    8,     24,     8,    16,    16,    8,    8,     30,   16,    96,   24) 
names(no.val) <- c("SiO2","Ge", "Si" ,"GaAs","LaBr3","CsI","BaF2","CaF2","NaI","SiGe","CZT","SrI2","YAG","YAP")

# molecular weight (g/mol): http://www.convertunits.com/molarmass
      M  <- c( 60.0843,72.64, 28.0855,144.6446, 378.6175, 259.80992,175.3238064,78.0748064, 149.894244, 100.7255,305.42, 341.42894, 593.61804,163.885588)
names(M) <- c("SiO2",  "Ge",  "Si",   "GaAs",   "LaBr3",  "CsI",    "BaF2",     "CaF2",     "NaI",      "SiGe",  "CZT",  "SrI2",    "YAG",    "YAP")

# page 40:
# ruo =  no.val * density * Avogadro Number / molecular Weight
  ruo <- no.val * d * C.Avogadro / M

# page 41: Rydberg Energy: Ry = (me^4) / (2h_hat^2) = 13.6 eV
R <- m.elec*e.elec^4 / (2*h_hat.eVs^2)
R <- 13.6				# hard coded as 13.6 eV

# page 40: (term) = (3*pi^2*ruo)^(2/3)
# term 3
term3 <- (3*pi^2*ruo)^(2/3)

# ------------------------------------
# Fermi Energy
# ------------------------------------
# page 40: Ef = ((h_hat^2 * C^2) / (2*m*c^2)) + (3*(pi^2)*ruo)^(2/3)
Ef <- ((h_hat.eVs^2) * ((100*v.light)^2) / (2*mc2.eV)) * term3

# chi-square
chi.square <- (13.6/(Ef * pi^2))^(1/2)
chi.square <- (R/(Ef * pi^2))^(1/2)

# ------------------------------------
# palsma energy
# ------------------------------------
Ep <- ((16/3)*chi.square)^(1/2) * Ef

# output
OUT <- cbind(d,M,no.val,ruo,term3, R,chi.square,Ef,Ep)

FL.OUT <- "Ef_Ep_Calculated.csv"
cat(paste("Ef and Ep,",sep=""),file=FL.OUT,append=TRUE)
write.table(OUT,file=FL.OUT,sep=",",row.names=TRUE,col.names=TRUE,append=TRUE)


# >  C <- 2.997925e10
# >  m <- 9.1093829e-31
# >  mc2 <- m*C*C
# >  mc2
# >  [1] 8.187107e-10
# >  m <- 9.1093829e-31
# >  
# >  C_m_s <- 2.997925e8
# >  m_kg <- 9.1093829e-31
# >  mc2_kgm2_s2 <- m_kg*C_m_s*C_m_s
# >  eV_2_J <- 1.602e-19
# >  mc2_ev <- mc2_kgm2_s2 / eV_2_J
# >  mc2_eV
# >  Error: object 'mc2_eV' not found
# >  mc2_ev
# >  [1] 511055.4
# >  
# >  
# >  plank_const_m2kg_s <- 6.62606957e-34
# >  plank_const_eVs <-plank_const_m2kg_s / eV_2_J
# >  plank_const_eVs
# >  [1] 4.136123e-15
# >  
# >  
# >  
# >  density_g_cm3 <- 2.33
# >  
# >  molWgt <- 29.086
# >  
# >  plank_const_eVs_hat <- plank_const_eVs / (2*pi)
# >  Ef <- plank_const_eVs^2 * (100*C_m_s^2) / (2*mc2_ev) * (3*pi^2*ruo)^(2/3)
# >  Error: object 'ruo' not found
# >  ruo <- 4* density_g_m3*Avogadro/molWgt
# >  Error: object 'density_g_m3' not found
# >  ruo <- 4* density_g_cm3*Avogadro/molWgt
# >  Ef <- plank_const_eVs^2 * (100*C_m_s^2) / (2*mc2_ev) * (3*pi^2*ruo)^(2/3)
# >  Ef
# >  [1] 4.807653
# >  ruo
# >  [1] 1.929669e+23
# >  plank_const_eVs
# >  [1] 4.136123e-15
# >  plank_const_eVs_hat
# >  [1] 6.582845e-16
# >  Ef <- plank_const_eVs_hat^2 * (100*C_m_s^2) / (2*mc2_ev) * (3*pi^2*ruo)^(2/3)
# >  Ef
# >  [1] 0.1217793
# >  
# >  
# >  
# >  C_m_s
# >  [1] 299792500
# >  C_cm_s <- C_m_s*100
#
# density in g
# >  Ef <- plank_const_eVs_hat^2 * (100*C_cm_s^2) / (2*mc2_ev) * (3*pi^2*ruo)^(2/3)
# >  Ef
# >  [1] 1217.793
# >  
# >  
# >  Ef <- plank_const_eVs_hat^2 * (100*C_cm_s^2) / (2*mc2_ev) * (3*(pi^2)*ruo)^(2/3)
# >  Ef
# >  [1] 1217.793
# >  Ef <- (plank_const_eVs_hat^2 * (C_cm_s^2) / (2*mc2_ev)) * (3*(pi^2)*ruo)^(2/3)
# >  Ef
# >  [1] 12.17793
# >  
# >  
# >  plank_const_eVs_hat
# >  [1] 6.582845e-16
# >  mc2_ev
# >  [1] 511055.4
# >  ruo
# >  [1] 1.929669e+23
# >  
# >  
# >  Ef <- (plank_const_eVs_hat^2 * (C_cm_s^2) / (2*mc2_ev)) * (3*(pi^2)*ruo)^(2/3)
# >  Ef
# >  [1] 12.17793
# >  
# >  
# >  
# >  
# >  
# >  
# >  
# >  Ef <- (plank_const_eVs_hat^2 * (C_cm_s^2) / (2*mc2_ev)) * (3*(pi^2)*ruo)^(2/3)
# >  Ef
# >  [1] 12.17793
# >  