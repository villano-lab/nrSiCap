#Detector and physical constants for R68 analysis

#Detector constants
V = 125.0 #R68 bias voltage
eps = 3.8 #eV for silicon
G_NTL = (1.+(V/eps)) #NTL gain
F = 0.1161 #silicon value taken from https://www.sciencedirect.com/science/article/pii/S0168900297009650

#Resolution model
#sigma_ee^2 = sigma0^2 + B*Eee + (A*Eee)^2 (from CDMSlite paper)
#B contains Fano contribution and maybe other stuff
# B = F*epsilon + B_1

#Params from Matt's Bkg resolution fit:
#https://zzz.physics.umn.edu/cdms/doku.php?id=cdms:k100:run_summary:run_68:run_68_panda:calibration#resolution_versus_energy
sigma0=10.27 #eV 
B=0.627*3.8 #This includes ER FANO, which we've already applied above
A=0 #TODO: This part has not been fit yet!