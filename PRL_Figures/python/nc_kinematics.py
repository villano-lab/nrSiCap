#######################
#### Dependencies #####
#######################
import numpy as np
import scipy.integrate as si

#######################
###### Constants ######
#######################

mn = 0.9395654133e9#Mass of neutron in [eV/c^2]

#could make a dictionary at some point? When more are being used.
#dictionary could be defined by Z,A or Element Name,isotope
m0 = 26.0603421e9    #Mass of 28Si                 in [eV/c^2]
mf = 26.9914339e9    #Mass of 29Si                 in [eV/c^2]
tau_stop = 0.3e-12   #this is made up. Stop time   in [s]

#(Excited) State Energies for 29Si (dictionaries)
Eg = {
    "energy": 0,
    "thalf" : 10000 #placeholder half-life should be longer than anything else.
}
E1 = {
    "energy": 1.3e6, #state's energy               in [eV]
    "thalf" : 291e-15#Half-life of excited state   in [s]
}
E2 = {
    "energy": 2.0e6,
    "thalf" : 306e-15
}
#in dictionary, multiple deltas (d1 - dX) would be listed, nan for invalid indices.

Q = m0 + mn - mf   #When variable m0 and mf are used, Q will have to be a function.
#It will stay in constants because it is constant for the specified isotope.

def mass(delta=Eg): #Calculate mass for energy level, defaulting to ground state
    return mf + delta["energy"]

# _  ___                            _   _          
#| |/ (_)_ __   ___ _ __ ___   __ _| |_(_) ___ ___ 
#| ' /| | '_ \ / _ \ '_ ` _ \ / _` | __| |/ __/ __|
#| . \| | | | |  __/ | | | | | (_| | |_| | (__\__ \
#|_|\_\_|_| |_|\___|_| |_| |_|\__,_|\__|_|\___|___/
#Kinematics

###########################################
############## Single-Step ################
###########################################

def D1s():     #Total energy deposit [eV]
    p1step = mf*(-1 + (1 + (2/mf)*Q)**(1/2))
    return p1step**2/(2*mf)
D1s = D1s()    #Return a value, not a function

###########################################
############## All 2-Step #################
###########################################

def p1(delta): #momentum just after capture [eV/c]
    return mass(delta)*(-1 + (1 + (2/mass(delta))*(Q-delta["energy"]))**(1/2))
def T1(delta): #energy just after capture   [eV]
    return p1(delta)**2/(2*mass(delta))
def v1(delta): #velocity just after capture [c]
    return p1(delta)/mass(delta)

#This uses a constant acceleration until stopping. This should change eventually.
def a(delta):                    #acceleration, from stop time[c/s]
    return v1(delta)/tau_stop
def v1star(tau,delta):           #velocity just before decay  [c]
    conditions = [tau <= tau_stop         , tau >= tau_stop]
    functions  = [v1(delta) - a(delta)*tau, 0              ]
    return np.piecewise(tau,conditions,functions)

def T1star(tau,delta):                      #Energy just before decay      [eV]
    return v1star(tau,delta)**2*mass(delta)/2

def D1(tau,delta):
    return T1(delta) - T1star(tau,delta)   #First energy deposit (from initial slowdown)

#######################
#### 2s Lab Frame #####
#######################

#The lab frame method is unstable. Use the CM Frame.

def D2_Lab(tau,theta,delta):                                             #Second energy deposit (decay) [eV]
    p2 = lambda beta: mass(Eg)*v1star*np.sin(theta - beta)/np.sin(theta) #The momentum just after decay [eV/c]
    rightside = lambda beta: 1/(2*mass(Eg))*p2(beta)**2 + mf*v1star*np.sin(beta)/np.sin(theta)                 
    
    #Computational solution for system of equations
    array = np.linspace(0,theta)
    z = 0
    for x in array:                                  #For every point in the array we generated,
        if rightside(x) >= T1star + delta["energy"]: #If it's past Q-d
            print(str(z) + ", " + str(rightside(z))) #Report the value just before it
            print(str(x) + ", " + str(rightside(x))) #And of course the curent value
            y = x                                    #Let's save this value to call it later,
            break                                    #Then kill the loop
        z = x #Remember this value during the next step, assuming we didn't get to the point we wanted.
    newarray = np.linspace(y,z)
    import scipy.interpolate as inter
    interp = inter.interp1d(rightside(newarray),newarray)
    beta = interp(T1star + delta["energy"])
    return p2(beta)**2/(2*mass(Eg))

def D2s_Lab(tau,theta,delta): #Total energy deposit
    return D1(tau,delta) + D2_Lab(tau,theta,delta)

#######################
###### CM Frame #######
#######################

#While labelled as CM, these *values* are in the lab frame (the calculations use the lab frame but don't end there)
def D2_CM(tau,beta,delta):                            #Second energy deposit (decay) [eV]
    v2starmag = -1 + (1 + (2/mass(Eg))*delta["energy"])**(1/2) #CM Frame velocity after decay [c]
    v2mag     = ((v1star(tau,delta) + v2starmag*np.cos(beta))**2 + (v2starmag*np.sin(beta))**2)**(1/2) #|v2| = (v2_x)^2 + (v2_y)^2
    return mass(Eg)*v2mag**2/2
D2 = lambda tau,beta,delta: D2_CM(tau,beta,delta)  #Since lab frame is deprecated, make D2 = D2_CM

def D2s_CM(tau,beta,delta): #Total energy deposit [eV]
    return D1(tau,delta) + D2_CM(tau,beta,delta)
D2s = lambda tau,beta,delta: D2s_CM   #Since lab frame is deprecated, make D2s = D2s_CM

#  ____  _     _        _ _           _   _                 
# |  _ \(_)___| |_ _ __(_) |__  _   _| |_(_) ___  _ __  ___ 
# | | | | / __| __| '__| | '_ \| | | | __| |/ _ \| '_ \/ __|
# | |_| | \__ \ |_| |  | | |_) | |_| | |_| | (_) | | | \__ \
# |____/|_|___/\__|_|  |_|_.__/ \__,_|\__|_|\___/|_| |_|___/
#Distributions

#######################
### 2-Step Capture ####
#######################

def deltafunc(x,delta): #Delta function-approximate centered at T1
    b = T1(delta)
    c = 0.07
    exponent = lambda x: np.exp(-np.power(x-b,2)/(2*np.power(c,2)))
    a = 1/si.quad(exponent,0,2*T1(delta))[0]
    return a*exponent(x)
                                

def PT1(delta):                                       #Probability that D1 = T1
    return np.exp(-np.log(2)*tau_stop/delta["thalf"])
def f1(D1,delta):                                     #PDF for D1 only
    conditions = [ (D1 < T1(delta)) & (D1 >= 0),
                    D1 == T1(delta) ]
    #v--- Main probability function ---v
    reg = lambda D1,delta: (np.log(2)/delta["thalf"])*1/(a(delta)*(2*mass(delta)*(T1(delta)-D1))**(1/2))*np.exp(-np.log(2)*(v1(delta)-(2*(T1(delta)-D1)/mass(delta))**(1/2))/(a(delta)*delta["thalf"])
                              )
    equations = [ reg(D1,delta),
                  PT1(delta)*deltafunc(D1,delta) + reg(D1,delta),
                  0 ]
    return np.piecewise(D1,conditions,equations)


def D2min(D1,delta):     #Minimum D2 when defining boundaries in piecewise
    return (delta["energy"]**2/(2*mass(delta)))*(2*mass(delta)*(T1(delta)-D1)/delta["energy"]**2 
                    - 2*np.sqrt(2*mass(delta)*(T1(delta)-D1)/delta["energy"]**2)
                    + 1)
def D2max(D1,delta):     #Maximum D2 when defining boundaries in piecewise
    return (delta["energy"]**2/(2*mass(delta)))*(2*mass(delta)*(T1(delta)-D1)/delta["energy"]**2 
                    + 2*np.sqrt(2*mass(delta)*(T1(delta)-D1)/delta["energy"]**2)
                    + 1)
def D2width(D1,delta):   #Width of the function D2
    return (2*delta["energy"])/np.sqrt(mass(delta)/(2*(T1(delta)-D1)))
def f2(D1,D2,delta):     #"PDF" for D2 (not true PDF; depends on D1)
    if D1 > T1(delta):
        return 0 #Prevent errors that occur when trying to check criteria at invalid values
    conditions = [ (D2 <= D2max(D1,delta)) & (D2 >= D2min(D1,delta)) & (D1 <= T1(delta)) & (D1 >= 0) ]
    equations  = [ 1/(2*delta["energy"])*np.sqrt(mass(delta)/(2*(T1(delta) - D1))),
                   0 ] #0 elsewhere
    return np.piecewise(D2,conditions,equations)

def g(D1,D2,delta):  #Joint PDF for D1 & D2
    return f1(D1,delta)*f2(D1,D2,delta)

def ftot(Dtot,delta):
    gnew = lambda D1,Dtot: g(D1,Dtot-D1,delta)
    return si.quad(gnew,0,np.inf,args=Dtot)[0]

# ---- Silicon Nucleus Simulator 2020 ---- #
def D2s_sim(delta):                      #Generates an energy deposit!
    taugenerator = np.random.rand()      #uniformly generate a random number from  0 to 1
    l = delta["thalf"]/np.log(2)
    tau = -np.log(taugenerator)*l
    betagenerator = np.random.rand()*2-1 #uniformly generate a random number from -1 to 1
    beta = np.arccos(betagenerator)
    return D2s_CM(tau,beta,delta)

#  ___            _          _   _             
# |_ _|___  _ __ (_)______ _| |_(_) ___  _ __  
#  | |/ _ \| '_ \| |_  / _` | __| |/ _ \| '_ \ 
#  | | (_) | | | | |/ / (_| | |_| | (_) | | | |
# |___\___/|_| |_|_/___\__,_|\__|_|\___/|_| |_|
# Ionization

#######################
###### Lindhard #######
#######################

def Y(E): #ionization for a *single* energy (assumes going to zero).
    Z = 14                     #Atomic number of silicon
    A = 29                     #AFTER capture (29Si)
    eps = 11.5*E/1e3*Z**-(7/3) #This expects an energy in keV, so divide by 1e3 to convert eV to keV. (Function now expects eV)
    #values found on pg 3 of Sorenson 2015 "Atomic limits in the search for galactic dark matter"
    a = 3.0
    gam = 0.15
    b = 0.7
    omg = 0.6
    #rest of the form given by Anthony in OneNote: Ionization -> Yield (Lindhard)
    g = a*eps**gam + b*eps**omg + eps
    k = 0.133*Z**(2/3)*A**(-1/2)
    #print('k value: ' + str(k))
    Y = k*g/(1+k*g)
    return Y

#Take an arbitrary number of energies! These are in eV and chronological order.
def Ei_tot(*Energies):     #This takes the energy LEVELS, not the deposits!
    Ei = lambda E: E*Y(E)  #Function for ionization energy from a state to ground.
    #threshold = 1e-4      #Minimum value for valid reading
    #print(type(Energies))
    if (type(Energies) == list) or (type(Energies) == np.ndarray) or (type(Energies) == tuple):
        dEi = lambda Ea,Eb: Ei(Ea) - Ei(Eb)  #Function for ionization energy from one state to another.
        LEnergies = []
        for item in Energies:
            LEnergies.extend(item)
        if len(LEnergies) % 2 != 0:
            if LEnergies[-1] == 0:
                print("Warning: For odd number of entries, expected last entry is nonzero. Please check for missing or extrenuous entries.")
                #return np.nan
            LEnergies += [0]
        i = 0
        sumval = 0
        while i < len(LEnergies):
            currentval = LEnergies[i],LEnergies[i+1]
            sumval += dEi(LEnergies[i],LEnergies[i+1])
            i += 2
        #if sumval > threshold: #threshold only effects the total, not the individual steps
        return sumval
        #else:
            #note that it didn't contribute to the spectrum, or maybe just
            #return 0
    elif (type(Energies) == int) or (type(Energies) == float):
        print('elif')
        return Ei(Energies)