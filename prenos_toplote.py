import numpy as np
import matplotlib.pyplot as plt
from funkcije import mapa, list_min

# Parametri podele asteroida
R = 1 # poluprecnik asteroida
Nr = 4 # podela u radijalnom pravcu
Nphi =12 # podela po longitudi
Ntheta = 6 # podela po latitudi

# Fizicki parametri
rho = 2500  # gustina asteoida (kg/m^3_
S0 = 1400  # Zracenje sa Sunca (W/m^2)
k = 1.5 # koeficijent konduktivnosti (W/(mK))
eps = 1. # emissivity of the surface element
sig = 5.6704e-8 # Stefan–Boltzmann constant (W/m^2K)
cp = 1000 # Toplotni kapacitet pri konstantnom pritisku (J/kg K)
albedo = 0

# ostali parametri
dt_faktor = 100 # koeficijent sa kojim se deli minimalni vremenski korak za koji kondukcija probije celu celiju
period = 15*3600 # ukupno vreme simulacije

 

# Ucitavanje mape
celije, okolina_indeksi, zapremine, okolina_povrsine, okolina_rastojanja, \
okolina_centralne_indeksi, okolina_centralne_povrsine, \
okolina_centralne_rastojanja, spoljne_celije_indeksi, spoljne_celije_povrsine, \
spoljne_celije_normale=mapa(R, Nr, Nphi, Ntheta)


d_min = np.min(okolina_rastojanja)
dt_critical=d_min**2*rho*cp/k # kriticni vremenski korak kada kondukcija probije celu celiju (jednacina 2)

dt = dt_critical/dt_faktor #s







dt=100
period = 3000000




N = len(celije) + 1 # dodaje se fiktivna celija






J = np.zeros(N) # ukupni fluks za svaku celiju (nije nas briga za svaku povrs posebno)

T = np.zeros(N)+0; # temperature svih celija
J = np.zeros(N)

r_sun = np.array([1,0,0])


vreme=np.arange(0, period+dt, dt) # vremenski trenuci u kojima se racunaju temperature

TT = np.zeros([N, len(vreme)]); #niz u kojem se pamte temperature u svakoj celiji i u svakom vremenskom koraku


for i in range(len(vreme)): # po vremenu
    
    if np.mod(i,1000)==0:
        print(vreme[i])
    
    # spolnjni fluks
    J[spoljne_celije_indeksi] = S0 * (1-albedo) * spoljne_celije_povrsine * np.dot(spoljne_celije_normale, r_sun) # dolazni fluks
    J[J < 0] = 0 # nulira se fluks sa nocne strane
    J[spoljne_celije_indeksi] -= spoljne_celije_povrsine*sig*eps*(T[spoljne_celije_indeksi])**4 # zracenje
    
    for j in range(N-1): # po celijama (-1 je zato sto ne racunamo za fiktivnu celiju)     
        # kondukcija
        if j==0: # centralna celija
            delta_T = T[okolina_centralne_indeksi] - T[0]
            J[j] += np.sum(k*delta_T/np.array(okolina_centralne_rastojanja)*okolina_centralne_povrsine)
        else:
            delta_T = T[okolina_indeksi[j-1]]-T[j] # razlika u odnosu na okolne celije (j-1 je samo za okolinu jer ona krece od 1 celije, a posebno je za centralnu)
            J[j] += np.sum(k*delta_T/np.array(okolina_rastojanja[j-1])*okolina_povrsine[j-1])
        
    
    dTdt = J[:-1]/(rho*zapremine*cp)
    T[:-1] += dt*dTdt
    
    # evolucija temperature za svaku celiju
    for kk in range(N-1):
        TT[kk][i]=T[kk]
#    
#
#plt.plot(vreme/3600, TT[0]-273.15, 'r', label='1. celija')
#plt.plot(vreme/3600, TT[1]-273.15, 'g', label='2. celija')
#plt.plot(vreme/3600, TT[2]-273.15, 'b', label='3. celija')
#plt.plot(vreme/3600, TT[3]-273.15, 'k', label='4. celija')
#plt.xlabel('vreme (h)',fontsize=18)
#plt.ylabel('T ($^{\circ}$C)', fontsize=18)
#plt.legend(fontsize=18)
#plt.grid()