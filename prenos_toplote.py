import numpy as np
import matplotlib.pyplot as plt
from funkcije import mesh, list_min, sun_position
import time

# Parametri podele asteroida
R = 6378e3 # poluprecnik asteroida (m)
Nr = 4 # podela u radijalnom pravcu
Nphi =12 # podela po longitudi
Ntheta = 6 # podela po latitudi

# Fizicki parametri
rho = 5500  # gustina asteoida (kg/m^3)
S0 = 1361  # Zracenje sa Sunca (W/m^2)
k = 30. # koeficijent konduktivnosti (W/(mK))
eps = 1. # emissivity of the surface element
sig = 5.6704e-8 # Stefan–Boltzmann constant (W/m^2K)
cp = 700. # Toplotni kapacitet pri konstantnom pritisku (J/kg K)
albedo = 0.3

# ostali parametri
dt_faktor = 100 # koeficijent sa kojim se deli minimalni vremenski korak za koji kondukcija probije celu celiju
period = 15*3600 # ukupno vreme simulacije

 

# Ucitavanje mape
celije, okolina_indeksi, zapremine, okolina_povrsine, okolina_rastojanja, \
okolina_centralne_indeksi, okolina_centralne_povrsine, \
okolina_centralne_rastojanja, spoljne_celije_indeksi, spoljne_celije_povrsine, \
spoljne_celije_normale=mesh(R, Nr, Nphi, Ntheta)


d_min = np.min(okolina_rastojanja)
dt_critical=d_min**2*rho*cp/k # kriticni vremenski korak kada kondukcija probije celu celiju (jednacina 2)

dt = dt_critical/dt_faktor #s


dt=3600.*6
period = 86400. * 10000000




N = len(celije) + 1 # dodaje se fiktivna celija



T = np.zeros(N)+260.; # temperature svih celija
J = np.zeros(N)
J1 = np.zeros(N)

r_sun = np.array([1,0,0])


vreme=np.arange(0, period+dt, dt) # vremenski trenuci u kojima se racunaju temperature

#TT = np.zeros([N, len(vreme)]); #niz u kojem se pamte temperature u svakoj celiji i u svakom vremenskom koraku



deltaT_proba = np.zeros([len(T), 6])


T_mean = np.zeros(int(len(vreme)/1000))

br=0

for i in range(len(vreme)): # po vremenu
    
    
    r_sun, aaa = sun_position(np.deg2rad(66.5), 0, 86164, 1, 0, vreme[i], M0 = 0, rotation_0 = 0)
    
    if np.mod(i,10000)==0:
        print(np.round(vreme[i]/86400), np.round(np.sum(T[-int(Nphi*Ntheta)-1:-1]*spoljne_celije_povrsine)/4/np.pi/R**2, 4))
        br+=1
        T_mean[br] = np.mean(T[:-1])
        
        
    
# =============================================================================
#     J[spoljne_celije_indeksi] = spoljne_celije_povrsine*(np.maximum( S0 * (1-albedo) * 
#      np.dot(spoljne_celije_normale, r_sun),0) 
#     - sig*eps*(T[spoljne_celije_indeksi])**4)
# 
#     pocetak = time.time()
#     for j in range(N-1): # po celijama (-1 je zato sto ne racunamo za fiktivnu celiju)     
#         # kondukcija
#         if j==0: # centralna celija
#             delta_T = T[okolina_centralne_indeksi] - T[0]
#             J[j] += np.sum(k*delta_T/np.array(okolina_centralne_rastojanja)*okolina_centralne_povrsine)
#         else:
#             delta_T = T[okolina_indeksi[j-1]]-T[j] # razlika u odnosu na okolne celije (j-1 je samo za okolinu jer ona krece od 1 celije, a posebno je za centralnu)
#             J[j] += np.sum(k*delta_T/np.array(okolina_rastojanja[j-1])*okolina_povrsine[j-1])
#             
#             deltaT_proba[j]=delta_T
# =============================================================================
            

    J1[spoljne_celije_indeksi] = spoljne_celije_povrsine*(np.maximum( S0 * (1-albedo) * 
     np.dot(spoljne_celije_normale, r_sun),0) 
    - sig*eps*(T[spoljne_celije_indeksi])**4)
    
    delta_T1 = T[okolina_centralne_indeksi] - T[0]
    J1[0] = J1[0] + np.sum(k*delta_T1/np.array(okolina_centralne_rastojanja)*okolina_centralne_povrsine)
    
    delta_T1 = T[okolina_indeksi] - T[1:-1, None] 
    
    J1[1:-1] = J1[1:-1] + np.sum(k*delta_T1/np.array(okolina_rastojanja)*okolina_povrsine, axis=1)


    dTdt = J1[:-1]/(rho*zapremine*cp)
    T[:-1] += dt*dTdt
    
    
    '''
    Racunanje Jarkovskog:
        prvo se sila koja se dobija iz Spitale (eq: 11) razlozi na komponente R (u pravcu radijus vekora) 
        i B (normalna na radijus vektor u ravni orbite) (Danby, 1992 eq 11.5.1 str. 323)
        
        Onda se da/dt dobije iz Danby eq:11.5.13 str. 327
    '''
    
    
#    # evolucija temperature za svaku celiju
#    for kk in range(N-1):
#        TT[kk][i]=T[kk]
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