import numpy as np
import matplotlib.pyplot as plt


N = 4; #Broj celija
a = 0.01; #m^2
V = a**3; #m^3
Period = 15*3600; #s

rho = 2810;  #kg/m^3
S0 = 1400; #W/m^2
k = 237; #W/(mK)
eps = (0.03+0.2)/2+0.03;
sig = 5.67*10**(-8); #W/m^2K
cp = 889.78; #J/kg K




dt_critical=a**2*rho*cp/k # kriticni vremenski korak kada kondukcija probije celu celiju (jednacina 2)

dt = dt_critical/5; #s


J = np.zeros(N) # ukupni fluks za svaku celiju (nije nas briga za svaku povrs posebno)
A = np.ones([N, 6])*a**2; # povrsine svih strana svake celije posebno
T = np.zeros(N); # temperature svih celija

# =============================================================================
# Sledece mape su napravljne po sistemu:
#     prednja-0, zadnja-1, leva-2, desna-3, donja-4, gornja-5
# =============================================================================
N_s=[[0,2,3,4,5],[2,3,4,5],[2,3,4,5], [1,2,3,4,5]] # slobodne povrsi svake celije posebno
okolina=[[1],[0,2], [1,3], [2]] # celije sa kojima je data celija u kontaktu (redosled odgovara redosledu povrsi u nizu N_k)
N_k=[[1],[0,1], [0,1], [0]] # povrsi u kontaktu sa drugim celijama (za svaku celiju posebno)
d=[[a],[a,a],[a,a],[a]] # rastojanja do centara susednih celija (potrebno za nablaT (jednacina 6)). Isti redosled kao u N_k

J_external=np.array([S0*A[0][0], 0, 0, 0]) # eksterni fluks na svaku celiju posebno

vreme=np.arange(0, Period+dt, dt)
TT = np.zeros([N, len(vreme)]);

for i in range(len(vreme)): # po vremenu
    
    for j in range(len(J)): # po celijama
        
        J[j]=J_external[j] # fluks spolja
        
        J[j] -= np.sum(A[j][N_s[j]]*sig*eps*(T[j])**4) # zracenje
        
        # kondukcija
        delta_T = T[okolina[j]]-T[j] # razlika u odnosu na okolne celije
        J[j] += np.sum(k*delta_T/np.array(d[j])*A[j][N_k[j]])
        

    dTdt = J/(rho*V*cp)
    T += dt*dTdt
    
    # evolucija temperature za svaku celiju
    for k in range(N):
        TT[k][i]=T[k]
    

plt.plot(vreme/3600, TT[0]-273.15, 'r', label='1. celija')
plt.plot(vreme/3600, TT[1]-273.15, 'g', label='2. celija')
plt.plot(vreme/3600, TT[2]-273.15, 'b', label='3. celija')
plt.plot(vreme/3600, TT[3]-273.15, 'k', label='4. celija')
plt.xlabel('vreme (h)',fontsize=18)
plt.ylabel('T ($^{\circ}$C)', fontsize=18)
plt.legend(fontsize=18)
plt.grid()