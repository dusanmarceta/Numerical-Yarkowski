import numpy as np
import matplotlib.pyplot as plt

from funkcije import mapa, centri, log_podela


R=1
Nr=3
Nphi=4
Ntheta=2


r = log_podela(Nr)*R

phi = np.linspace(0, 2*np.pi, Nphi+1)
theta = np.linspace(-np.pi/2, np.pi/2, Ntheta+1)

dphi = phi[1] - phi[0]
dtheta = theta[1] - theta[0]

# Odredjivanje centara elemenata
phi_centri = centri(phi)
theta_centri = centri(theta)
r_centri = centri(r)
    

celije, okolina_adrese, okolina_indeksi, spoljne_celije_indeksi, okolina_povrsine, okolina_rastojanja, zapremine, povrsine_spoljne, normale =mapa(R, Nr, Nphi, Ntheta)



# provera unutrasnjih povrsi od svake celije (prema centru)
p1=0
sfera=1 # redni broj sfere oko centralne

for i in range(len(celije)):
    
    if celije[i][0]==sfera: # dati nivo mereno od centra

        p1+=okolina_povrsine[i][0]
    
provera1=4*r[sfera]**2*np.pi


# provera spoljasnjih povrsi od svake celije (od centra)
p2=0
sfera=2 # redni broj sfere oko centralne

for i in range(len(celije)):
    
    if celije[i][0]==sfera: # dati nivo mereno od centra

        p2+=okolina_povrsine[i][1]
    
provera2=4*r[sfera+1]**2*np.pi


# provera bocnih (po longitudi)
p3=0
lon=32 # redni broj longitude (mada treba da je po svim isto)
br=0
for i in range(1,len(celije)-1):
    
    if celije[i][1]==lon: # dati nivo mereno od centra
        br+=1

        
        p3+=okolina_povrsine[i][3]
    
provera3=(R**2-r[1]**2)*np.pi/2



# provera konusnih povrsi (p4)
p4=0
lat=0 # redni broj latitude
br=0
for i in range(1,len(celije)-1):
    
    if celije[i][2]==lat: # dati nivo mereno od centra
        br+=1

        
        p4+=okolina_povrsine[i][4]
        
        

    
provera4=(R**2-r[1]**2)*np.cos(theta[lat])*np.pi


# provera konusnih povrsi (p4)
p5=0
lat=1 # redni broj latitude
br=0
for i in range(1,len(celije)-1):
    
    if celije[i][2]==lat: # dati nivo mereno od centra
        br+=1

        
        p5+=okolina_povrsine[i][5]
        
        

    
provera5=(R**2-r[1]**2)*np.cos(theta[lat+1])*np.pi


# provera rastojanja

# po latitudi

r1=0
long=1 # redni broj latitude
br=0
sfera==1
for i in range(1,len(celije)-1):
    if celije[i][0]==sfera:
        if celije[i][1]==long and celije[i][2]!=0: #
            br+=1
    
            r1+=okolina_rastojanja[i][4]
        
provera5=r_centri[sfera]*np.pi
print((provera5-r1)/provera5*100)


r1=0
long=3 # redni broj latitude
br=0
sfera==1
for i in range(1,len(celije)-1):
    if celije[i][0]==sfera:
        if celije[i][1]==long and celije[i][2]!=Ntheta-1: #
            br+=1
    
            r1+=okolina_rastojanja[i][5]
        
provera5=r_centri[sfera]*np.pi
print((provera5-r1)/provera5*100)



# u radijalnom pravcu

r1=0
long=3 # redni broj latitude
lat=4
br=0
sfera==1
rr=[]
for i in range(1,len(celije)-1):
    if celije[i][1]==long and celije[i][2]==lat:

        br+=1
    
        rr.append(okolina_rastojanja[i][1])
        
provera5=r_centri[sfera]*np.pi
print((provera5-r1)/provera5*100)
