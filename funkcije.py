import numpy as np

def list_min(lista):
    
    '''
    nalazi minimum iz liste koja u sebi ima druge liste
    '''
    
    d_min = 1e6

    for i in range(len(lista)):
        for j in range(len(lista[i])):
            if lista[i][j] < d_min:
                d_min = lista[i][j]
    return d_min

def log_podela(broj): # radi
    '''
    Deli opseg na broj delova po logaritamskoj skali. 
    Koristimo za podelu asteroid au radijalnom pravcu kako bi celije blize povrsini
    bile manje u odnosu na one u unutrasnjosti asteroida.
    input: broj celija
    output: podela na intervalu [0,1]
    
    Mnozimo output sa radijusom da bismo dobili podelu na intervalu [0,r]
    '''
    
    podela=[0]
    skala= 1 / np.log10(1.0 + broj);
    
    for i in range(broj):
        
        granica=np.log10(2.0 + i) * skala
        podela.append(granica)
        
    return np.array(podela)

def centri(a, nula=0): # radi
    '''
    racuna koordinate centara intervala unutar nekog opsega. Ovo nam treba za
    racunanuje rastojanje izmedju susednih celija
    
    input: opseg podeljen na intervale (koji dobijamo pozivanjem funkcije log_podela)
    
    output: koordinate centara
    '''
    
    
    x=np.zeros(len(a)-1)
    
    for i in range(1,len(a)):
        x[i-1]=a[i-1]+(a[i]-a[i-1])/2
        
    if nula==1: # kada se racunaju centri po r prva celija treba da se nulira
        x[0]=0
        
    return x


def metrika(r1,r2):
    return np.sqrt((r1[0]-r2[0])**2 + (r1[1]-r2[1])**2 + (r1[2]-r2[2])**2)

def mapa(R, Nr, Nphi, Ntheta):
    
    '''
    Daje mapu podele asteroida na celije, kao i karakteristike celija koje su potrebne
    za racun prenosa toplote
    
    input:
        R: radijus asteroida
        Nr: podela po radijusu (broj celija u radijalnom pravcu)
        Nphi: podela po longitudi (broj celija po longitudi na intervalu [0, 2pi))
        Ntheta: podela po latitudi (broj celija po latitudi na intervalu [-pi/2, +pi/2))
    
    output:
        celije: niz adresa svake celije u formatu [a,b,c]. a je adresa duz radijalnog pravca, b duz longitude, c duz latitude
        npr. celija [3,3,0] je cetvrta od centra (zato sto je prva sferna celija u centru), cetvrta po longitudi i prva po latitudi (na juznoj strani ose) 
        -----------------
        okolina_adrese: adrese okolnih celija za svaku celiju. 
        okolina_indeksi: indeksi okolnih celija za svaku celiju.
    '''

    # Formiranje nizova celija
    r = log_podela(Nr)*R
#    r=np.array([0,0,R])
    

    phi = np.linspace(0, 2*np.pi, Nphi+1)
    theta = np.linspace(-np.pi/2, np.pi/2, Ntheta+1)

    
    # Odredjivanje centara elemenata
    phi_centri = centri(phi)
    r_centri = centri(r, nula=1)
    theta_centri = centri(theta)
    
    # Razlika dva susedna ugla (phi, theta)
    dphi = phi[1] - phi[0]
    dtheta = theta[1] - theta[0]
    
    # ------ Formiranje liste elemenata ------ #
    
    celije = [] # (adrese svake celije) bez centralne celija
    zapremine = [] # zapremine celija
    koordinate_CM=[[0,0,0]]


    for i in range(1, Nr): # ide od 1 da bi se izbacila centralna celija [0,0,0]
        for j in range(Nphi):
            for k in range(Ntheta):

                celije.append([i, j, k]) 
                zapremine.append((r[i+1]**3-r[i]**3)/3*dphi*
                                 (np.sin(theta[k+1])-np.sin(theta[k])))
                
                '''
                koordinate centara masa svake celije. proveriti ovo!
                '''

                C=3/8*(r[i+1]**4-r[i]**4)/(r[i+1]**3-r[i]**3)/(np.sin(theta[k+1])-np.sin(theta[k]))/dphi
                x_CM=C*(np.sin(phi[j+1])-np.sin(phi[j]))*(dtheta+1/2*(np.sin(2*theta[k+1])-np.sin(2*theta[k])))
                y_CM=C*(np.cos(phi[j])-np.cos(phi[j+1]))*(dtheta+1/2*(np.sin(2*theta[k+1])-np.sin(2*theta[k])))
                z_CM=C*(np.sin(theta[k+1])**2-np.sin(theta[k])**2)*dphi
                
                koordinate_CM.append([x_CM, y_CM, z_CM])

    # ------ Formiranje liste susednih celija------- #
    okolina_adrese = []
    okolina_indeksi = []
    okolina_povrsine = []
    okolina_rastojanja=[]
    fiktivna_celija=[len(celije)+1, 0, 0] # dodaje se na kraj i sluzi samo da bi sve celije imale po 6 okolnih. One koje imaju manje, dodaje im se
    # ova fiktivna celija koja ima zapreminu nula, kao i sve povrsine. Na ovaj nacin ne utice na racun, ali omogucuje laksi rad sa nizovima jer ce biti iste duzine.
    rastojanja_kontrola=[]
    
    for i in range(0, len(celije)):
   
        if celije[i][0] == Nr-1: # ovo su spoljne celije
            
            if celije[i][2] == 0: # naslanjaju se na juznu stranu ose rotacije pa imaju 4 susedne celije
                
                # Adrese okolnih celija
                c0 = [celije[i][0]-1, celije[i][1], celije[i][2]] # prethodna u radijalnom pravcu (sledece nema jer je na povrsini)
                c1 = fiktivna_celija
                c2 = [celije[i][0], (celije[i][1]-1)%Nphi, celije[i][2]] # prethodna po longitudi
                c3 = [celije[i][0], (celije[i][1]+1)%Nphi, celije[i][2]] # sledeca po longitudi
                c4 = fiktivna_celija
                c5 = [celije[i][0], celije[i][1], celije[i][2]+1] # sledeca po latitudi (prethodne nema jer je na juznoj strani ose)
                
                # povrsi preko kojih se granici sa susednim celijama
                p0 = r[-2]**2*dphi*(np.sin(-np.pi/2 + dtheta) + 1) # prema centru
                p1 = 0
                p2 = dtheta*0.5*(r[-1]**2-r[-2]**2) # prethodni po longitudi
                p3 = p2 # naredni po longitudi
                p4 = 0
                p5 = np.cos(-np.pi/2 + dtheta)*dphi*0.5*(r[-1]**2-r[-2]**2) # naredni po latitudi
                      
                # rastojanja do susednih celija
                r0 = r_centri[-1] - r_centri[-2]
                r1 = np.inf # fiktivno
                r2 = 2*r_centri[-1]*np.sin(dphi/2)*np.sin(dtheta/2)
                r3 = r2
                r4 = np.inf # fiktivno
                r5 = 2*r_centri[-1]*np.sin(dtheta/2)
#                
    
            elif celije[i][2] == Ntheta-1: # naslanjaju se na severnu stranu ose rotacije pa imaju 4 susedne celije
                
                # Adrese okolnih celija
                c0 = [celije[i][0]-1, celije[i][1], celije[i][2]] # prethodna u radijalnom pravcu (sledece nema jer je na povrsini)
                c1 = fiktivna_celija
                c2 = [celije[i][0], (celije[i][1]-1)%Nphi, celije[i][2]] # prethodna po longitudi
                c3 = [celije[i][0], (celije[i][1]+1)%Nphi, celije[i][2]] # sledeca po longitudi
                c4 = [celije[i][0], celije[i][1], celije[i][2]-1] # prethodna po latitudi (sledece nema jer je na severnoj strani ose)
                c5 = fiktivna_celija
                
                # povrsi preko kojih se granici sa susednim celijama
                p0 = r[-2]**2*dphi*(1 - np.sin(np.pi/2 - dtheta)) # prema centru
                p1 = 0
                p2 = dtheta*0.5*(r[-1]**2-r[-2]**2) # prethodni po longitudi
                p3 = p2 # naredni po longitudi
                p4 = np.cos(np.pi/2 - dtheta)*dphi*0.5*(r[-1]**2-r[-2]**2) # prethodni po latitudi
                p5 = 0
                
                # rastojanja do susednih celija
                r0 = r_centri[-1] - r_centri[-2]
                r1 = np.inf
                r2 = 2*r_centri[-1]*np.sin(dphi/2)*np.sin(dtheta/2)
                r3 = r2
                r4 = 2*r_centri[-1]*np.sin(dtheta/2)
                r5 = np.inf

            else: # ne naslanjaju se na osu rotacije pa imaju 5 susednih celija (a spoljasnje su)
                   
                # Adrese okolnih celija
                c0 = [celije[i][0]-1, celije[i][1], celije[i][2]] # prethodna u radijalnom pravcu (sledece nema jer je na povrsini)
                c1 = fiktivna_celija
                c2 = [celije[i][0], (celije[i][1]-1)%Nphi, celije[i][2]] # prethodna po longitudi
                c3 = [celije[i][0], (celije[i][1]+1)%Nphi, celije[i][2]] # sledeca po longitudi
                c4 = [celije[i][0], celije[i][1], (celije[i][2]-1)] # prethodna po latitudi
                c5 = [celije[i][0], celije[i][1], (celije[i][2]+1)] # sledeca po latitudi
                
                # povrsi preko kojih se granici sa susednim celijama
                p0 = r[-2]**2*dphi*(np.sin(theta[celije[i][2] + 1]) - np.sin(theta[celije[i][2]])) # prema centru
                p1 = 0
                p2 = dtheta*0.5*(r[-1]**2-r[-2]**2) # prethodni po longitudi
                p3 = p2 # naredni po longitudi
                p4 = np.cos(theta[celije[i][2]])*dphi*0.5*(r[-1]**2-r[-2]**2) # prethodni po latitudi
                p5 = np.cos(theta[celije[i][2] + 1])*dphi*0.5*(r[-1]**2-r[-2]**2) # prethodni po latitudi
                
                # rastojanja do susednih celija
                r0 = r_centri[-1] - r_centri[-2]
                r1 = np.inf
                r2 = 2*r_centri[-1]*np.sin(dphi/2)*np.cos(theta_centri[celije[i][2]])
                r3 = r2
                r4 = 2*r_centri[-1]*np.sin(dtheta/2)
                r5 = r4
                    
        else: # ovo su unutrasnje celije

            if celije[i][2] == 0: # naslanjaju se na juznu stranu ose pa imaju 5 susednih celija

                # Adrese okolnih celija
                c0 = [celije[i][0]-1, celije[i][1], celije[i][2]] # prethodna u radijalnom pravcu
                c1 = [celije[i][0]+1, celije[i][1], celije[i][2]] # sledeca u radijalnom pravcu
                c2 = [celije[i][0], (celije[i][1]-1)%Nphi, celije[i][2]] # prethodna po longitudi
                c3 = [celije[i][0], (celije[i][1]+1)%Nphi, celije[i][2]] # sledeca po longitudi
                c4 = fiktivna_celija
                c5 = [celije[i][0], celije[i][1], celije[i][2]+1] # sledeca po latitudi (prethodne nema jer je na juznoj strani ose)
                
                # povrsi preko kojih se granici sa susednim celijama
                p0 = r[celije[i][0]]**2*dphi*(np.sin(-np.pi/2 + dtheta) + 1) # prema centru
                p1 = r[celije[i][0] + 1]**2*dphi*(np.sin(-np.pi/2 + dtheta) + 1) # od centra
                p2 = dtheta*0.5*(r[celije[i][0] + 1]**2-r[celije[i][0]]**2) # prethodni po longitudi
                p3 = p2 # naredni po longitudi
                p4 = 0
                p5 = np.cos(-np.pi/2 + dtheta)*dphi*0.5*(r[celije[i][0] + 1]**2-r[celije[i][0]]**2) # sledeci po latitudi
                
                # rastojanja do susednih celija
                r0 = r_centri[celije[i][0]] - r_centri[celije[i][0] - 1]
                r1 = r_centri[celije[i][0] + 1] - r_centri[celije[i][0]]
                r2 = 2*r_centri[celije[i][0]]*np.sin(dphi/2)*np.sin(dtheta/2)
                r3 = r2
                r4 = np.inf
                r5 = 2*r_centri[celije[i][0]]*np.sin(dtheta/2)

            elif celije[i][2] == Ntheta-1: # naslanjaju se na severun stranu ose pa imaju 5 susednih celija
               
                # Adrese okolnih celija
                c0 = [celije[i][0]-1, celije[i][1], celije[i][2]] # prethodna u radijalnom pravcu
                c1 = [celije[i][0]+1, celije[i][1], celije[i][2]] # sledeca u radijalnom pravcu
                c2 = [celije[i][0], (celije[i][1]-1)%Nphi, celije[i][2]] # prethodna po longitudi
                c3 = [celije[i][0], (celije[i][1]+1)%Nphi, celije[i][2]] # sledeca po longitudi
                c4 = [celije[i][0], celije[i][1], celije[i][2]-1] # prethodna po latitudi (sledece nema jer je na severnoj strani ose)
                c5 = fiktivna_celija
                
                # povrsi preko kojih se granici sa susednim celijama
                p0 = r[celije[i][0]]**2*dphi*(1 - np.sin(np.pi/2 - dtheta)) # prema centru
                p1 = r[celije[i][0] +1 ]**2*dphi*(1 - np.sin(np.pi/2 - dtheta)) # od centra
                p2 = dtheta*0.5*(r[celije[i][0] + 1]**2-r[celije[i][0]]**2) # prethodni po longitudi
                p3 = p2 # naredni po longitudi
                p4 = np.cos(np.pi/2 - dtheta)*dphi*0.5*(r[celije[i][0] + 1]**2-r[celije[i][0]]**2) # prethodni po latitudi
                p5 = 0
                
                # rastojanja do susednih celija
                r0 = r_centri[celije[i][0]] - r_centri[celije[i][0] - 1]
                r1 = r_centri[celije[i][0] + 1] - r_centri[celije[i][0]]
                r2 = 2*r_centri[celije[i][0]]*np.sin(dphi/2)*np.sin(dtheta/2)
                r3 = r2
                r4 = 2*r_centri[celije[i][0]]*np.sin(dtheta/2)
                r5 = np.inf
            
            else: # unutrasnje koje nisu na osi pa imaju 6 susednih celija
                
                # Adrese okolnih celija
                c0 = [celije[i][0]-1, celije[i][1], celije[i][2]] # prethodna u radijalnom pravcu
                c1 = [celije[i][0]+1, celije[i][1], celije[i][2]] # sledeca u radijalnom pravcu
                c2 = [celije[i][0], (celije[i][1]-1)%Nphi, celije[i][2]] # prethodna po longitudi
                c3 = [celije[i][0], (celije[i][1]+1)%Nphi, celije[i][2]] # sledeca po longitudi
                c4 = [celije[i][0], celije[i][1], celije[i][2]-1] # prethodna po latitudi 
                c5 = [celije[i][0], celije[i][1], celije[i][2]+1] # sledeca po latitudi
                
                # povrsi preko kojih se granici sa susednim celijama
                p0 = r[celije[i][0]]**2*dphi*(np.sin(theta[celije[i][2] + 1]) - np.sin(theta[celije[i][2]])) # prema centru
                p1 = r[celije[i][0] + 1]**2*dphi*(np.sin(theta[celije[i][2] + 1]) - np.sin(theta[celije[i][2]])) # od centra
                p2 = dtheta*0.5*(r[celije[i][0] + 1]**2-r[celije[i][0]]**2) # prethodni po longitudi
                p3 = p2 # naredni po longitudi
                p4 = np.cos(theta[celije[i][2]])*dphi*0.5*(r[celije[i][0] + 1]**2-r[celije[i][0]]**2) # prethodni po latitudi
                p5 = np.cos(theta[celije[i][2] + 1])*dphi*0.5*(r[celije[i][0] + 1]**2-r[celije[i][0]]**2) # sledeci po latitudi
                
                # rastojanja do susednih celija
                r0 = r_centri[celije[i][0]] - r_centri[celije[i][0] - 1]
                r1 = r_centri[celije[i][0] + 1] - r_centri[celije[i][0]]
                r2 = 2*r_centri[celije[i][0]]*np.sin(dphi/2)*np.cos(theta_centri[celije[i][2]])
                r3 = r2
                r4 = 2*r_centri[celije[i][0]]*np.sin(dtheta/2)
                r5 = r4

        okolina_adrese.append([c0, c1, c2, c3, c4, c5])
        okolina_povrsine.append([p0, p1, p2, p3, p4, p5])
        rastojanja_kontrola.append([r0, r1, r2, r3, r4, r5])


    # nuliranje svih koordinata centralne celije. Prethodni algoritam ce dati da 
    #je prethodna celija od [1, 2, 3] da je prethodna po radijalnom pravcu 
    # [0, 2, 3] a to je zapravo centralna celija [0,0,0]
    for i in range(len(celije)):
        if okolina_adrese[i][0][0] == 0:
            okolina_adrese[i][0] = [0, 0, 0]
            
    # formiranje okoline centralne celije
    okolina_centralne_adrese = []
    
    for i in range(Nphi):
        for j in range(Ntheta):
            okolina_centralne_adrese.append([1, i, j])
            
    # dodavanje centralne na pocetak liste celija
    celije.insert(0, [0,0,0])
    
    # formiranje niza u kojem su za svaku celiju dati indeksi okolnih celija, a ne njihove adrese po r, fi, theta
    okolina_indeksi = []
    for i in range(len(okolina_adrese)):
        trenutna_okolina = []
        for j in range(len(okolina_adrese[i])):
            try:
                trenutna_okolina.append(celije.index(okolina_adrese[i][j]))
            except:
                trenutna_okolina.append(len(celije))
        okolina_indeksi.append(trenutna_okolina)
        
        
    okolina_centralne_indeksi=[]
    for i in range(len(okolina_centralne_adrese)):
        okolina_centralne_indeksi.append(celije.index(okolina_centralne_adrese[i]))

    '''
    rastojanja do susednih celija (bez centralne, ona ide posebno)
    '''
    
    okolina_rastojanja=np.zeros([len(celije)-1,6])
    
    for i in range(0,len(okolina_rastojanja)): 
        for j in range(6):            
            if okolina_indeksi[i][j] == len(celije): # fiktivna celija
                okolina_rastojanja[i][j]=np.inf   
            else: 
                okolina_rastojanja[i][j]=metrika(koordinate_CM[i+1], koordinate_CM[okolina_indeksi[i][j]])
    
    '''
    CENTRALNA CELIJA
    '''
    

    # rastojanja do susednih od centralne

    okolina_centralne_rastojanja=np.zeros(len(okolina_centralne_indeksi))
    for i in range(len( okolina_centralne_rastojanja)):
         okolina_centralne_rastojanja[i]=metrika([0,0,0], koordinate_CM[okolina_centralne_indeksi[i]] )
    
    
    

    # povrsine za susede od centralne

    okolina_centralne_povrsine=np.zeros(len(okolina_centralne_indeksi))
    for i in range(len(okolina_centralne_indeksi)):
        okolina_centralne_povrsine[i]=r[1]**2*dphi*(np.sin(theta[okolina_centralne_adrese[i][2] + 1]) - np.sin(theta[okolina_centralne_adrese[i][2]]))
            
    # zapremina centralne celije       
    zapremina_centralne=4/3*r[1]**3*np.pi

    '''
    SPOLJNE CELIJE
    
    OVO NE VALJA. Treba normale da se smeste u tezista povrsina, a ne u sredine po koordinatama fi i teta. Priblizno resenje moze biti da se nadje
    teziste za ravanski trougao ili trapez i da onda to bude aproksimacija tezista za sferni poligon.
    '''  
    redni_broj=np.arange((Nr-1)*Nphi*Ntheta + 1) # redni broj svake celije
    spoljne_celije_indeksi=redni_broj[np.transpose(np.array(celije))[0] == Nr-1] # redni broj spoljnih celija
    spoljne_celije_normale=[] # samo za spoljne celije
    spoljne_celije_povrsine=[]

    for i in range(len(spoljne_celije_indeksi)):
        spoljne_celije_normale.append([np.cos(phi_centri[celije[spoljne_celije_indeksi[i]][1]])*np.cos(theta_centri[celije[spoljne_celije_indeksi[i]][2]]),
                          np.sin(phi_centri[celije[spoljne_celije_indeksi[i]][1]])*np.cos(theta_centri[celije[spoljne_celije_indeksi[i]][2]]),
                          np.sin(theta_centri[celije[spoljne_celije_indeksi[i]][2]])])

        spoljne_celije_povrsine.append(R**2*dphi*(np.sin(theta[celije[spoljne_celije_indeksi[i]][2] + 1]) - 
        np.sin(theta[celije[spoljne_celije_indeksi[i]][2]])))


    return (celije, okolina_indeksi, np.array(zapremine), np.array(okolina_povrsine), np.array(okolina_rastojanja), 
            okolina_centralne_indeksi, zapremina_centralne, okolina_centralne_povrsine, okolina_centralne_rastojanja,
            spoljne_celije_indeksi, np.array(spoljne_celije_povrsine), np.array(spoljne_celije_normale))

                


                
                