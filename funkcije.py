import numpy as np

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

def centri(a): # radi
    '''
    racuna koordinate centara intervala unutar nekog opsega. Ovo nam treba za
    racunanuje rastojanje izmedju susednih celija
    
    input: opseg podeljen na intervale (koji dobijamo pozivanjem funkcije log_podela)
    
    output: koordinate centara
    '''
    
    
    x=np.zeros(len(a)-1)
    
    for i in range(1,len(a)):
        x[i-1]=a[i-1]+(a[i]-a[i-1])/2
        
    return x

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
    phi = np.linspace(0, 2*np.pi, Nphi+1)
    theta = np.linspace(-np.pi/2, np.pi/2, Ntheta+1)
    
    # Odredjivanje centara elemenata
    phi_centri = centri(phi)
    theta_centri = centri(theta)
    r_centri = centri(r)
    
    # Razlika dva susedna ugla (phi, theta)
    dphi = phi[1] - phi[0];
    dtheta = theta[1] - theta[0];
    
    # ------ Formiranje liste elemenata ------ #
    
    celije = [[0,0,0]] # (adrese svake celije) pocinje se od centralne celije
    redni_broj=np.arange((Nr-1)*Nphi*Ntheta+1) # redni broj svake celije
    zapremine = [] # zapremine celija

    for i in range(1, Nr): # ide od 1 da bi se izbacila unutrasnja celija [0,0,0]
        for j in range(Nphi):
            for k in range(Ntheta):

                celije.append([i, j, k]) 
                zapremine.append((r[i+1]**3-r[i]**3)/3*dphi*
                                 (np.sin(theta[k+1])-np.sin(theta[k])))

    
    celije_num = np.array(celije)            
    spoljne_celije_num = celije_num[np.transpose(celije_num)[0] == Nr-1] # sve one kod kojih je prva koordinata ima najveci moguci broj, tj. Nr-1. To znaci da su najdalje po radijalnom pravcu   
    spoljne_celije_indeks=redni_broj[np.transpose(celije_num)[0] == Nr-1] # redni broj spoljnih celija
    
    

    # Povrsine spoljnih celija
    S = np.zeros(len(spoljne_celije_num))
    for i in range(len(spoljne_celije_num)):
        S[i] = R**2*dphi*(np.sin(spoljne_celije_num[i][2] + 1) - 
        np.sin(spoljne_celije_num[i][2]))
        

    # ------ Formiranje liste susednih celija------- #
    okolina_adrese = []
    okolina_indeksi=[]
    
    for i in range(len(celije)):

        
        if celije[i][0] == Nr-1: # ovo su spoljne celije
            
            if celije[i][2] == 0: # naslanjaju se na juznu stranu ose rotacije pa imaju 4 susedne celije
                okolina_adrese.append([[celije[i][0]-1, celije[i][1], celije[i][2]], # prethodna u radijalnom pravcu (sledece nema jer je na povrsini)
                            [celije[i][0], (celije[i][1]-1)%Nphi, celije[i][2]], # prethodna po longitudi
                            [celije[i][0], (celije[i][1]+1)%Nphi, celije[i][2]], # sledeca po longitudi
                           [celije[i][0], celije[i][1], celije[i][2]+1]]) # sledeca po latitudi (prethodne nema jer je na juznoj strani ose)
    
    
            elif celije[i][2] == Ntheta-1: # naslanjaju se na severnu stranu ose rotacije pa imaju 4 susedne celije
                okolina_adrese.append([[celije[i][0]-1, celije[i][1], celije[i][2]], # prethodna u radijalnom pravcu (sledece nema jer je na povrsini)
                            [celije[i][0], (celije[i][1]-1)%Nphi, celije[i][2]], # prethodna po longitudi
                            [celije[i][0], (celije[i][1]+1)%Nphi, celije[i][2]], # sledeca po longitudi
                           [celije[i][0], celije[i][1], celije[i][2]-1]]) # prethodna po latitudi (sledece nema jer je na severnoj strani ose)

            else: # ne naslanjaju se na osu rotacije pa imaju 5 susednih celija
                okolina_adrese.append([[celije[i][0]-1, celije[i][1], celije[i][2]],  # prethodna u radijalnom pravcu (sledece nema jer je na povrsini)
                            [celije[i][0], (celije[i][1]-1)%Nphi, celije[i][2]], # prethodna po longitudi
                            [celije[i][0], (celije[i][1]+1)%Nphi, celije[i][2]], # sledeca po longitudi
                            [celije[i][0], celije[i][1], (celije[i][2]-1)], # prethodna po latitudi
                            [celije[i][0], celije[i][1], (celije[i][2]+1)]]) # sledeca po latitudi

    
        else: # ovo su unutrasnje celije

            if celije[i][2] == 0: # naslanjaju se na juznu stranu ose pa imaju 5 susednih celija
                okolina_adrese.append([[celije[i][0]-1, celije[i][1], celije[i][2]], # prethodna u radijalnom pravcu
                                    [celije[i][0]+1, celije[i][1], celije[i][2]], # sledeca u radijalnom pravcu
                                    [celije[i][0], (celije[i][1]-1)%Nphi, celije[i][2]], # prethodna po longitudi
                                    [celije[i][0], (celije[i][1]+1)%Nphi, celije[i][2]], # sledeca po longitudi
                                    [celije[i][0], celije[i][1], celije[i][2]+1]]) # sledeca po latitudi (prethodne nema jer je na juznoj strani ose)
    
            elif celije[i][2] == Ntheta-1:
                okolina_adrese.append([[celije[i][0]-1, celije[i][1], celije[i][2]], # prethodna u radijalnom pravcu
                                    [celije[i][0]+1, celije[i][1], celije[i][2]], # sledeca u radijalnom pravcu
                                    [celije[i][0], (celije[i][1]-1)%Nphi, celije[i][2]], # prethodna po longitudi
                                    [celije[i][0], (celije[i][1]+1)%Nphi, celije[i][2]], # sledeca po longitudi
                                    [celije[i][0], celije[i][1], celije[i][2]-1]]) # prethodna po latitudi (sledece nema jer je na severnoj strani ose)
            
            else: # unutrasnje koje nisu na osi pa imaju 6 susednih celija
                okolina_adrese.append([[celije[i][0]-1, celije[i][1], celije[i][2]], 
                                    [celije[i][0]+1, celije[i][1], celije[i][2]], 
                                    [celije[i][0], (celije[i][1]-1)%Nphi, celije[i][2]], 
                                    [celije[i][0], (celije[i][1]+1)%Nphi, celije[i][2]],
                                    [celije[i][0], celije[i][1], celije[i][2]-1], 
                                    [celije[i][0], celije[i][1], celije[i][2]+1]])   
                
    
    # nuliranje svih koordinata centralne celije. Prethodni algoritam ce dati da 
    #je prethodna celija od [1, 2, 3] da je prethodna po radijalnom pravcu 
    # [0, 2, 3] a to je zapravo centralna celija [0,0,0]
    for i in range(len(celije)):
        if okolina_adrese[i][0][0] == 0:
            okolina_adrese[i][0] = [0, 0, 0]
            
    # formiranje okoline centralne celije
    okolina_centralne=[]
    
    for i in range(Nphi):
        for j in range(Ntheta):
            okolina_centralne.append([1, i, j])
            
    okolina_adrese[0]=okolina_centralne
    
    # formiranje niza u kojem su za svaku celiju dati indeksi okolnih celija, a ne njihove adrese po r, fi, theta
    okolina_indeksi=[]
    for i in range(len(okolina_adrese)):
        trenutna_okolina=[]
        for j in range(len(okolina_adrese[i])):
            trenutna_okolina.append(celije.index(okolina_adrese[i][j]))
        okolina_indeksi.append(trenutna_okolina)
    
    """
    mapa je sada dobra, 
    treba proveriti zapremine i povrsine, ostalo je OK
    """
    
    return celije, okolina_adrese, spoljne_celije_num, spoljne_celije_indeks, okolina_indeksi

#    
#    
#    # ------- Odredjivanje jedinicnih vektora
#    # normala povrsinskih elemenata -------- #
#    
#    normale = [[0, 0, 0]] # Samo prva celija
#    
#    for i in range(1, len(celije)):
#        if celije[i][0] == Nr-1:
#            normale.append([np.cos(phi_centri[celije[i][1]])*np.cos(theta_centri[celije[i][2]]),
#                           np.sin(phi_centri[celije[i][1]])*np.cos(theta_centri[celije[i][2]]),
#                           np.sin(theta_centri[celije[i][2]])])
#        else:
#            normale.append([0, 0, 0]) 
#            
#    # ------- Rastojanja celija ------- #
#    # ------- Povrsi susednih celija ------- #
#    
#    rastojanja = [[]]
#    povrsi = [[]]
#    for i in range(1, len(celije)):
#        pomocni_rastojanja = []
#        pomocni_povrsi = []
#        for j in range(len(okolina_adrese[i])):  
#            # Unutrasnje i spoljasnje po r
#            if celije[i][1] == okolina_adrese[i][j][1] and celije[i][2] == okolina_adrese[i][j][2]: 
#                pomocni_rastojanja.append(abs(r_centri[celije[i][0]]
#                - r_centri[okolina_adrese[i][j][0]]))
#                
#                # Slucaj ukoliko je celija spoljna u odnosu na okolinu
#                if celije[i][0] > okolina_adrese[i][j][0]:
#                    pomocni_povrsi.append(r[celije[i][0]]**2*dphi*(np.sin(theta[celije[i][2] + 1]) 
#                    - np.sin(theta[celije[i][2]])))
#                    
#                # Slucaj ako je celija unutrasnja u odnosu na okolnu
#                else: 
#                    pomocni_povrsi.append(r[celije[i][0] + 1]**2*dphi*(np.sin(theta[celije[i][2] + 1])
#                    - np.sin(theta[celije[i][2]])))
#            # Levo i desno po phi
#            elif celije[i][0] == okolina_adrese[i][j][0] and celije[i][2] == okolina_adrese[i][j][2]:
#                pomocni_rastojanja.append(2*r_centri[celije[i][0]]*np.sin(dphi/2)
#                *np.cos(theta_centri[celije[i][2]]))
#                pomocni_povrsi.append((r[celije[i][0]+1]**2 - r[celije[i][0]]**2)/2*dtheta)
#                
#            # Gore i dole po theta
#            elif celije[i][0] == okolina_adrese[i][j][0] and celije[i][1] == okolina_adrese[i][j][1]:
#                pomocni_rastojanja.append(2*r_centri[celije[i][0]]*np.sin(dtheta/2))
#                
#                # Slucaj ako je celija iznad u odnosu na okolnu
#                if celije[i][2] > okolina_adrese[i][j][2]: 
#                    pomocni_povrsi.append((r[celije[i][0]+1]**2 - 
#                    r[celije[i][0]]**2)/2*dphi*np.cos(theta[celije[i][2]]))
#                    
#                # Slucaj ako je celija ispod u odnosu na okolnu
#                else: 
#                    pomocni_povrsi.append((r[celije[i][0]+1]**2 -
#                    r[celije[i][0]]**2)/2*dphi*np.cos(theta[celije[i][2] + 1]))
#        rastojanja.append(pomocni_rastojanja)
#        povrsi.append(pomocni_povrsi)
#    
#    # Formiranje liste rastojanja susednih celija
#    # Formiranje liste povrsi susednih celija
#    
#    rastojanja[0] = list([r_centri[1]])*18
#    povrsi_centar = []
#    for i in range(Nphi):
#        for j in range(Ntheta):
#            povrsi_centar.append(r[1]**2*dphi*(np.sin(theta[j+1]) - np.sin(theta[j])))
#    
#    povrsi[0] = povrsi_centar
#    
#    return (celije, okolina_adrese, okolina_indeksi, povrsi, zapremine, rastojanja, normale, spoljne_celije, S)


# sta je ovo???
def list_min(lista):
    
    d_min = 1e6

    for i in range(len(lista)):
        for j in range(len(lista[i])):
            if lista[i][j] < d_min:
                d_min = lista[i][j]
    return d_min
                


                
                