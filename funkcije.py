import numpy as np

def log_podela(broj):
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

def centri(a):
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
        npr. celija [3,3,0] je cetvrta od centra (zato sto je prva sferna celija u centru), treca po longitudi i prva po latitudi (na juznom polu) 
        -----------------
        okolina: adrese okolnih celija za svaku celiju. 
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

    
    celije = np.array(celije)            
    spoljne_celije = celije[np.transpose(celije)[0] == Nr-1] # sve one kod kojih je prva koordinata ima najveci moguci broj, tj. Nr-1. To znaci da su najdalje po radijalnom pravcu   
    spoljne_celije_indeks=redni_broj[np.transpose(celije)[0] == Nr-1] # redni broj spoljnih celija
    
    

    # Povrsine spoljnih celija
    S = np.zeros(len(spoljne_celije))
    for i in range(len(spoljne_celije)):
        S[i] = R**2*dphi*(np.sin(spoljne_celije[i][2] + 1) - 
        np.sin(spoljne_celije[i][2]))
        
    '''
    treba proveriti zapremine i povrsine, ostalo je OK
    '''
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # ------ Formiranje liste susednih elemenata ------- #
    okolina = []
    
    for i in range(len(celije)):
        # Izbacivanje spoljne/naredne celije u radijalnom pravcu, ukoliko je 
        # trenutna celija zadnja u nizu/spoljna
        if celije[i][0] == Nr-1:
            # Celije u theta pravcu koje su u dodiru sa osom rotacije i plus su 
            # spoljne celije nemaju 6 susednih celija vec 4
            if celije[i][2] == 0 or celije[i][2] == Ntheta-1:
                okolina.append([[celije[i][0]-1, celije[i][1], celije[i][2]],  
                            [celije[i][0], (celije[i][1]-1)%Nphi, celije[i][2]], 
                            [celije[i][0], (celije[i][1]+1)%Nphi, celije[i][2]],
                           [celije[i][0], celije[i][1], (celije[i][2]+1)%Ntheta]])
            # Celije u theta pravcu koje su spoljne ali nisu u kontaktu sa osom
            # rotacije sfere nemaju 6 susednih celija vec 5
            else:
                okolina.append([[celije[i][0]-1, celije[i][1], celije[i][2]],  
                            [celije[i][0], (celije[i][1]-1)%Nphi, celije[i][2]], 
                            [celije[i][0], (celije[i][1]+1)%Nphi, celije[i][2]],
                            [celije[i][0], celije[i][1], (celije[i][2]-1)%Ntheta], 
                            [celije[i][0], celije[i][1], (celije[i][2]+1)%Ntheta]])
        # Slucaj ukoliko celija nije poslednja u nizu u radijalnom pravcu/poseduje
        # narednu celiju/nije spoljna
        else:
            # Celije u theta pravcu koje su u dodiru sa osom rotacije i nisu 
            # spoljne celije nemaju 6 susednih celija vec 5
            if celije[i][2] == 0 or celije[i][2] == Ntheta-1:
                okolina.append([[celije[i][0]-1, celije[i][1], celije[i][2]], 
                                    [celije[i][0]+1, celije[i][1], celije[i][2]], 
                                    [celije[i][0], (celije[i][1]-1)%Nphi, celije[i][2]], 
                                    [celije[i][0], (celije[i][1]+1)%Nphi, celije[i][2]],
                                    [celije[i][0], celije[i][1], (celije[i][2]-1)%Ntheta]])
            # Celije u theta pravcu koje nisu spoljne i nisu u kontaktu sa osom
            # rotacije sfere imaju svih 6 susdenih elemenata
            else:
                okolina.append([[celije[i][0]-1, celije[i][1], celije[i][2]], 
                                    [celije[i][0]+1, celije[i][1], celije[i][2]], 
                                    [celije[i][0], (celije[i][1]-1)%Nphi, celije[i][2]], 
                                    [celije[i][0], (celije[i][1]+1)%Nphi, celije[i][2]],
                                    [celije[i][0], celije[i][1], (celije[i][2]-1)%Ntheta], 
                                    [celije[i][0], celije[i][1], (celije[i][2]+1)%Ntheta]])   
                
    # ------- Formiranje centralne celije ------- #
    # Formiranje liste okolina sa centralnom celijom
    for i in range(len(celije)):
        if okolina[i][0][0] == 0:
            okolina[i][0] = [0, 0, 0]
            

            
    
    return celije, okolina, spoljne_celije, spoljne_celije_indeks
#    okolina1 = []
#    for i in range(len(okolina)):
#        for j in range(len(okolina[i])):
#            for k in range(len(celije)):
#                print(np.array(okolina[i][j]), "AAA") # Treba da se napravi novi niz za okolinu u kom ce clanovi niza predstvaljati redni broj celije iz niza celije a ne samo njenu adresu (poziciju)
#                print(np.array(celije[k]))
#                if np.array(okolina[i][j]) == np.array(celije[k]):
#                    okolina1.append(k)
#            
#    # Formiranje okoline centralne celije
#    okolina_centra = []
#    for i in range(Nphi):
#        for j in range(Ntheta):
#            okolina_centra.append([1, i, j])
#            
#    okolina[0] = okolina_centra
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
#        for j in range(len(okolina[i])):  
#            # Unutrasnje i spoljasnje po r
#            if celije[i][1] == okolina[i][j][1] and celije[i][2] == okolina[i][j][2]: 
#                pomocni_rastojanja.append(abs(r_centri[celije[i][0]]
#                - r_centri[okolina[i][j][0]]))
#                
#                # Slucaj ukoliko je celija spoljna u odnosu na okolinu
#                if celije[i][0] > okolina[i][j][0]:
#                    pomocni_povrsi.append(r[celije[i][0]]**2*dphi*(np.sin(theta[celije[i][2] + 1]) 
#                    - np.sin(theta[celije[i][2]])))
#                    
#                # Slucaj ako je celija unutrasnja u odnosu na okolnu
#                else: 
#                    pomocni_povrsi.append(r[celije[i][0] + 1]**2*dphi*(np.sin(theta[celije[i][2] + 1])
#                    - np.sin(theta[celije[i][2]])))
#            # Levo i desno po phi
#            elif celije[i][0] == okolina[i][j][0] and celije[i][2] == okolina[i][j][2]:
#                pomocni_rastojanja.append(2*r_centri[celije[i][0]]*np.sin(dphi/2)
#                *np.cos(theta_centri[celije[i][2]]))
#                pomocni_povrsi.append((r[celije[i][0]+1]**2 - r[celije[i][0]]**2)/2*dtheta)
#                
#            # Gore i dole po theta
#            elif celije[i][0] == okolina[i][j][0] and celije[i][1] == okolina[i][j][1]:
#                pomocni_rastojanja.append(2*r_centri[celije[i][0]]*np.sin(dtheta/2))
#                
#                # Slucaj ako je celija iznad u odnosu na okolnu
#                if celije[i][2] > okolina[i][j][2]: 
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
#    return (celije, okolina, okolina1, povrsi, zapremine, rastojanja, normale, spoljne_celije, S)

def list_min(lista):
    
    d_min = 1e6

    for i in range(len(lista)):
        for j in range(len(lista[i])):
            if lista[i][j] < d_min:
                d_min = lista[i][j]
    return d_min
                


                
                