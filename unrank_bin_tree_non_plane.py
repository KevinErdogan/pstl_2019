
import time


class Arbre: 
        
    def __init__(self, cle, F):
        self.cle = cle
        self.fils = F
       
    def dot_ch(self):
        dot = ''
        index = self.cle
        for f in self.fils:
            if f.cle>0:
                dot += str(index) + ' -> ' + str(f.cle) + ';\n'
        for f in self.fils:
            if f.cle>0:
                dot += f.dot_ch()
        return dot

def gener_feuille(cle):
    return Arbre(cle, [])

def gener_noeud(cle, F):
    return Arbre(cle, F)

def affiche(A):
    dot = 'digraph { node [shape=point]; edge [arrowhead=none];\n'
    return dot + A.dot_ch() + '}'
#%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%


# stocke dans la variable T le vecteur du nombre d'arbres
# binaires de chaque taille
T = [0,1]
def enum(n):
    global T

    if len(T)<n:
        for i in range(len(T),n+1):
            if i%2 == 1:
                t = sum([T[k]*T[i-1-k]/2 for k in range(1,i-1,2)])
                if i%2 == 1:
                    t += T[(i-1)//2]/2
            else:
                t = 0
            T.append( int(t) )
    return T[n]


Etiq = set()
# construit l'arbre de rang r (entre 0 et Tn-1) parmi les arbres de taille n
def gen(n, r):
    global T, Etiq

    if n==1:
        t = time.time()
        while t in Etiq:
            t = time.time()
        Etiq.add(t)
        return gener_feuille(t)

    k = 1
    while r>=0:
        temp = enum(k)*enum(n-1-k)
        if n%2 == 1 and k == (n-1)//2:
            temp += enum((n-1)//2)/2
        r -= temp
        k += 2
    k -= 2
    r += temp
    print(k,r)
    t = time.time()
    while t in Etiq:
        t = time.time()
    Etiq.add(t)

    if k < n-1-k:
        f1 = gen(k, r % enum(k))
        f2 = gen(n-1-k, r // enum(k))
    elif r < enum((n-1)//2):
        f1 = gen((n-1)//2, r % enum((n-1)//2))
        f2 = gen((n-1)//2, r % enum((n-1)//2))
    else:
        r -= enum((n-1)//2)
        f1 = gen((n-1)//2, r % enum((n-1)//2))
        f2 = gen((n-1)//2, r // enum((n-1)//2))
        
    return gener_noeud(t, [f1, f2])


A = gen(21, 101)
## copier la chaine suivante dans un fichier : test.dot
## puis lancer graphviz :
## dot -Tpdf test.dot -o test.pdf

print(affiche(A))

