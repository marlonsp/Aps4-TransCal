from funcoesTermosol import importa, plota, geraSaida

#nn -> nº de nós
#N -> matriz de nós
#nm -> nº de membros
#Inc -> matriz de incidência
#nc -> nº de cargas
#F -> vetor de carregamento
#nr -> nº de restrições
#R -> vetor de restrições

[nn, N, nm, Inc, nc, F, nr, R] = importa('entrada.xlsx')

def matriz_rigidez(x1,y1,x2,y2,e,a):
    import numpy as np
    import math as mt
    # Calcula o comprimento do elemento
    L = mt.sqrt((x2-x1)**2+(y2-y1)**2)
    
    # Calcula o seno e cosseno
    c = (x2-x1)/L
    s = (y2-y1)/L
    
    # calcula matriz de rigidez do elemento, no sistema global
    K = np.array([[c**2, c*s, -c**2, -c*s],[c*s, s**2, -c*s, -s**2],[-c**2, -c*s, c**2, c*s],[-c*s, -s**2, c*s, s**2]])
    K = (e*a/L)*K
    
    return K

for i in range(len(N[0])):
    no = [N[0][i], N[1][i]]
    print('Nó: ', no)

K1 = matriz_rigidez(N[0][int(Inc[0][0])+1], N[1][int(Inc[0][0])+1], N[0][1], N[1][1], Inc[0][2], Inc[0][3])
print('K1: ', K1)