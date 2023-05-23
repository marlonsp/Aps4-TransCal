from funcoesTermosol import importa, plota, geraSaida, getDofsIndices
import numpy as np
#nn -> nº de nós
#N -> matriz de nós
#nm -> nº de membros
#Inc -> matriz de incidência
#nc -> nº de cargas
#F -> vetor de carregamento
#nr -> nº de restrições
#R -> vetor de restrições

[nn, N, nm, Inc, nc, F, nr, R] = importa('entrada.xlsx')

def inverte_coluna_linha_matrix(matriz):
    import numpy as np
    matriz_invertida = np.zeros((len(matriz[0]),len(matriz)))
    for i in range(len(matriz)):
        for j in range(len(matriz[0])):
            matriz_invertida[j][i] = matriz[i][j]
    return matriz_invertida
N = inverte_coluna_linha_matrix(N)

def matriz_rigidez(x1,y1,x2,y2,e,a):
    import numpy as np
    import math as mt
    # Calcula o comprimento do elemento
    L = mt.sqrt((x2-x1)**2+(y2-y1)**2)
    # print('L: ', L)
    # Calcula o seno e cosseno
    c = (x2-x1)/L
    s = (y2-y1)/L
    # print('c: ', c)
    # print('s: ', s)
    # calcula matriz de rigidez do elemento, no sistema global
    K = np.array([[c**2, c*s, -c**2, -c*s],[c*s, s**2, -c*s, -s**2],[-c**2, -c*s, c**2, c*s],[-c*s, -s**2, c*s, s**2]])
    K = (e*a/L)*K
    
    return K

# print('N: ', N)
# print('Inc: ', Inc)

K = []
ndof = 2*nn
superK = np.zeros((ndof,ndof))
for i in range(len(Inc)):
    Kx = matriz_rigidez(N[int(Inc[i][0])-1][0], N[int(Inc[i][0])-1][1], N[int(Inc[i][1]-1)][0], N[int(Inc[i][1])-1][1], Inc[i][2], Inc[i][3])
    e_dofs = getDofsIndices([int(Inc[i][0]), int(Inc[i][1])])
    # print('e_dofs: ', e_dofs)
    superK[np.ix_(e_dofs, e_dofs)] += Kx
    K.append(Kx)


print(superK)

# Aplicando as condições de contorno (drop de linhas e colunas com restrições)
for i in range(len(R)-1, -1, -1):
    print(i)
    superK = np.delete(superK, int(R[i][0]), 0)
    superK = np.delete(superK, int(R[i][0]), 1)

print(superK)
