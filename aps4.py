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
N = np.transpose(N)
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


# print(superK)
# Aplicando as condições de contorno (drop de linhas e colunas com restrições)
superK = np.delete(superK, list(map(int, (R[:,0]))), 0)
superK = np.delete(superK, list(map(int, (R[:,0]))), 1)

newF = np.delete(F, list(map(int, (R[:,0]))), 0)
# print(superK)

def gauss_seidel(ite, tol, K, F): 
    """
    Args: 
        ite: numero maximo de iteracoes
        tol: tolerancia
        K: matriz de rigidez
        F: vetor de forcas
    Returns:
        u: vetor de deslocamentos
        ei: erro maximo
    """
    # Metodo de Gauss-Seidel
    x1, x2, x3 = 0, 0, 0
    U = np.array([x1, x2, x3])
    
    for i in range(ite):
        x1 = (F[0] - K[0,1]*x2 - K[0,2]*x3)/K[0,0]
        x2 = (F[1] - K[1,0]*x1 - K[1,2]*x3)/K[1,1]
        x3 = (F[2] - K[2,0]*x1 - K[2,1]*x2)/K[2,2]
        U = np.array([x1, x2, x3])
        ei = np.linalg.norm(np.dot(K, U) - F)/np.linalg.norm(F)
        if ei < tol:
            break

    return U, ei