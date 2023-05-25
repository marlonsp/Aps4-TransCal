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
superK_antes = np.zeros((ndof,ndof))
for i in range(len(Inc)):
    Kx = matriz_rigidez(N[int(Inc[i][0])-1][0], N[int(Inc[i][0])-1][1], N[int(Inc[i][1]-1)][0], N[int(Inc[i][1])-1][1], Inc[i][2], Inc[i][3])
    e_dofs = getDofsIndices([int(Inc[i][0]), int(Inc[i][1])])
    # print('e_dofs: ', e_dofs)
    superK_antes[np.ix_(e_dofs, e_dofs)] += Kx
    K.append(Kx)


# Aplicando as condições de contorno (drop de linhas e colunas com restrições)
superK = np.delete(superK_antes, list(map(int, (R[:,0]))), 0)
superK = np.delete(superK, list(map(int, (R[:,0]))), 1)
print("superK:\n",superK)

# newF = np.delete(F, list(map(int, (R[:,0]))), 0)
newF = np.delete(F, list(map(int, (R[:,0]))), 0).squeeze()
print("newF:\n",newF)

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
    variaveis = {}

    for i in range(len(F)):
        variaveis["x"+str(i+1)] = 0
    
    U = np.array([])
    for i in range(len(F)):
        U = np.append(U, variaveis["x"+str(i+1)])
    
    for i in range(ite):
        for j in range(len(F)):
            variaveis["x"+str(j+1)] = F[j]
            for k in range(len(F)):
                if k != j:
                    variaveis["x"+str(j+1)] -= K[j,k]*variaveis["x"+str(k+1)]

            variaveis["x"+str(j+1)] /= K[j,j]
            U[j] = variaveis["x"+str(j+1)]

        ei = np.linalg.norm(np.dot(K, U) - F)/np.linalg.norm(F)
        if ei < tol:
            break

    return U, ei

def jacobi(ite, tol, K, F): 
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
    # Metodo de Jacobi
    variaveis_ant = {}
    variaveis_pos = {}
    for i in range(len(F)):
        variaveis_ant["x"+str(i+1)] = 0
        variaveis_pos["x"+str(i+1)] = 0

    U = np.array([])
    for i in range(len(F)):
        U = np.append(U, variaveis_ant["x"+str(i+1)])

    for i in range(ite):
        for j in range(len(F)):
            variaveis_pos["x"+str(j+1)] = F[j]
            for k in range(len(F)):
                if k != j:
                    variaveis_pos["x"+str(j+1)] -= K[j,k]*variaveis_ant["x"+str(k+1)]

            variaveis_pos["x"+str(j+1)] /= K[j,j]
            U[j] = variaveis_pos["x"+str(j+1)]

            if j == len(F)-1:
                for k in range(len(F)):
                    variaveis_ant["x"+str(k+1)] = variaveis_pos["x"+str(k+1)]

        ei = np.linalg.norm(np.dot(K, U) - F)/np.linalg.norm(F)
        if ei < tol:
            break
    
    return U, ei

ite = 100
tol = 1e-6

print("Gauss-Seidel:\n", gauss_seidel(ite, tol, superK, newF))
print("Jacobi:\n", jacobi(ite, tol, superK, newF))

deslocamentos_global = gauss_seidel(ite, tol, superK, newF)[0]

#preencher o vetor de deslocamentos_global com os valores de deslocamentos_global nulos
for i in range(len(R)):
    deslocamentos_global = np.insert(deslocamentos_global, int(R[i][0]), 0)
#dividir o vetor de deslocamentos_global em listas de 2 em 2
deslocamentos = [deslocamentos_global[i:i+2] for i in range(0, len(deslocamentos_global), 2)]

def deformação_específica_elemento(x1, y1, x2, y2, u,v):
    import math as mt
    # Calcula o comprimento do elemento
    L = mt.sqrt((x2-x1)**2+(y2-y1)**2)
    # Calcula o seno e cosseno
    c = (x2-x1)/L
    s = (y2-y1)/L
    # Calcula a deformação específica
    e = (c*(v[0]-u[0])+s*(v[1]-u[1]))/L
    
    return e



deformações = []
for i in range(len(Inc)):
    deformações.append(deformação_específica_elemento(N[int(Inc[i][0])-1][0], N[int(Inc[i][0])-1][1], N[int(Inc[i][1]-1)][0], N[int(Inc[i][1])-1][1],deslocamentos[int(Inc[i][0])-1],deslocamentos[int(Inc[i][1])-1]))

print("Deformações específicas:\n", deformações)

def tensão_elemento(x1, y1, x2, y2, e, u,v):
    # Calcula a tensão
    import math as mt
    # Calcula o comprimento do elemento
    L = mt.sqrt((x2-x1)**2+(y2-y1)**2)
    # Calcula o seno e cosseno
    c = (x2-x1)/L
    s = (y2-y1)/L
    # Calcula a tensão
    t = (e/L)*(c*(v[0]-u[0])+s*(v[1]-u[1]))
    
    return t

tensões = []
for i in range(len(Inc)):
    tensões.append(tensão_elemento(N[int(Inc[i][0])-1][0], N[int(Inc[i][0])-1][1], N[int(Inc[i][1]-1)][0], N[int(Inc[i][1])-1][1],Inc[i][2], deslocamentos[int(Inc[i][0])-1],deslocamentos[int(Inc[i][1])-1]))

print("Tensões:\n", tensões)

# print f
reacoes = np.dot(superK_antes, deslocamentos_global)
print("Reações de apoio:\n", reacoes)