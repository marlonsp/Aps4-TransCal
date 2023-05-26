from funcoesTermosol import *
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

#matriz de rigidez global
N = np.transpose(N)
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

# Aplicando as condições de contorno (drop de linhas com restrições)
newF = np.delete(F, list(map(int, (R[:,0]))), 0).squeeze()

#iterações e tolerância
ite = 1000
tol = 1e-8

#deslocamentos globais
deslocamentos_global = gauss_seidel(ite, tol, superK, newF)[0]

#preencher o vetor de deslocamentos_global com os valores de deslocamentos_global nulos
for i in range(len(R)):
    deslocamentos_global = np.insert(deslocamentos_global, int(R[i][0]), 0)

#dividir o vetor de deslocamentos_global em listas de 2 em 2
deslocamentos = [deslocamentos_global[i:i+2] for i in range(0, len(deslocamentos_global), 2)]

#deformações específicas
deformações = []
for i in range(len(Inc)):
    deformações.append(deformação_específica_elemento(N[int(Inc[i][0])-1][0], N[int(Inc[i][0])-1][1], N[int(Inc[i][1]-1)][0], N[int(Inc[i][1])-1][1],deslocamentos[int(Inc[i][0])-1],deslocamentos[int(Inc[i][1])-1]))

#tensões
tensões = []
for i in range(len(Inc)):
    tensões.append(tensão_elemento(N[int(Inc[i][0])-1][0], N[int(Inc[i][0])-1][1], N[int(Inc[i][1]-1)][0], N[int(Inc[i][1])-1][1],Inc[i][2], deslocamentos[int(Inc[i][0])-1],deslocamentos[int(Inc[i][1])-1]))


#reacoes de apoio
reacoes = np.dot(superK_antes, deslocamentos_global)
reacoes = reacoes[list(map(int, (R[:,0])))] 

#Forças internas
areas = Inc[:,3]
forças = areas*tensões

#plotar o gráfico antes e depois da deformação
N = np.transpose(N)
deslocamentos = np.transpose(deslocamentos)
N_pos = N + (deslocamentos*10000)

plota(N, Inc)
plota(N_pos+deslocamentos, Inc)

with open('output.txt', 'w') as file:
    file.write("Reacoes de apoio [N]\n" + str(reacoes.reshape(-1, 1)) + "\n\n")
    file.write("Deslocamentos [m]\n" + str(deslocamentos_global.reshape(-1, 1)) + "\n\n")
    file.write("Deformacoes []\n" + str(np.array(deformações).reshape(-1, 1)) + "\n\n")
    file.write("Forcas internas [N]\n" + str(np.array(forças).reshape(-1, 1)) + "\n\n")
    file.write("Tensoes internas [Pa]\n" + str(np.array(tensões).reshape(-1, 1)) + "\n")

