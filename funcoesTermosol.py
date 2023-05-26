    # -*- coding: utf-8 -*-
"""
A funcao 'plota' produz um gráfico da estrutura definida pela matriz de nos N 
e pela incidencia Inc.

Sugestao de uso:

from funcoesTermosol import plota
plota(N,Inc)
-------------------------------------------------------------------------------
A funcao 'importa' retorna o numero de nos [nn], a matriz dos nos [N], o numero
de membros [nm], a matriz de incidencia [Inc], o numero de cargas [nc], o vetor
carregamento [F], o numero de restricoes [nr] e o vetor de restricoes [R] 
contidos no arquivo de entrada.

Sugestao de uso:
    
from funcoesTermosol import importa
[nn,N,nm,Inc,nc,F,nr,R] = importa('entrada.xlsx')
-------------------------------------------------------------------------------
A funcao 'geraSaida' cria um arquivo nome.txt contendo as reacoes de apoio Ft, 
deslocamentos Ut, deformacoes Epsi, forcas Fi e tensoes Ti internas. 
As entradas devem ser vetores coluna.

Sugestao de uso:
    
from funcoesTermosol import geraSaida
geraSaida(nome,Ft,Ut,Epsi,Fi,Ti)
-------------------------------------------------------------------------------

"""
def plota(N,Inc):
    # Numero de membros
    nm = len(Inc[:,0])
    
    import matplotlib as mpl
    import matplotlib.pyplot as plt

#    plt.show()
    fig = plt.figure()
    # Passa por todos os membros
    for i in range(nm):
        
        # encontra no inicial [n1] e final [n2] 
        n1 = int(Inc[i,0])
        n2 = int(Inc[i,1])        

        plt.plot([N[0,n1-1],N[0,n2-1]],[N[1,n1-1],N[1,n2-1]],color='r',linewidth=3)


    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.grid(True)
    plt.axis('equal')
    plt.show()
    
def importa(entradaNome):
    
    import numpy as np
    import xlrd
    
    arquivo = xlrd.open_workbook(entradaNome)
    
    ################################################## Ler os nos
    nos = arquivo.sheet_by_name('Nos')
    
    # Numero de nos
    nn = int(nos.cell(1,3).value)
                 
    # Matriz dos nós
    N = np.zeros((2,nn))
    
    for c in range(nn):
        N[0,c] = nos.cell(c+1,0).value
        N[1,c] = nos.cell(c+1,1).value
    
    ################################################## Ler a incidencia
    incid = arquivo.sheet_by_name('Incidencia')
    
    # Numero de membros
    nm = int(incid.cell(1,5).value)
                 
    # Matriz de incidencia
    Inc = np.zeros((nm,4))
    
    for c in range(nm):
        Inc[c,0] = int(incid.cell(c+1,0).value)
        Inc[c,1] = int(incid.cell(c+1,1).value)
        Inc[c,2] = incid.cell(c+1,2).value
        Inc[c,3] = incid.cell(c+1,3).value
    
    ################################################## Ler as cargas
    carg = arquivo.sheet_by_name('Carregamento')
    
    # Numero de cargas
    nc = int(carg.cell(1,4).value)
                 
    # Vetor carregamento
    F = np.zeros((nn*2,1))
    
    for c in range(nc):
        no = carg.cell(c+1,0).value
        xouy = carg.cell(c+1,1).value
        GDL = int(no*2-(2-xouy)) 
        F[GDL-1,0] = carg.cell(c+1,2).value
         
    ################################################## Ler restricoes
    restr = arquivo.sheet_by_name('Restricao')
    
    # Numero de restricoes
    nr = int(restr.cell(1,3).value)
                 
    # Vetor com os graus de liberdade restritos
    R = np.zeros((nr,1))
    
    for c in range(nr):
        no = restr.cell(c+1,0).value
        xouy = restr.cell(c+1,1).value
        GDL = no*2-(2-xouy) 
        R[c,0] = GDL-1

    return nn,N,nm,Inc,nc,F,nr,R

def geraSaida(nome,Ft,Ut,Epsi,Fi,Ti):
    nome = nome + '.txt'
    f = open("saida.txt","w+")
    f.write('Reacoes de apoio [N]\n')
    f.write(str(Ft))
    f.write('\n\nDeslocamentos [m]\n')
    f.write(str(Ut))
    f.write('\n\nDeformacoes []\n')
    f.write(str(Epsi))
    f.write('\n\nForcas internas [N]\n')
    f.write(str(Fi))
    f.write('\n\nTensoes internas [Pa]\n')
    f.write(str(Ti))
    f.close()

def getDofsIndices(array):
    # Retorna os indices dos graus de liberdade de um vetor de nos
    '''
    Args:
        array: vetor de nos
    Returns:
        indices: vetor com os indices dos graus de liberdade
    '''
    import numpy as np
    indices = np.zeros(len(array)*2,dtype=int)
    for i in range(len(array)):
        indices[i*2] = array[i]*2-2
        indices[i*2+1] = array[i]*2-1
    return indices

def matriz_rigidez(x1,y1,x2,y2,e,a):
    # Calcula a matriz de rigidez de um elemento de barra
    '''
    Args:
        x1,y1,x2,y2: coordenadas dos nos do elemento
        e: modulo de elasticidade
        a: area da secao transversal
    Returns:
        K: matriz de rigidez do elemento
    '''
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
def gauss_seidel(ite, tol, K, F): 
    # Metodo de Gauss-Seidel para resolver sistemas lineares
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
    import numpy as np
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
    # Metodo de Jacobi para resolver sistemas lineares
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
    import numpy as np
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

def deformação_específica_elemento(x1, y1, x2, y2, u,v):
    # Calcula a deformação específica
    '''
    Args:
        x1,y1,x2,y2: coordenadas dos nos do elemento
        u,v: deslocamentos dos nos do elemento
    Returns:
        e: deformação específica
    '''
    import math as mt
    # Calcula o comprimento do elemento
    L = mt.sqrt((x2-x1)**2+(y2-y1)**2)
    # Calcula o seno e cosseno
    c = (x2-x1)/L
    s = (y2-y1)/L
    # Calcula a deformação específica
    e = (c*(v[0]-u[0])+s*(v[1]-u[1]))/L
    
    return e

def tensão_elemento(x1, y1, x2, y2, e, u,v):
    # Calcula a tensão
    '''
    Args:
        x1,y1,x2,y2: coordenadas dos nos do elemento
        e: modulo de elasticidade
        u,v: deslocamentos dos nos do elemento
    Returns:
        t: tensão
    '''
    import math as mt
    # Calcula o comprimento do elemento
    L = mt.sqrt((x2-x1)**2+(y2-y1)**2)
    # Calcula o seno e cosseno
    c = (x2-x1)/L
    s = (y2-y1)/L
    # Calcula a tensão
    t = (e/L)*(c*(v[0]-u[0])+s*(v[1]-u[1]))
    
    return t