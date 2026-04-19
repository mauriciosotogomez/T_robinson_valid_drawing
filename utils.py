import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from itertools import combinations
import io
import csv
import math

#Solve LP
from pulp import * # solve LP

import random
from typing import List, Tuple, Dict, Optional

def get_example(k) :
    examples = [
        [[5, 2, 2, 1, 1],
         [2, 5, 3, 2, 1],
         [2, 3, 5, 4, 1],
         [1, 2, 4, 5, 5],
         [1, 1, 1, 5, 5]],
        
        [[10, 4, 3, 2, 2],
         [4, 10, 6, 4, 2],
         [3, 6, 10, 8, 2],
         [2, 4, 8, 10, 10],
         [2, 2, 2, 10, 10]],

        [[10, 4, 3, 2],  
         [4, 10, 6, 4],  
         [3, 6, 10, 8],  
         [2, 4, 8, 10]],        
         
        [[2, 2, 1, 0, 0],
         [2, 2, 2, 1, 1],
         [1, 2, 2, 2, 1],
         [0, 1, 2, 2, 2],
         [0, 1, 1, 2, 2]],

        [[10,7, 6, 0, 0, 0, 0],
         [7,10, 7, 3, 2, 1, 1],
         [6, 7,10, 7, 2, 2, 1],
         [0, 3, 7, 10,3, 3, 3],
         [0, 2, 2, 3, 10,7, 5],
         [0, 1, 2, 3, 7, 10,6],
         [0, 1, 1, 3, 5, 6,10]],
    
        [[9,4,5,4,1], 
         [4,9,5,4,1], 
         [5,5,9,8,7], 
         [4,4,8,9,7],
         [1,1,7,7,9]]
        ]
    return(examples[k])

def compute_max_closer(distance_matrix) :
    n= len(distance_matrix)
    max_closer=np.zeros((n,n))
    # Compute hanging node distances
    hang_distance = [0]*n
    for i in range(n):
        for j in range(n):
            max_closer[i][j]=i
            increment=1 # left -> right
            if (i>j) :
                increment=-1 # right -> left
            for k in range(i,j,increment):
                if distance_matrix[i][k] < distance_matrix[k][j]:
                    max_closer[i][j] = k
                
    return max_closer

def compute_hang_inter_distance(distance_matrix) :
    n= len(distance_matrix)
    # Compute hanging node distances
    hang_distance = [0]*n
    for i in range(0,n-1):
        hang_distance[i]= (distance_matrix[0][i]+distance_matrix[i][n-1]-distance_matrix[0][n-1])/2    
    # Compute inter distances with previous 
        inter_distance = [0]*n
    for i in range(1,n):
        inter_distance[i]= distance_matrix[i][i-1]-hang_distance[i]-hang_distance[i-1]    
    return(hang_distance,inter_distance)
    
def plot_caterpillar(distance_matrix, nodelabel, is_correct=True) :
    n= len(distance_matrix)
    tol=0.000001
    xhop = 2
    yhop = -2
    pos = {}
    edge_labels={}
    edges= []
    hang_distance,inter_distance=compute_hang_inter_distance(distance_matrix)
    edge_color='red'
    if is_correct:
        edge_color='black'
    
    # Create edges
    prev_auxnode=nodelabel[0] # previous node in the backbone    
    for i in range(0,n) :
        # put the node in the backbone
        pos[nodelabel[i]]=(i*xhop,0)
        auxnode=nodelabel[i]
        # if hang is positive create the auxnode an the hang edge
        if hang_distance[i]>tol :
            auxnode="aux"+str(i)
            pos[auxnode]=(i*xhop,0)
            pos[nodelabel[i]]=(i*xhop,yhop)
            # create hang edge
            edge=[auxnode, nodelabel[i]]
            edges.append(edge)
            edge_labels[tuple(edge)]=str(round(hang_distance[i],2))
        # create the edge in the backbone (with the previous auxnode)    
        if i>0 :
            edge=[auxnode, prev_auxnode]
            edges.append(edge)
            edge_labels[tuple(edge)]=str(round(inter_distance[i],2))
            prev_auxnode=auxnode

    # Create graph
    G = nx.Graph()
    G.add_edges_from(edges)
    #nx.write_network_text(G, vertical_chains=True, sources=["S0"])
    plt.figure()
    # Plot nodes
    nodes=nx.draw(
        G,
        pos,
        nodelist=nodelabel,
        edgecolors=edge_color,
        node_size=700,
        node_color='gainsboro',
        alpha=1.0,
        labels={node: node for node in nodelabel}
        )
    # Plot edges
    nx.draw_networkx_edge_labels(
        G, pos,
        edge_labels=edge_labels,
        font_color=edge_color,
        rotate=False,
        font_size=14
    )
    plt.axis('off')    # Create edges
    prev_auxnode=nodelabel[0] # previous node in the backbone    
    for i in range(0,n) :
        # put the node in the backbone
        pos[nodelabel[i]]=(i*xhop,0)
        auxnode=nodelabel[i]
        # if hang is positive create the auxnode an the hang edge
        if hang_distance[i]>tol :
            auxnode="aux"+str(i)
            pos[auxnode]=(i*xhop,0)
            pos[nodelabel[i]]=(i*xhop,yhop)
            # create hang edge
            edge=[auxnode, nodelabel[i]]
            edges.append(edge)
            edge_labels[tuple(edge)]=str(round(hang_distance[i],2))
        # create the edge in the backbone (with the previous auxnode)    
        if i>0 :
            edge=[auxnode, prev_auxnode]
            edges.append(edge)
            edge_labels[tuple(edge)]=str(round(inter_distance[i],2))
            prev_auxnode=auxnode
    plt.axis('equal')
    plt.show()

def plot_text_caterpillar(distance_matrix, nodelabel) :
    print("\n Tree structure")    
    n= len(distance_matrix)
    tol=0.000001
    hang_distance,inter_distance=compute_hang_inter_distance(distance_matrix)
    f = io.StringIO()

    print("("+str(nodelabel[0])+")", file=f)
    for i in range(1,n) :
        print(" |", file=f)
        print(str(round(inter_distance[i],2)), file=f)
        # put the node in the backbone
        print(" |", file=f)
        # if hang is positive create the auxnode an the hang edge
        if hang_distance[i]>tol :
            print(" |- "+str(round(hang_distance[i],2))+" - ("+str(nodelabel[i])+")", file=f)
        else : 
            print("("+str(nodelabel[i])+")", file=f)
    # Print in stand.out
    out = f.getvalue()
    print(out)
    f.close()
    

def check_solution (distance_matrix, outmatrix):
    is_correct=True
    n=len(distance_matrix)
    for j in range(n):
        for i in range(n): # center
            for k in range(n):
                if (distance_matrix[i][j] > distance_matrix[i][k]) and  (outmatrix[i][j] < outmatrix[i][k]):  
                    is_correct=False
                    print(f'delta_sim:{distance_matrix[i][j]} - {distance_matrix[i][k]}, delta_dist:{outmatrix[i][j]} - {outmatrix[i][k]}')
    return(is_correct)


def all_kernel_bases(matrix, output_file):
    rows, cols = matrix.shape
    
    # Get all combinations of column indices with cardinality equal to the number of rows
    column_combinations = combinations(range(cols), rows)

    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for column_indices in column_combinations:
            subset_cols = matrix[:, column_indices]
            remaining_cols = np.delete(matrix, column_indices, axis=1)
        
            # Check if the submatrix is invertible
            if np.linalg.det(subset_cols) != 0:        
                # Compute the inverse of the submatrix
                subset_cols_inv = np.linalg.inv(subset_cols)
        
                # Compute the multiplication of the inverse with the remaining columns
                result = -np.dot(subset_cols_inv, remaining_cols)
                
                # Save the result to the CSV file
                column_indices_str = ','.join(map(str, column_indices))
                for row in result:
                    writer.writerow([column_indices_str] + row.tolist())
                print(column_indices_str)
                writer.writerow("")

def generate_strict_center_matrix(n, seed=None):
    if seed is not None:
        np.random.seed(seed)

    M = np.zeros((n, n), dtype=int)

    # Main diagonal
    for i in range(n):
        M[i, i] = i

    # Secondary upper diagonal
    for i in range(n - 1):
        M[i, i + 1] = i

    # Upper half
    for d in range(2, n):  # Range using the distance from the main diagonal
        for i in range(n - d):
            j = i + d
            a = M[i, j - 1]
            b = M[i + 1, j]
            low = min(a, b)
            high = max(a, b)
            M[i, j] = np.random.randint(low, high + 1)

    # Parte inferior izquierda (simétrica + 1)
    for i in range(n):
        for j in range(i):
            M[i, j] = M[j, i] + 1

    # Lower Secondary diagonal 
    for i in range(1,n):
        M[i, i-1] = i

    # # Lower part (other than the secondary diagonal)
    # for d in range(2, n):  # distance from the diagonal
    #     for i in range(n - d):
    #         j = i + d
    #         a = M[j - 1, i]
    #         b = M[j, i + 1]
    #         low = max(a, M[i,j]+1)
    #         high = b
    #         M[j, i] = np.random.randint(low, high + 1)
            
    return M
                

def generate_center_matrix(n, seed=None):
    if seed is not None:
        np.random.seed(seed)

    M = np.zeros((n, n), dtype=int)

    # Main diagonal
    for i in range(n):
        M[i, i] = i

    # Secondary upper diagonal
    for i in range(n - 1):
        M[i, i + 1] = i

    # Upper half
    for d in range(2, n):  # Range using the distance from the main diagonal
        for i in range(n - d):
            j = i + d
            a = M[i, j - 1]
            b = M[i + 1, j]
            low = min(a, b)
            high = max(a, b)
            M[i, j] = np.random.randint(low, high + 1)

    # Lower Secondary diagonal 
    for i in range(1,n):
        M[i, i-1] = i

    # Lower part (other than the secondary diagonal)
    for d in range(2, n):  # distance from the diagonal
        for i in range(n - d):
             j = i + d
             a = M[j - 1, i]
             b = M[j, i + 1]
             low = max(a, M[i,j]+1)
             high = b
             M[j, i] = np.random.randint(low, high + 1)

    #print(M)        
    return M

def generate_simple_matrix(n, seed=None):
    if seed is not None:
        np.random.seed(seed)

    M = np.zeros((n, n), dtype=int)

    # Diagonal principal
    for i in range(n):
        M[i, i] = i

    # Diagonal justo encima de la principal
    for i in range(n - 1):
        M[i, i + 1] = i

    # Parte superior derecha (más allá de la diagonal justo encima)
    for d in range(2, n):  # distancia desde la diagonal principal
        for i in range(n - d):
            j = i + d
            M[i, j] = math.floor((j+i)/2)

    # Parte inferior izquierda (simétrica + 1)
    for i in range(n):
        for j in range(i):
            M[i, j] = M[j, i] + 1

    return M


# ------------------------------------------------------------
# 1) Generaci'f3n de matriz M
# ------------------------------------------------------------
def generar_matriz(n: int, seed: Optional[int] = None) -> List[List[int]]:
    """
    Genera una matriz M de tama'f1o n'd7n cumpliendo:
      1) M(i,i) = i
      2) M(i,i+1) = i
      3) Para j - i >= 2: M(i,j) es entero aleatorio en [ M(i,j-1), M(i+1,j) ]
      4) Para i > j: M(j,i) = M(i,j) + 1

    Notaci'f3n matem'e1tica 1-indexada en esta especificaci'f3n; implementaci'f3n 0-indexada.
    """
    if n <= 0:
        raise ValueError("n debe ser un entero positivo")

    if seed is not None:
        random.seed(seed)

    M: List[List[Optional[int]]] = [[None for _ in range(n)] for _ in range(n)]

    # Diagonal principal: M(i,i) = i
    for i in range(n):
        M[i][i] = i + 1

    # Segunda diagonal: M(i,i+1) = i
    for i in range(n - 1):
        M[i][i + 1] = i + 1

    # Diagonales superiores restantes (j - i >= 2)
    for d in range(2, n):  # d = j - i
        for i in range(0, n - d):
            j = i + d
            low = M[i][j - 1]
            up  = M[i + 1][j]
            if low is None or up is None:
                raise RuntimeError("L'edmites no definidos al construir M.")
            if low > up:
                raise ValueError(f"Intervalo vac'edo al fijar M({i+1},{j+1}): [{low}, {up}]")
            M[i][j] = random.randint(low, up)

    # Tri'e1ngulo inferior: M(j,i) = M(i,j) + 1 para i < j
    for i in range(n):
        for j in range(i + 1, n):
            M[j][i] = M[i][j] + 1

    # type: ignore (ya no hay None)
    return [list(map(int, fila)) for fila in M]  # convertir a int por seguridad


def generate_robinson_from_center_matrix(M) :
    
    #######################
    # Problem description #
    #######################
    n = M.shape[1]
    ROWS = COLS = range(n)

    ###############################
    #### PROBLEM FOR MATRIX ENTRIES ###
    ###############################
    problem = LpProblem("Robinson_from_Center_Matrix", LpMinimize)
    # Define decision variables
    distance = LpVariable.dicts("distance", (ROWS, COLS), cat='Continuous', lowBound=0)
    # We do not define an objective function since none is needed

    #problem += distance[0][n-1] == 100
    
    for i in ROWS :
        problem += distance[i][i] == 0
        for j in range(i) :
            problem += distance[i][j] - distance[j][i] == 0  
        
    # UPPER HALF i<j
    for j in range(1,n):
        for i in range(j):
            # Robinson condition
            problem += distance[i][j-1] - distance[i][j] <= 0
            problem += distance[i+1][j] - distance[i][j] <= 0
            #Left Center 
            left_center= M[i][j]
            problem += distance[i][left_center] - distance[left_center][j] <= -1            
            problem += distance[left_center+1][j] - distance[i][left_center+1] <= 0


    # LOWER HALF i>j
    for i in range(1,n):
        for j in range(i):
            # Robinson condition
            problem += distance[i][j+1] - distance[i][j] <= 0
            problem += distance[i-1][j] - distance[i][j] <= 0
            #Right Center 
            right_center= M[i][j]
            problem += distance[right_center][i] - distance[j][right_center] <= -1            
            problem += distance[j][right_center-1] - distance[right_center-1][i] <= 0

    # The problem data is written to an .lp file
    problem.writeLP("robinson_from_center.lp")
    # Solve the aux_problem
    problem.solve(PULP_CBC_CMD(msg=0))

    has_solution = (LpStatus[problem.status]=="Optimal")
    # print(f"Has solution {has_solution}")

    # Compute output matrix
    outmatrix=np.zeros((n,n))
    
    for r in ROWS:
        for c in COLS:
            outmatrix[r][c] = distance[r][c].varValue
   
    if has_solution :
        return outmatrix
    return False

# ------------------------------------------------------------
# 2) Construcci'f3n del grafo y verificaci'f3n DAG
# ------------------------------------------------------------
def construir_grafo(M: List[List[int]]) -> Dict[Tuple[int, int], List[Tuple[int, int]]]:
    """
    Construye el grafo dirigido a partir de M seg'fan:
      - V = {(i,j) : 1 <= i <= j <= n}
      - Regla 1: (i,j) -> (i,j+1)           para 1 <= i <= j < n
      - Regla 2: (i,j) -> (i-1,j)           para 1 <  i <= j <= n
      - Regla 3: (i, M(i,j)) -> (M(i,j), j) para 1 <= i <  j <= n   [corregida: solo i<j]
      - Regla 4: (M(i,j)+1, j) -> (i, M(i,j)+1) para 1 <= i < j <= n [corregida: solo i<j]

    Usa 1-indexado en los nombres de v'e9rtices.
    """
    n = len(M)
    V = {(i, j) for i in range(1, n + 1) for j in range(i, n + 1)}
    adj: Dict[Tuple[int, int], List[Tuple[int, int]]] = {v: [] for v in V}

    def add_edge(u: Tuple[int, int], v: Tuple[int, int]) -> None:
        if u in V and v in V:
            adj[u].append(v)

    # Regla 1
    for i in range(1, n + 1):
        for j in range(i, n):
            add_edge((i, j), (i, j + 1))

    # Regla 2
    for i in range(2, n + 1):
        for j in range(i, n + 1):
            add_edge((i, j), (i - 1, j))

    # Reglas 3 y 4 (solo si i < j)
    for i in range(1, n + 1):
        for j in range(i, n + 1):
            if i < j:
                mij = M[i-1][j-1]
                # Regla 3
                add_edge((i, mij), (mij, j))
                # Regla 4
                add_edge((mij + 1, j), (i, mij + 1))

    return adj


def es_DAG(M: List[List[int]]) -> str:
    """
    Devuelve "YES" si el grafo dirigido inducido por M es ac'edclico (DAG),
    y "NO" en caso contrario.
    """
    n = len(M)
    V = {(i, j) for i in range(1, n + 1) for j in range(i, n + 1)}
    adj = construir_grafo(M)

    WHITE, GRAY, BLACK = 0, 1, 2
    color: Dict[Tuple[int, int], int] = {v: WHITE for v in V}

    def dfs(v: Tuple[int, int]) -> bool:
        color[v] = GRAY
        for u in adj[v]:
            if color[u] == GRAY:
                return False  # ciclo
            if color[u] == WHITE and not dfs(u):
                return False
        color[v] = BLACK
        return True

    for v in V:
        if color[v] == WHITE:
            if not dfs(v):
                return "NO"
    return "YES"
