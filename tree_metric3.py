import numpy as np # represent
from pulp import * # solve LP

# compute the tree 
from skbio import DistanceMatrix 
from skbio.tree import nj

# plot the tree
from Bio import Phylo
from io import StringIO

# arguments
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-e", "--example number",
                   dest="n_example", default=0,
                   help="use this example from the list ")
args = parser.parse_args()

n_example = 0
n_example = int(args.n_example)

all_labels=["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]

all_input_matrix = [[[5, 2, 2, 1, 1],
                         [2, 5, 3, 2, 1],
                         [2, 3, 5, 4, 1],
                         [1, 2, 4, 5, 5],
                         [1, 1, 1, 5, 5]],
                        [[10, 4, 3, 2, 2],
                         [4, 10, 6, 4, 2],
                         [3, 6, 10, 8, 2],
                         [2, 4, 8, 10, 10],
                         [2, 2, 2, 10, 10]],
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
                         [3,3,8,9,7],
                         [1,1,7,7,9]]]

input_matrix = np.array(all_input_matrix[n_example])

n = len(input_matrix)
ROWS = COLS = range(n)
nodelabel = all_labels[:n]
nodelabel = ["S"+str(i) for i in range(n)]


# The boxes list is created, with the row and column index of each square in each box
four_points = [
    (i,j,k,l) 
    for l in range(3,n)
    for k in range(2,l)
    for j in range(1,k)
    for i in range(j)
]

three_points = [
    (i,j,k) 
    for k in range(n)
    for j in range(n)
    for i in range(n)
]

two_points = [
    (i,j) 
    for j in range(1,n)
    for i in range(j)
]

# CREATE THE PROBLEM  ##################################
problem = LpProblem("Has a tree diagram Problem", LpMinimize)

# Define decision variables
distance = LpVariable.dicts("distance", (ROWS, COLS), cat='Continuous', lowBound=0)
max = LpVariable.dicts("max", (four_points), cat='Continuous', lowBound=0)

# We do not define an objective function since none is needed
problem += lpSum([max[q] for q in four_points] )

# Constraints ensuring the 4-point condition
#   i  j   k   l
# i * Dij Dik Dil
# j    *  Djk Djl
# k        *  Dkl
# l            *
# Assuming d(i,k)+d(j,l) is the middle sum
# d(i,j)+d(k,l) <= d(i,k)+d(j,l)
# d(i,l)+d(j,k) <=  d(i,k)+d(j,l)

for (i,j,k,l) in four_points :
    problem += distance[i][k] + distance[j][l] == max[(i,j,k,l)]
    problem += distance[i][j] + distance[k][l] <= max[(i,j,k,l)]
    problem += distance[i][l] + distance[j][k] <= max[(i,j,k,l)] 

# Constraints ensuring the
#   - 3-point condition (consistent with the matrix)
#   - triangular inequality (metric)
for (i,j,k) in three_points:
    if input_matrix[i][j] > input_matrix[i][k] :
        problem += distance[i][j] - distance[i][k] <= -1
    problem += distance[i][j] + distance[j][k] - distance[i][k] >= 0        

# Symmetry 2-points condition (metric)
for (i,j) in two_points:
    problem += distance[i][j] - distance[j][i] == 0

# The problem data is written to an .lp file
problem.writeLP("tree_metric.lp")

# Solve the problem
problem.solve()

# The status and value of the solution is printed to the screen
print(f'Status: {LpStatus[problem.status]}')
print(f'Value : {value(problem.objective)}')

#PRINT INPUT/OUTPUT INFO ##################################

# Print input matrix
print(f"\n Input similarity matrix \n{input_matrix}")

# Print output matrix
outmatrix=np.zeros((n,n))
for r in ROWS:
    for c in COLS:
       outmatrix[r][c] = value(distance[r][c])

print(f"\n Output distance matrix \n{outmatrix}")

# A file called matrixout.txt is created/overwritten for writing to
matrixout = open("tree_metric_out.txt", "w")

# The solution is written to the matrixout.txt file
for c in COLS:
    matrixout.write(f"{nodelabel[c]}, " )
matrixout.write("\n")
for r in ROWS:
    matrixout.write(f"{nodelabel[r]}")
    for c in COLS:
        matrixout.write(f", {value(distance[r][c])}")
    matrixout.write("\n")
matrixout.close()

# Print the tree
dm = DistanceMatrix(outmatrix, nodelabel)
tree = nj(dm)
newick_str = nj(dm, result_constructor=str)

print("\n Structure of the tree")
print(tree.ascii_art())

print("\n Newick notation of the tree")
print(newick_str)

for (i,j,k,l) in four_points :    
     # if  value(distance[i][j]) + value(distance[k][l]) != value(max[(i,j,k,l)]) and value(distance[i][k]) + value(distance[j][l]) != value(max[(i,j,k,l)]) and value(distance[i][l]) + value(distance[j][k]) != value(max[(i,j,k,l)]) :
    print(f"[{(i,j,k,l)}] max={value(max[(i,j,k,l)])}: d({i},{j})+d({k},{l})={value(distance[i][j])+value(distance[k][l])}, d({i},{k})+d({j},{l})={value(distance[i][k]) + value(distance[j][l])}, d({i},{l})+d({j},{k})={value(distance[i][l]) + value(distance[j][k])}") 


# Plot the tree
handle = io.StringIO(newick_str)
tree = Phylo.read(handle, "newick")
tree.ladderize()  # Flip branches so deeper clades are displayed at top
Phylo.draw(tree,branch_labels=lambda c: c.branch_length)
