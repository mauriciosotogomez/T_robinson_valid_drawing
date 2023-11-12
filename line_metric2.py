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

all_input_matrix = [
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
 [1,1,7,7,9]]
]

input_matrix = np.array(all_input_matrix[n_example])

n = len(input_matrix)
ROWS = COLS = range(n)
nodelabel = all_labels[:n]
nodelabel = ["S"+str(i) for i in range(n)]


# The boxes list is created, with the row and column index of each square in each box
# four_points = [
#     (i,j,k,l) 
#     for l in range(3,n)
#     for k in range(2,l)
#     for j in range(1,k)
#     for i in range(j)
# ]

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

# Create the problem
problem = LpProblem("Has a line diagram Problem",LpMinimize)

# Define decision variables
distance = LpVariable.dicts("distance", (ROWS, COLS), cat='Continuous', lowBound=0)
max = LpVariable.dicts("max", (three_points), cat='Continuous',  lowBound=0)

# We define an objective function since none is needed
problem += lpSum([distance[i][j] + distance[j][k] + distance[i][k]- 2*max[(i,j,k)] for (i,j,k) in three_points])

# Constraints ensuring the 3-point condition
# consistent with the matrix and line metric
#  d(i,j) + d(j,k) + d(i,k) <= 2*max{ d(i,k), d(k,j), d(i,j) } 
for (i,j,k) in three_points:
    if input_matrix[i][j] > input_matrix[i][k] :
        problem += distance[i][j] <= distance[i][k] - 1
    problem += distance[i][j] + distance[j][k] + distance[i][k] == 2*max[(i,j,k)] 
    problem += distance[i][j] <= max[(i,j,k)] 
    problem += distance[j][k] <= max[(i,j,k)] 
    problem += distance[i][k] <= max[(i,j,k)]

# Symmetry 2-points condition
for (i,j) in two_points:
    problem += distance[i][j] == distance[j][i]

# Positive 1-points condition
for i in COLS:
    problem += distance[i][i] == 0
    
# The problem data is written to an .lp file
problem.writeLP("line_metric.lp")

# Solve the problem
problem.solve()

#PRINT INPUT/OUTPUT INFO

# The status and value of the solution is printed to the screen
print("Status:", LpStatus[problem.status])
print(f'Value: {value(problem.objective)}')

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

# The solution is written to the sudokuout.txt file
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

for (i,j,k) in three_points :
     if value(distance[i][j]) != value(max[(i,j,k)]) and value(distance[j][k]) != value(max[(i,j,k)]) and  value(distance[i][k]) != value(max[(i,j,k)]) :
         print(f"[{(i,j,k)}] max={value(max[(i,j,k)])}: d(i,j)={value(distance[i][j])}, d(j,k)={value(distance[j][k])}, d(i,k)={value(distance[i][k])}") 

# Plot the tree
handle = io.StringIO(newick_str)
tree = Phylo.read(handle, "newick")
tree.ladderize()  # Flip branches so deeper clades are displayed at top
Phylo.draw(tree,branch_labels=lambda c: c.branch_length)
