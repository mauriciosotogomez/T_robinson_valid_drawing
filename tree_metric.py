# import utils functions
import utils

# Arguments
from argparse import ArgumentParser

import numpy as np # represent
#Solve LP
from pulp import * # solve LP

# Construct the tree (install scikit-bio)
from skbio import DistanceMatrix 
from skbio.tree import nj

# Plot the tree
## from Bio import Phylo
## from io import StringIO

#subprocess.call ("/pathto/MyrScript.r")

parser = ArgumentParser()
parser.add_argument("-e", "--example_number",
                   dest="n_example", default=0,
                   help="use this example from the list ")

parser.add_argument("-f", "--input_file",
                   dest="input_file", default=None,
                   help="use this .csv file as input")

parser.add_argument("-g", "--generate_random",
                    dest="generate_random_size", default=None,
                    help="generate a random matrix with this size")

parser.add_argument("-v",
                    action="store_true",
                    help="verbose")

parser.add_argument("-d",
                    action="store_true",
                    help="use distance matrix")
args = parser.parse_args()

# Get problem input matrix
# from the predefined list
input_matrix = np.array(utils.get_example(int(args.n_example)))
# from a file
if args.input_file is not None:
    input_matrix = np.genfromtxt(args.input_file, delimiter=',')
# random generated
if args.generate_random_size is not None:
    subprocess.call(["/usr/bin/Rscript", "create_Robinson.r",args.generate_random_size])
    input_matrix = np.genfromtxt("random_robinson_matrix.csv", delimiter=',')
    # it is a distance matrix

#######################
# Problem description #
#######################

n = len(input_matrix)
ROWS = COLS = range(n)
# Node labels
nodelabel = ["S"+str(i) for i in range(n)]
#all_labels=["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
#nodelabel = all_labels[:n]

input_type="distance"
distance_matrix=input_matrix 
if (not args.d) : # input similarity matrix
    input_type ="similarity"
    distance_matrix = input_matrix.max()-input_matrix

# Create the n-tuples  created, with the row and column index
four_points = [ # order i<j<k<l 
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

two_points = [ # i<j
    (i,j) 
    for j in range(1,n)
    for i in range(j)
]

###########################
# Create and solve the LP #
###########################
problem = LpProblem("Has_a_tree_diagram_Problem")

# Define decision variables
distance = LpVariable.dicts("distance", (ROWS, COLS), cat='Continuous', lowBound=0)

# We do not define an objective function since none is needed
#problem += lpSum([distance[i][j] for (i,j) in two_points])

# Constraints ensuring the 4-point condition
# WE ASSUME distance[i][l] + distance[j][k] >= distance[i][j] + distance[k][l] 
for (i,j,k,l) in four_points:
    problem += distance[i][k] + distance[j][l] - distance[i][l] - distance[j][k] == 0

# Constraints ensuring the 3-point condition
for (i,j,k) in three_points:
    # Triangular inequality
    problem += distance[i][j] + distance[j][k] - distance[i][k] >= 0
    #  Consistent with the input matrix
    if input_matrix[i][j] > input_matrix[i][k] :
        if (not args.d) :  # we are using a similarity matrix
            problem += distance[i][j] - distance[i][k] <= -1
        else : # we are using a distance matrix
            problem += distance[i][j] - distance[i][k] >= 1
        
# Symmetry 2-points condition
for (i,j) in two_points:
    problem += distance[i][j] - distance[j][i] == 0

# The problem data is written to an .lp file
problem.writeLP("tree_metric.lp")

# Solve the problem
problem.solve(PULP_CBC_CMD(msg=0))

# Compute output matrix
outmatrix=np.zeros((n,n))
for r in ROWS:
    for c in COLS:
       outmatrix[r][c] = value(distance[r][c])

# The solution is is created/overwritten to the outmatrix.txt file
matrixout = open("tree_metric_out.txt", "w")
for c in COLS:
    matrixout.write(f"{nodelabel[c]}, " )
    matrixout.write("\n")
for r in ROWS:
    matrixout.write(f"{nodelabel[r]}")
    for c in COLS:
        matrixout.write(f", {value(distance[r][c])}")
        matrixout.write("\n")
matrixout.close()

# The status of the solution is printed to the screen
print("Status:", LpStatus[problem.status])

# Check the solution
print("Correct solution:",end="")
is_correct=utils.check_solution(distance_matrix,outmatrix)
if is_correct :
    print("\033[1m\033[92m {}\033[00m" .format(is_correct)) # Green
else :
    print("\033[91m {}\033[00m" .format(is_correct)) # Red

############################
# PRINT INPUT/OUTPUT INFO  #
############################

if args.v :
    if (not args.d) : # input similarity matrix
        print(f"\n Input matrix ({input_type}) \n{input_matrix}")

    print(f"\n Input distance matrix \n{distance_matrix}")
     
    ## # Print the violations of triangular inequality
    print(f"\n Violations of triangular inequality")
    for (i,j,k) in three_points:
        # Triangular inequality
        delta  =distance_matrix[i][j] + distance_matrix[j][k] - distance_matrix[i][k] 
        if delta<0 :
            print(f"nodes: {(i,j,k)} delta: {delta}")

    # Print the output matrix
    np.set_printoptions(precision=1,floatmode='fixed', suppress=True)
    print(f"\n Output distance matrix \n{outmatrix}")

    # Print the difference between input and output
    print(f"\n Output difference matrix \n{outmatrix-distance_matrix}")

    ## # Print the tree text mode
    ## dm = DistanceMatrix(outmatrix, nodelabel)
    ## tree = nj(dm)
    ## newick_str = nj(dm, result_constructor=str)

    ## print("\n Newick notation of the tree")
    ## print(newick_str)

    ## print("\n Structure of the tree")
    ## print(tree.ascii_art())
    
    # Plot the tree
    utils.plot_text_caterpillar(outmatrix, nodelabel)
    utils.plot_caterpillar(outmatrix,nodelabel)
    
    ## handle = io.StringIO(newick_str)
    ## tree = Phylo.read(handle, "newick")
    ## tree.ladderize()  # Flip branches so deeper clades are displayed at top
    ## Phylo.draw(tree,branch_labels=lambda c: c.branch_length)

    
    
