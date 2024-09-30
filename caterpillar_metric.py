# import utils functions
import utils
import numpy as np
import time

# Arguments
from argparse import ArgumentParser

#Solve LP
from pulp import * # solve LP

# Construct the tree (install scikit-bio)
## from skbio import DistanceMatrix 
## from skbio.tree import nj

# Plot the tree
## from Bio import Phylo
## from io import StringIO


#######################
# Read Input  #
#######################

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

parser.add_argument("-p",
                    action="store_true",
                    help="use a path instead of centipede")

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

# Compute distance matrix, if necessary
input_type="distance"
distance_matrix=input_matrix 
if (not args.d) : # input similarity matrix
    input_type ="similarity"
    distance_matrix = input_matrix.max()-input_matrix

max_closer=utils.compute_max_closer(distance_matrix)
    
#######################
# Problem description #
#######################
n = len(input_matrix)
ROWS = COLS = range(n)
tol=0.11

# Node labels
nodelabel = ["S"+str(i) for i in range(n)]

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

# Pairs
two_points = [ 
    (i,j) 
    for j in range(n)
    for i in range(n)
]

###########################
# Create and solve the LP #
###########################
tic = time()
problem = LpProblem("Has_a_caterpillar_diagram_Problem", LpMinimize)

# Define decision variablesdistance 
distance = LpVariable.dicts("distance", ROWS, cat='Continuous', lowBound=0)
leg      = LpVariable.dicts("leg", ROWS, cat='Continuous', lowBound=0)

# We do not define an objective function since none is needed
problem += lpSum([leg[i] for i in range(n)])
#problem += distance[n-1]

# Constraints for each pair
for (i,j) in two_points:
    if (i<j-1) :
        problem += LpConstraint( distance[i] + distance[j] - 2*distance[max_closer[i,j]] - leg[i] + leg[j], sense=LpConstraintGE, rhs=1 , name='C'+str(i)+str(j))
        #problem += distance[i] + distance[j] - 2*distance[max_closer[i,j]] - leg[i] + leg[j] >= 1
    elif (i>j+1) :
        problem += LpConstraint( -distance[i] - distance[j] + 2*distance[max_closer[i,j]] - leg[i] + leg[j], sense=LpConstraintGE, rhs=1 , name='C'+str(i)+str(j))
        #problem += - distance[i] - distance[j] + 2*distance[max_closer[i,j]] - leg[i] + leg[j] >= 1
    elif (i==j-1) :
        problem += LpConstraint( distance[i] + distance[j] - 2*distance[max_closer[i,j]] - leg[i] + leg[j], sense=LpConstraintGE, rhs=0 , name='C'+str(i)+str(j))
        #problem += distance[i] + distance[j] - 2*distance[max_closer[i,j]] - leg[i] + leg[j] >= 0
        #problem += -distance[i] + distance[j] - leg[i] + leg[j] >= 0
    elif (i==j+1) :
        problem += LpConstraint( -distance[i] - distance[j] + 2*distance[max_closer[i,j]] - leg[i] + leg[j], sense=LpConstraintGE, rhs=0 , name='C'+str(i)+str(j))
        #problem += - distance[i] - distance[j] + 2*distance[max_closer[i,j]] - leg[i] + leg[j] >= 0
        #problem += distance[i] - distance[j] - leg[i] + leg[j] >= 0

# Forse to be in a path        
if args.p:
    for i in ROWS:
        problem += LpConstraint(leg[i], sense=LpConstraintEQ, rhs=0 , name='L'+str(i))
# The problem data is written to an .lp file
problem.writeLP("caterpillar_metric.lp")

# Solve the problemvalue(
problem.solve(PULP_CBC_CMD(msg=0))

toc = time()

# Compute output matrix
outmatrix=np.zeros((n,n))
for r in ROWS:
    for c in COLS:
        if (c>r) :
            outmatrix[r][c] = value(distance[c]) - value(distance[r]) + value(leg[c]) + value(leg[r])
        elif (c<r) :
            outmatrix[r][c] = value(distance[r]) - value(distance[c]) + value(leg[c]) + value(leg[r])

# The solution is is created/overwritten to the outmatrix.txt file
matrixout = open("caterpillar_metric_out.txt", "w")
for c in COLS:
    matrixout.write(f"{nodelabel[c]}, " )
    matrixout.write("\n")
for r in ROWS:
    matrixout.write(f"{nodelabel[r]}")
    for c in COLS:
        matrixout.write(f", {value(outmatrix[r][c])}")
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
print(f'Time: {toc - tic} seconds')

if args.v :     
    # # Print the violations of triangular inequality
    # print(f"\n Violations of triangular inequality - input distance matrix")
    # for (i,j,k) in three_points:
    #     # Triangular inequality
    #     delta  =distance_matrix[i][j] + distance_matrix[j][k] - distance_matrix[i][k] 
    #     if delta<0 :
    #         print(f"nodes: {(i,j,k)} delta: {delta}")
            
    # # Print the violations of 4-point condition
    # print(f"\n Violations of  4-point condition - input distance matrix")
    # for (i,j,k,l) in four_points:
    #     d1 = distance_matrix[i][l] + distance_matrix[j][k]
    #     d2 = distance_matrix[i][k] + distance_matrix[j][l]
    #     if (d1 - d2) > tol :
    #         print(f"nodes: {(i,j,k,l)} d1={d1} d2={d2}")

    if (not args.d) : # input similarity matrix
        print(f"\n Input matrix ({input_type}) \n{input_matrix}")

    print(f"\n Input distance matrix \n{distance_matrix}")

    # Print the max_closer matrix
    print(f"\n Max closer matrix\n{max_closer}")
    
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

    # Plot Primal/Dual Solution
    print("\nPrimal Variables")
    for v in problem.variables():
         print(v.name, ":" "\t", v.varValue)
        #varsdict[v.name] = v.varValue
    #print(varsdict)
       
    print("\nDual Variables")
    for name, c in list(problem.constraints.items()):
        print(name, ":" "\t", c.pi, "\t slack:", c.slack)

     
    # Plot the tree
    utils.plot_text_caterpillar(outmatrix, nodelabel)
    utils.plot_caterpillar(outmatrix,nodelabel)
   
    ## handle = io.StringIO(newick_str)
    ## tree = Phylo.read(handle, "newick")
    ## tree.ladderize()  # Flip branches so deeper clades are displayed at top
    ## Phylo.draw(tree,branch_labels=lambda c: c.branch_length)

    
    
