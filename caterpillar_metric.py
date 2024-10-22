# import utils functions
import utils
import numpy as np
import time
import pandas as pd
import re

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

#### AUX PROBLEM #################################################################################################################################################################################
# Time aux
tic_aux = time()
#Aux problem                                                                                                                                                                                   #
aux_problem = LpProblem("A_has_a_positive_linear_combination", LpMinimize)                                                                                                                     #
aux_distance = LpVariable.dicts("aux_distance", ROWS, cat='Continuous')                                                                                                                        #
aux_leg      = LpVariable.dicts("aux_leg", ROWS, cat='Continuous')                                                                                                                             #
aux_problem += 0 # feasibility                                                                                                                                                                 #
# Constraints for each pair                                                                                                                                                                    #
for (i,j) in two_points:                                                                                                                                                                       #
    if (i<j-1) :                                                                                                                                                                               #
        aux_problem += LpConstraint( aux_distance[i] + aux_distance[j] - 2*aux_distance[max_closer[i,j]] - aux_leg[i] + aux_leg[j], sense=LpConstraintGE, rhs=1 , name='A'+str(i)+'_'+str(j))  #
    elif (i>j+1) :                                                                                                                                                                             #
        aux_problem += LpConstraint( -aux_distance[i] - aux_distance[j] + 2*aux_distance[max_closer[i,j]] - aux_leg[i] + aux_leg[j], sense=LpConstraintGE, rhs=1 , name='A'+str(i)+'_'+str(j)) #
    elif (i==j-1) :                                                                                                                                                                            #
        aux_problem += LpConstraint( aux_distance[i] + aux_distance[j] - 2*aux_distance[max_closer[i,j]] - aux_leg[i] + aux_leg[j], sense=LpConstraintGE, rhs=1 , name='A'+str(i)+'_'+str(j))  #
    elif (i==j+1) :                                                                                                                                                                            #
        aux_problem += LpConstraint( -aux_distance[i] - aux_distance[j] + 2*aux_distance[max_closer[i,j]] - aux_leg[i] + aux_leg[j], sense=LpConstraintGE, rhs=1 , name='A'+str(i)+'_'+str(j)) #
# Force to be in a path        
if args.p:
    for i in ROWS:
        aux_problem += LpConstraint(aux_leg[i], sense=LpConstraintEQ, rhs=0 , name='AL'+str(i))

# Solve the aux_problem                                                                                                                                                                  #
aux_problem.solve(PULP_CBC_CMD(msg=0))
# time it
toc_aux = time()
# The status of the solution is printed to the screen                                                                                                                                          #
print("AUX Status:", LpStatus[aux_problem.status])#
# Print time
print(f'Time AUX: {toc_aux - tic_aux} seconds')
# Plot Primal/Dual Solution                                                                                                                                                                    #
print("\nAUX Primal Variables")                                                                                                                                                                #
for v in aux_problem.variables():                                                                                                                                                              #
    print(v.name, ":" "\t", v.varValue)                                                                                                                                                        #

# Create Matrix A
# create DataFrame
list_constraints =  list(aux_problem.constraints)
list_variables = list()
for var in aux_problem.variables():
    list_variables.append(var.name)
# remove dummy variable
list_variables=list_variables[1:]

# create empty matrix A
A = pd.DataFrame(np.zeros(shape=(len(list_constraints),len(list_variables))), columns=list_variables, index=list_constraints)
# problem dictionary
dict = aux_problem.to_dict()
# fill matrix A
for constraint in dict['constraints']:
    for coefficient in constraint['coefficients']:
        A.loc[constraint['name'], coefficient['name']]=coefficient['value']
A.to_csv('A.csv')

# Create matrix B
l2 = int(len(list_variables)/2) # two parts
columns_B = ['d+l_'+str(i) for i in range(l2)]+['d-l_'+str(i) for i in range(l2)]
B = pd.DataFrame(np.zeros(shape=(len(list_constraints),len(list_variables))), columns=columns_B, index=list_constraints)

for i in range(l2):
    B.iloc[:,i] = (A.iloc[:,i]+A.iloc[:,i+l2])/2
    B.iloc[:,i+l2] = (A.iloc[:,i]-A.iloc[:,i+l2])/2
B.astype(int).to_csv('B.csv')

# Order rows of B 
# add type constraint

type_constraint = []
type_constraint_order = []
groups_size =[0,0,0]

for (i,j) in two_points:
    if i != j :
        type = 0
        order = i
        if (i<j) :                                                                                                                                                                               #
            if max_closer[i,j] == i:
                type = 0
                order = i
            else :
                type = 1
                order = j
        elif (i>j) :                                                                                                                                                                         
            if max_closer[i,j] == i:
                type = 0
                order = i
            else :
                type = 2
                order = max_closer[i,j]
        type_constraint.append(int(type))
        type_constraint_order.append(int(order))
        groups_size[type] = groups_size[type]+1

# Add type and order to sort the rows
B['type_constraint']=type_constraint
B['type_constraint_order']=type_constraint_order
B.sort_values(by=['type_constraint', 'type_constraint_order'], ascending=[True, True], inplace=True)
print(B)

# Pop the new order rows
type_constraint=[int(x) for x in list(B['type_constraint'])]
type_constraint_order=[int(x) for x in list(B['type_constraint_order'])]
B.drop(columns=['type_constraint', 'type_constraint_order'],inplace=True)

# REPAIR COEFFICIENTS
coeff = np.array([i for i in range(l2)]+[i for i in range(l2)])
B.to_numpy()
result = np.dot(B,coeff)

print("================"+" 0 "+"=====================")
print(coeff)
print(list(map(int,result)))
#print(type_constraint)
print()

for i in range(len(result)) :
     if result[i] <= 0 :
        if type_constraint[i] == 1 :
            k = type_constraint_order[i]
            while k < l2 :
                coeff[k] = coeff[k]+groups_size[2]-result[i]+1
                k=k+1
        if type_constraint[i] == 2 :
            k = l2+type_constraint_order[i]
            while k < len(coeff) :
                coeff[k] = coeff[k]-result[i]+1
                k=k+1
        result = np.dot(B,coeff)

        print("================ "+"row "+list_constraints[i]+" type "+str(type_constraint[i])+" =====================")
        print(coeff)
        print(list(map(int,result)))
        print()

# Print solution from Aux
legs = (coeff[:l2]-coeff[l2:])/2
distance2 = (coeff[:l2]+coeff[l2:])/2  
distance = coeff[:l2]-legs  
print(f'distances : {distance}')
#print(f'distances2: {distance2}')
print(f'legs      : {legs}')
print()


# # Compute the SVD of the matrix
# U, sigma, VT = np.linalg.svd(B)

# np.savetxt("U.csv", U, delimiter=",")
# np.savetxt("Sigma.csv", sigma, delimiter=",")
# np.savetxt("VT.csv", VT, delimiter=",")

 ##################################################################################################################################################################################################

###########################
# Create and solve the LP #
###########################
tic = time()

problem = LpProblem("Has_a_caterpillar_diagram_Problem", LpMinimize)
#problem = LpProblem("Has_a_caterpillar_diagram_Problem", LpMaximize)

# Define decision variablesdistance 
distance = LpVariable.dicts("distance", ROWS, cat='Continuous', lowBound=0)
leg      = LpVariable.dicts("leg", ROWS, cat='Continuous', lowBound=0)

# We do not define an objective function since none is needed
problem += lpSum([0*leg[i] for i in range(n)])
#problem += distance[n-1]

# Constraints for each pair
for (i,j) in two_points:
    if (i<j-1) :
        problem += LpConstraint( distance[i] + distance[j] - 2*distance[max_closer[i,j]] - leg[i] + leg[j], sense=LpConstraintGE, rhs=1 , name='C'+str(i)+'.'+str(j))
        #problem += distance[i] + distance[j] - 2*distance[max_closer[i,j]] - leg[i] + leg[j] >= 1
    elif (i>j+1) :
        problem += LpConstraint( -distance[i] - distance[j] + 2*distance[max_closer[i,j]] - leg[i] + leg[j], sense=LpConstraintGE, rhs=1 , name='C'+str(i)+'.'+str(j))
        #problem += - distance[i] - distance[j] + 2*distance[max_closer[i,j]] - leg[i] + leg[j] >= 1
    elif (i==j-1) :
        problem += LpConstraint( distance[i] + distance[j] - 2*distance[max_closer[i,j]] - leg[i] + leg[j], sense=LpConstraintGE, rhs=1 , name='C'+str(i)+'.'+str(j))
        #problem += distance[i] + distance[j] - 2*distance[max_closer[i,j]] - leg[i] + leg[j] >= 0
        #problem += -distance[i] + distance[j] - leg[i] + leg[j] >= 0
    elif (i==j+1) :
        problem += LpConstraint( -distance[i] - distance[j] + 2*distance[max_closer[i,j]] - leg[i] + leg[j], sense=LpConstraintGE, rhs=1 , name='C'+str(i)+'.'+str(j))
        #problem += - distance[i] - distance[j] + 2*distance[max_closer[i,j]] - leg[i] + leg[j] >= 0
        #problem += distance[i] - distance[j] - leg[i] + leg[j] >= 0

# Force to be in a path        
if args.p:
    for i in ROWS:
        problem += LpConstraint(leg[i], sense=LpConstraintEQ, rhs=0 , name='L'+str(i))

        # The problem data is written to an .lp file
problem.writeLP("caterpillar_metric.lp")

# Solve the LP
problem.solve(PULP_CBC_CMD(msg=0))

# Stop clock
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

# Check solution
is_correct=utils.check_solution(distance_matrix,outmatrix)

############################
# PRINT REPORT INFO  #
############################

# The status of the solution is printed to the screen
print("Status:", LpStatus[problem.status])

# Check the solution
print("Correct solution:",end="")
if is_correct :
    print("\033[1m\033[92m {}\033[00m" .format(is_correct)) # Green
else :
    print("\033[91m {}\033[00m" .format(is_correct)) # Red

# Print time
print(f'Time: {toc - tic} seconds')

# Print info (if -v/verbose )
if args.v :     
    if (not args.d) : # input similarity matrix
        print(f"\n Input matrix ({input_type}) \n{input_matrix}")

    print(f"\n Input distance matrix \n{distance_matrix}")

    # Print the max_closer matrix
    print(f"\n Max closer matrix\n{max_closer}")
    
    # Print the output matrix
    np.set_printoptions(precision=1,floatmode='fixed', suppress=True)
    print(f"\n Output distance matrix \n{outmatrix}")

    # Print the difference between input and output
    #print(f"\n Output difference matrix \n{outmatrix-distance_matrix}")

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
    # utils.plot_caterpillar(outmatrix,nodelabel,is_correct)
   
    ## handle = io.StringIO(newick_str)
    ## tree = Phylo.read(handle, "newick")
    ## tree.ladderize()  # Flip branches so deeper clades are displayed at top
    ## Phylo.draw(tree,branch_labels=lambda c: c.branch_length)

    

    ############################
    # PRINT REPORT INFO  #
    ############################

    # The status of the solution is printed to the screen
    print("Status:", LpStatus[problem.status])

    # Check the solution
    print("Correct solution:",end="")
    if is_correct :
        print("\033[1m\033[92m {}\033[00m" .format(is_correct)) # Green
    else :
        print("\033[91m {}\033[00m" .format(is_correct)) # Red

        # Print time
        print(f'Time: {toc - tic} seconds')

