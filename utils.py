import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import io

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
    
def plot_caterpillar(distance_matrix, nodelabel, is_correct) :
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
        for i in range(n):
            for k in range(n):
                if (distance_matrix[i][j] > distance_matrix[i][k]) and  (outmatrix[i][j] < outmatrix[i][k]):  
                    is_correct=False 
    return(is_correct)
