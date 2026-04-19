from typing import List
import utils

n=100
n_dag=0
is_equal=False

for n in range(10,200,10):
    n_dag=0
    for i in range (1000):
        M = utils.generar_matriz(n)
        M1 = [[1, 1, 1, 1, 3, 3], [2, 2, 2, 3, 3, 3], [2, 3, 3, 3, 3, 5], [2, 4, 4, 4, 4, 5], [4, 4, 4, 5, 5, 5], [4, 4, 6, 6, 6, 6]] 
        #graph=utils.construir_grafo(M)
        # print("============")
        # print(M)
        # print("============")
        #is_equal=(M==M1)
        output=utils.es_DAG(M)
        #print(output)
        if (output=="YES"):
            n_dag+=1
        #print(is_equal)
    print(n,":",n_dag)
