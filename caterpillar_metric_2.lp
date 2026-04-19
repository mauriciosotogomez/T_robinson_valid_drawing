\* Has_a_caterpillar_diagram_Problem *\
Minimize
OBJ: __dummy
Subject To
C0.1: - distance_0 + distance_1  >= 1
C0.2: - distance_0 + distance_2  >= 1
C0.3: - distance_0 + distance_3  >= 1
C0.4:   distance_0 - 2 distance_2 + distance_4  >= 1
C1.0: - distance_0 + distance_1  >= 1
C1.2: - distance_1 + distance_2  >= 1
C1.3: - distance_1 + distance_3  >= 1
C1.4:   distance_1 - 2 distance_2 + distance_4  >= 1
C2.0: - distance_0 + 2 distance_1 - distance_2  >= 1
C2.1: - distance_1 + distance_2  >= 1
C2.3: - distance_2 + distance_3  >= 1
C2.4: - distance_2 + distance_4  >= 1
C3.0: - distance_0 + 2 distance_2 - distance_3  >= 1
C3.1: - distance_1 + 2 distance_2 - distance_3  >= 1
C3.2: - distance_2 + distance_3  >= 1
C3.4: - distance_3 + distance_4  >= 1
C4.0: - distance_0 + 2 distance_3 - distance_4  >= 1
C4.1: - distance_1 + 2 distance_3 - distance_4  >= 1
C4.2: - distance_2 + 2 distance_3 - distance_4  >= 1
C4.3: - distance_3 + distance_4  >= 1
Bounds
 __dummy = 0
End
