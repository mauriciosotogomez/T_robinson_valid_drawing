"""
Generate and validate a Robinson matrix from a random center matrix.

This script generates a center matrix of size n x n (strict or non-strict),
attempts to construct a Robinson matrix from it, and verifies correctness by
recomputing the center matrix and comparing.

Usage:
    python test.py <n> [s]

Arguments:
    n       Size of the matrix (integer).
    s       Optional flag. Pass 's' for strict Robinson matrix; omit for non-strict.

Examples:
    python test.py 6        # Generate a non-strict Robinson matrix of size 6
    python test.py 6 s      # Generate a strict Robinson matrix of size 6
"""

import utils
import sys


def main(n, strict="n"):
    
    n = int(n)
    is_strict = (strict=="s")

    if is_strict :
        C = utils.generate_strict_center_matrix(n)
    else :
        C = utils.generate_center_matrix(n)

    M = utils.generate_robinson_from_center_matrix(C) 

    if M is not False :
        Cm = utils.compute_max_closer(M)
        # Check solution
        print((Cm.astype(int)== C).all())
    else :
        print(False)

    print(M)
        # print(Cm)
    print(C)    

if __name__ == "__main__":
    main(*sys.argv[1:])



    
 