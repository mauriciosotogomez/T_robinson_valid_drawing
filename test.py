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



    
 