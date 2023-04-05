"""parameters are, in order: M, b, d0, dM, r_um, r_mu, p"""


TOY_DYING_PARAMS = [
    4, 
    2, 
    3, 
    1.5, 
    1, 
    1, 
    1, 
  ]

TOY_LIVING_PARAMS = [
    4, 
    2.2, 
    3, 
    1, 
    0.2, 
    0.2, 
    1, 
  ]

def TOY_INITIAL(M : int):
    n = [0] * (M + 1)
    n[0] = 100
    return n