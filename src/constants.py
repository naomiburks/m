"""parameters are, in order: M, b, d0, dM, r_um, r_mu, p"""


BASIC_PARAMS = {
    'b_0': 1.0,
    'b_M': 1.0, 
    'd_0': 1.0,
    'd_M': 1.0, 
    'r_um': 1.0,
    'r_mu': 1.0, 
    'M': 10,
    'p': 1.0,
  }

BASIC_PARAMS_COLL = {
    'b_0': 1.0,
    'b_M': 1.0, 
    'd_0': 1.0,
    'd_M': 1.0, 
    'r_um': 1.0,
    'r_mu': 1.0, 
    'r_mu_u': 10.0,
    'r_um_m': 10.0,
    'M': 10,
    'p': 1.0,
  }


DYING_PARAMS = {
    'b_0': 2.0,
    'b_M': 2.0, 
    'd_0': 3.0,
    'd_M': 1.5, 
    'r_um': 1.0,
    'r_mu': 1.0, 
    'M': 4,
    'p': 1.0,
  }

LIVING_DEATHRATE_PARAMS = {
    'b_0': 1.8,
    'b_M': 1.8, 
    'd_0': 3,
    'd_M': 1, 
    'r_um': 0.2,
    'r_mu': 0.1, 
    'M': 10,
    'p': 1,
  }


LIVING_BIRTHRATE_PARAMS = {
    'b_0': 0.8,
    'b_M': 2.8, 
    'd_0': 2,
    'd_M': 2, 
    'r_um': 0.2,
    'r_mu': 0.1, 
    'M': 10,
    'p': 1,
  }

COLLABORATIVE_PARAMS = {
    'b_0': 1.0,
    'b_M': 1.0, 
    'd_0': 1.0,
    'd_M': 1.0, 
    'r_um': 1.0,
    'r_mu': 1.0,
    'r_um_m': 1.0,
    'r_mu_u': 1.0,
    'M': 10,
    'p': 1,
    
}

BARELY_LIVING_PARAMS = {
    'b_0': 1.8,
    'b_M': 1.8, 
    'd_0': 3,
    'd_M': 1, 
    'r_um': 0.036,
    'r_mu': 0.036, 
    'M': 100,
    'p': 1,
}

def UNMETHYLATED_INITIAL(M : int, cell_count = 100):
    n = [0] * (M + 1)
    n[0] = cell_count
    return n

def HALF_METHYLATED_INITIAL(M: int, cell_count = 100):
    n = [0] * (M + 1)
    n[M // 2] = cell_count
    return n

def ALL_EQUAL(M: int, cell_count = 100):
    n = [cell_count] * (M + 1)
    return n