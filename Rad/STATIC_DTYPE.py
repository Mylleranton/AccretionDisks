import numpy as np

def dtype_newtonian():
    dtype = np.dtype([
        ("ix", int),
        ("iy", int),
        ("iz", int),
        ("r", float),
        ("theta", float),
        ("phi", float),
        ("rho", float),
        ("temp", float),
        ("u_t", float),
        ("u_1", float),
        ("u_2", float),
        ("u_3", float),
        ("volume", float),
        ("bsq", float),
        ("b_1", float),
        ("b_2", float),
        ("b_3", float)
        ])
    return dtype

def dtype_newtonian_rad():
    dtype = np.dtype([
        ("ix", int),
        ("iy", int),
        ("iz", int),
        ("r", float),
        ("theta", float),
        ("phi", float),
        ("rho", float),
        ("temp", float),
        ("u_t", float),
        ("u_1", float),
        ("u_2", float),
        ("u_3", float),
        ("volume", float),
        ("rad_ehat", float),
        ("radflux_1", float),
        ("radflux_2", float),
        ("radflux_3", float),
        ("rad_tot_en", float),
        ("rad_compt_en", float),
        ("rad_bb_temp", float),
        ("rad_bb_temp2", float),
        ("bsq", float),
        ("b_1", float),
        ("b_2", float),
        ("b_3", float)
        ])
    return dtype
    
def dtype_rel_radiation():
    dtype = np.dtype([
        ("ix", int),
        ("iy", int),
        ("iz", int),
        ("r", float),
        ("theta", float),
        ("phi", float),
        ("rho", float),
        ("u_internal", float),
        ("u_t", float),
        ("u_1", float),
        ("u_2", float),
        ("u_3", float),
        ("volume", float),
        ("bsq", float),
        ("b_1", float),
        ("b_2", float),
        ("b_3", float),
        ("T_00", float),
        ("T_01", float),
        ("T_02", float),
        ("T_03", float),
        ("T_10", float),
        ("T_11", float),
        ("T_12", float),
        ("T_13", float),
        ("T_20", float),
        ("T_21", float),
        ("T_22", float),
        ("T_23", float),
        ("T_30", float),
        ("T_31", float),
        ("T_32", float),
        ("T_33", float),
        ("rad_ehat", float),
        ("R_00", float),
        ("R_01", float),
        ("R_02", float),
        ("R_03", float),
        ("R_10", float),
        ("R_11", float),
        ("R_12", float),
        ("R_13", float),
        ("R_20", float),
        ("R_21", float),
        ("R_22", float),
        ("R_23", float),
        ("R_30", float),
        ("R_31", float),
        ("R_32", float),
        ("R_33", float),
        ("G_time", float),
        ("G_radial", float),
        ("G_theta", float),
        ("G_phi", float)
        ])
    return dtype
    
def dtype_rel_hydro():
    dtype = np.dtype([
        ("ix", int),
        ("iy", int),
        ("iz", int),
        ("r", float),
        ("theta", float),
        ("phi", float),
        ("rho", float),
        ("u_internal", float),
        ("u_t", float),
        ("u_1", float),
        ("u_2", float),
        ("u_3", float),
        ("volume", float),
        ("bsq", float),
        ("b_1", float),
        ("b_2", float),
        ("b_3", float),
        ("T_00", float),
        ("T_01", float),
        ("T_02", float),
        ("T_03", float),
        ("T_10", float),
        ("T_11", float),
        ("T_12", float),
        ("T_13", float),
        ("T_20", float),
        ("T_21", float),
        ("T_22", float),
        ("T_23", float),
        ("T_30", float),
        ("T_31", float),
        ("T_32", float),
        ("T_33", float)
        ])
    return dtype
    