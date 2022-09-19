import itertools
import pathlib
import csv
import os



def SaveParamsMCMC():
    param_value_coarse = [0.25, 0.5, 0.75]
    p = {'xi_AA': param_value_coarse,
         'xi_AR': param_value_coarse,
         'xi_RA': param_value_coarse,
         'xi_RR': param_value_coarse,
         'eta': param_value_coarse,
         'dname': ('df_regulon.csv',)}
    createParamsMCMC(p)
    p = {'xi_AA': param_value_coarse,
         'xi_AR': param_value_coarse,
         'xi_RA': param_value_coarse,
         'xi_RR': param_value_coarse,
         'eta': param_value_coarse,
         'dname': ('df_subti.csv',)}

    createParamsMCMC(p)
    p = {'xi_AA': [0.7, 0.75, 0.8],
         'xi_AR': [0.7, 0.75, 0.8],
         'xi_RA': [0.1, 0.15, 0.2, 0.25],
         'xi_RR': [0.1, 0.15, 0.2, 0.25],
         'eta': [0.4, 0.5, 0.6],
         'dname': ('df_regulon.csv',)}
    createParamsMCMC(p)
    #SUbti
    p = {'xi_AA': [0.9, 0.95, 0.99],
            'xi_AR': [0.9, 0.95, 0.99],
            'xi_RA': [0.05, 0.1, 0.15, 0.2],
            'xi_RR': [0.05, 0.1, 0.15, 0.2],
            'eta': [0.4, 0.5, 0.6],
            'dname': ('df_subti.csv',)}
    createParamsMCMC(p)



def createParamsMCMC(params):
    params_names = ['xi_A_A', 'xi_R_R', 'xi_A_R', 'xi_R_A', 'eta', 'dname']
    xi_AA_val = params['xi_AA']
    xi_AR_val = params['xi_AR']
    xi_RA_val = params['xi_RA']
    xi_RR_val = params['xi_RR']
    eta_val = params['eta']
    dname = params['dname']
    all_params = list(itertools.product(xi_AA_val,xi_RR_val,xi_AR_val,xi_RA_val, eta_val, dname))
    # Produce a list of dictionaries of parameters
    l = [{key: val for (key, val) in zip(params_names, ele)} for ele in all_params]
    path = pathlib.Path().absolute()
    for f in l:
        name = f"{path}/parameters/{f['xi_A_A']}_{f['xi_A_R']}_{f['xi_R_A']}_{f['xi_R_R']}_{f['eta']}_{f['dname'][:-4]}.txt"
        with open(name, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=params_names)
            writer.writeheader()
            writer.writerows([f])


if __name__=="__main__":
    root = "./"  # <-- set this
    path = os.path.join(root, 'parameters/')
    if not os.path.exists(path):
        # Create a parameters directory
        os.makedirs(path)
    path = os.path.join(root, 'output/')

    if not os.path.exists(path):
        # Create an output directory
        os.makedirs(path)
    SaveParamsMCMC()