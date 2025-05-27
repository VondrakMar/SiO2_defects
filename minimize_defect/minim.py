from mace.calculators import mace_mp
import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase import build, units
from ase.md import Langevin
from ase.io.trajectory import Trajectory
from ase.md.velocitydistribution import (
    MaxwellBoltzmannDistribution,
    Stationary,
    ZeroRotation)
from ase.optimize import QuasiNewton, MDMin
from ase.filters import FrechetCellFilter
from ase.io import read,write
from ase.constraints import FixAtoms, Hookean
from ase.calculators.mixing import SumCalculator
from ase.visualize import view

def run_minimization(mol,
                    # calculator, # calculator should be attached before
                    stationary=True,
                    trj_name = "bfgs_ls",
                    fmax=0.05,steps=100000000000):
    if stationary:
        Stationary(mol)
    # mol.calc = calculator
    dyn = QuasiNewton(atoms=mol, trajectory=f'{trj_name}.traj')#, restart=f'{trj_name}.pckl')
    dyn.run(fmax=fmax,steps=steps)


def run_minimization_cell(mol,
                    # calculator, # calculator should be attached before
                    stationary=True,
                    trj_name = "bfgs_ls",
                    fmax=0.05,
                    steps=100000000000):
    if stationary:
        Stationary(mol)
    # mol.calc = calculator
    ecf = FrechetCellFilter(mol)
    qn = QuasiNewton(ecf)
    traj = Trajectory(f"{trj_name}.traj","w",mol)
    qn.attach(traj)
    qn.run(fmax=fmax,steps=steps) 
    # dyn = QuasiNewton(atoms=mol, trajectory=f'{trj_name}.traj')#, restart=f'{trj_name}.pckl') qn.run(fmax=fmax,steps=steps)


def get_distances(mol,structure):
    '''
        Indexed I am using here:
        Relevant Si are Si atoms of defects
        relevant O are oxygen bonded to these Si atoms
        secondery Os are oxygen behind 127 Si.
        From the paper table (I copied it here as is in the paper, but they have switched SiA-O and SiB-O distances (I checked for pucker they sent me, but Ebp also is not correct)
        ___________________________________________________________________________________
                                              Data for alpha-Quartz:
        ___________________________________________________________________________________
                   240ind          |  V_O   |  E_di    | Epuck |    Ebp     |    V_Om     |  Pristin
        SiA-SiB| 127-117           | 2.41   |  3.00    | 4.45  |    5.38    |   2.47      |    3.09 
        SiA-OB | 127-156           | 4.10   |  3.71    | 1.83  |    1.84    |   4.07      |    3.66 
        SiA-O  | 127-(124,130,131) | 1.63   |  1.58    | 1.62  |  1.61–1.64 |  1.69–1.73  |    1.61 
        SiB-O  | 117-(96,97,120)   | 1.62   |  1.57    | 1.58  |  1.57–1.59 |  1.66–1.68  |    1.61 
        ___________________________________________________________________________________
    '''
    tmp_pos = mol.get_positions()
    if structure == "big":
        '''
        Not all line here are relevant, but you can use them to pick a defect site
        '''
        relevant_Si = np.array([117,127])
        relevant_O = np.array([97, 98, 120, 124, 130, 131, 156])
        secondery_O = np.array([93,94,95,96,177])
        defect_site = np.union1d(relevant_Si,relevant_O)
        defect_site = np.union1d(defect_site,secondery_O)
        SiA_neighOs = np.array([124,130,131])
        SiB_neighOs = np.array([97,98,120])
        SiA_Ob = np.array([156])
        SiAO_dists = []
        SiBO_dists = []
        SiASiB_dist = np.linalg.norm(tmp_pos[117]-tmp_pos[127])
        SiAOb_dist = np.linalg.norm(tmp_pos[156]-tmp_pos[127])
        for Oind in SiA_neighOs:
            SiAO_dists.append(np.linalg.norm(tmp_pos[127]-tmp_pos[Oind]))
        for Oind in SiB_neighOs:
            SiBO_dists.append(np.linalg.norm(tmp_pos[117]-tmp_pos[Oind]))
        # defect_site_mol = mol[defect_site]
        # view(defect_site_mol)
        return {"SiAOs":SiAO_dists,"SiBOs":SiBO_dists,"SiOb":SiAOb_dist,"SiSi":SiASiB_dist}

def compare_defects(mol,
                    structure="big", # big for 240 atom structure, "small" for 144 atoms
                    charge=None,
                    use_paper=False):
    mol_dists = get_distances(mol,structure)
    if structure == "big":
        Ebp = read("Ebp.xyz",format="extxyz")
        Ebp_dists = get_distances(Ebp,structure)
        Edi = read("Edi.xyz",format="extxyz")
        Edi_dists = get_distances(Edi,structure)
        Epuc = read("Epuc.xyz",format="extxyz")
        Epuc_dists = get_distances(Epuc,structure)
        V0 = read("V0.xyz",format="extxyz")
        V0_dists = get_distances(V0,structure)
        Vm = read("Vm.xyz",format="extxyz")
        Vm_dists = get_distances(Vm,structure)
            
    elif structure == "small":
        print("Has to be implemented")
        return 0
    else:
        print("only big (240 atoms) and small (144 atoms) is doable")
        return -1
    dist_Si = {
        "Ebp":np.linalg.norm(mol_dists["SiSi"]-Ebp_dists["SiSi"]),
        "Edi":np.linalg.norm(mol_dists["SiSi"]-Edi_dists["SiSi"]),
        "Epuc":np.linalg.norm(mol_dists["SiSi"]-Epuc_dists["SiSi"]),
        "V0":np.linalg.norm(mol_dists["SiSi"]-V0_dists["SiSi"]),
        "Vm":np.linalg.norm(mol_dists["SiSi"]-Vm_dists["SiSi"])
    }
    dist_SiOb = {
        "Ebp":np.linalg.norm(mol_dists["SiOb"]-Ebp_dists["SiOb"]),
        "Edi":np.linalg.norm(mol_dists["SiOb"]-Edi_dists["SiOb"]),
        "Epuc":np.linalg.norm(mol_dists["SiOb"]-Epuc_dists["SiOb"]),
        "V0":np.linalg.norm(mol_dists["SiOb"]-V0_dists["SiOb"]),
        "Vm":np.linalg.norm(mol_dists["SiOb"]-Vm_dists["SiOb"])
    }
    min_SiSidiff = min(dist_Si, key=dist_Si.get)
    min_SiObdiff = min(dist_SiOb, key=dist_SiOb.get)
    if min_SiSidiff in ["V0","Vm","Edi"] and min_SiObdiff in ["V0","Vm","Edi"]:
        print(f"Current defect is probably a dimer with Si-Si distance {mol_dists['SiSi']:.2f} and Si-Ob distance {mol_dists['SiOb']:.2f}")
    elif min_SiSidiff in ["Ebp"] and min_SiSidiff in ["Ebp","Epuc"]:
        print(f"Current defect is probably a back-project pucker Ebp with Si-Si distance {mol_dists['SiSi']:.2f} and Si-Ob distance {mol_dists['SiOb']:.2f}")
    elif min_SiSidiff in ["Epuc"] and min_SiSidiff in ["Ebp","Epuc"]:
        print(f"Current defect is probably a forward project pucker Epuc with Si-Si distance {mol_dists['SiSi']:.2f} and Si-Ob distance {mol_dists['SiOb']:.2f}")
    return min_SiObdiff
    
if __name__ == "__main__":
    mol = read("struc0020.xyz",format="extxyz")
    # macemp = mace_mp(model="MACE-matpes-pbe-omat-ft.model",device="cuda",dispersion=False,enable_cueq=True,dtype="float64")
    # run_minimization(mol,fmax=0.001,steps=1000,trj_name=f"struc_bfgs")
    defect = compare_defects(mol,"big")
    print("Probable defect",defect)
