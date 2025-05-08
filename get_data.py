from ase.io import read,write                                                                                                     
import numpy as np                                                                                                                
import matplotlib.pyplot as plt                                                                                                   
import os                                                                                                                                  

def get_data():
    with open("data_here.data", "r") as f:
        lines = [line.strip() for line in f]
    i=int(lines[0]) 
    j=int(lines[1]) 
    k=int(lines[2]) 
    r0=float(lines[3])
    file_name=lines[4]
    return i,j,k,r0,file_name

def get_plot(i,j,k,r0,file_name):
    mols = read("mdNVT.traj@:")                                                                                                       
    folder_number = os.path.basename(os.path.dirname(os.getcwd()))[-4:]
    if i==117 and j == 127:
        distSiSi = []                                                                                                                         
        for mol in mols:                                                                                                                  
            tmp_pos = mol.get_positions()                                                                                                 
            distSiSi.append(np.linalg.norm(tmp_pos[i]-tmp_pos[j]))                                                                        
        plt.figure(figsize=(6, 4),dpi=100)  
        plt.title(f"Si-Si{i}, {j}, {k}, {r0}, {file_name} {folder_number}")
        plt.scatter(range(len(mols)),distSiSi)                                                                                                
        plt.plot(range(len(mols)),[r0 for _ in mols],color="black")                                                                     
        plt.xlabel("step")
        plt.ylabel("r (A)") 
        plt.savefig("data_graphSiSi.png", bbox_inches='tight')#, dpi=300)  
        np.save("SiSi.npy",distSiSi)
        # plt.show() 
    else:
        distSiSi = []                                                                                                                         
        distij = []                                                                                                                         
        for mol in mols:                                                                                                                  
            tmp_pos = mol.get_positions()                                                                                                 
            distSiSi.append(np.linalg.norm(tmp_pos[117]-tmp_pos[127]))                                                                        
            distij.append(np.linalg.norm(tmp_pos[i]-tmp_pos[j]))                                                                        
        plt.figure(figsize=(6, 4),dpi=100)  
        plt.title(f"SiSi{i}, {j}, {k}, {r0}, {file_name} {folder_number}")
        plt.scatter(range(len(mols)),distSiSi)                                                                                                
        # plt.plot(range(len(mols)),[r0 for _ in mols],color="black")                                                                     
        plt.xlabel("step")
        plt.ylabel("r (A)") 
        plt.savefig("data_graphSiSi.png", bbox_inches='tight')#, dpi=300)  
        #plt.show()
        plt.close() 
        plt.figure(figsize=(6, 4),dpi=100)  
        plt.title(f"{i}, {j}, {k}, {r0}, {file_name} {folder_number}")
        plt.scatter(range(len(mols)),distij)                                                                                                
        plt.plot(range(len(mols)),[r0 for _ in mols],color="black")                                                                     
        plt.xlabel("step")
        plt.ylabel("r (A)") 
        plt.savefig("data_graphij.png", bbox_inches='tight')#, dpi=300)  
        #plt.show()
        plt.close() 
        np.save("SiSi.npy",distSiSi)
        np.save("ij.npy",distij)





if __name__=="__main__":
    i,j,k,r0,file_name = get_data()
    get_plot(i,j,k,r0,file_name)
