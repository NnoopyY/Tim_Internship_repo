import numpy as np 
import matplotlib.pyplot as plt 
import re
import itertools
import matplotlib.lines as mlines
import matplotlib.colors as mcolors


#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


#**********************************************************************************
#**************************** Open the files **************************************
#**********************************************************************************
f1 = open('PROCAR', 'r')
PROCAR = f1.readlines()
f2 = open('DOSCAR', 'r')
DOSCAR = f2.readlines()
f3 = open('POSCAR', 'r')
POSCAR = f3.readlines()
#-------------- Extract info on the system ---------------------
#Extract the numlber of kpt band and ions
nkpoint=int(PROCAR[1].strip().split()[3])
nband=int(PROCAR[1].strip().split()[7])
nions=int(PROCAR[1].strip().split()[11])
E_fermi = float(DOSCAR[5].strip().split()[3])
blocksize_DOS = int(DOSCAR[5].strip().split()[2])
bandsize = 5+nions
blocksize_PRO = nband*bandsize
linesize = len(DOSCAR[blocksize_DOS + 8].strip().split())

#Extract ions list and ions numbers
ions_list = []
ions_number = []
for i in range(len(POSCAR[5].strip().split())):
    ions_list.append(str(POSCAR[5].strip().split()[i]))
    ions_number.append(int(POSCAR[6].strip().split()[i]))
projections_list = PROCAR[7].strip().split()[1:-1]
#Create an array for the index of the sum used later
ions_numbersum = np.zeros_like(ions_number)
for i in range(len(ions_number)-1):
    ions_numbersum[i+1] = ions_numbersum[i]+ions_number[i]

#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

#********************************************************************************
#******************* Extracting kcoord + HS points ******************************
#********************************************************************************
Energy_dos = np.zeros(blocksize_DOS)
Dosp = np.zeros(blocksize_DOS)
Dosm = np.zeros(blocksize_DOS)

# Extract the energy and DOS values from the DOSCAR file
for i in range(blocksize_DOS):
    index=6+i
    Energy_dos[i] = float(DOSCAR[index].strip().split()[0])
    Dosp[i]= float(DOSCAR[index].strip().split()[1])
    Dosm[i]= float(DOSCAR[index].strip().split()[2])

# Extract the energy par kpoints du PROCAR 
kbandenergyup=np.zeros([nkpoint, nband])
kbandenergydown = np.zeros([nkpoint, nband])
for i in range(nkpoint):
    for j in range(nband):
        kbandenergyup[i,j]=float(PROCAR[5+i*(blocksize_PRO+3)+j*bandsize].strip().split()[4])
for i in range(nkpoint):
    for j in range(nband):
        kbandenergydown[i,j]=float(PROCAR[6+nkpoint*(blocksize_PRO+3)+i*(blocksize_PRO+3)+j*bandsize].strip().split()[4])




#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


#********************************************************************************
#******************* Extracting kcoord + HS points ******************************
#********************************************************************************

#******************Extract the k coordinate for each kpoints*********************
kcoord = []
for i in range(nkpoint):
    if (nions == 1):
        line = PROCAR[i*((nband*5)+3)+3].strip()
    else :
        line = PROCAR[i*(blocksize_PRO+3)+3].strip()
    #Insert a space before a '-' that follows a digit (but not at the beginning)
    line = re.sub(r'(?<=\d)-', ' -', line)
    parts = line.split()
    # Now extract the desired elements. Adjust the indices if needed.
    kpt = [float(val) for val in parts[3:6]]
    kcoord.append(kpt)
kcoord = np.array(kcoord)

# Compute the cumulative k-point distance along the path
kDis = np.zeros(len(kcoord))
kDis[0] = 0
for i in range(len(kcoord)-1):
    if (np.sqrt((kcoord[i+1,0] - kcoord[i,0])**2 + (kcoord[i+1,1] - kcoord[i,1])**2 + (kcoord[i+1,2] - kcoord[i,2])**2)>0.4):
        kDis[i+1] = kDis[i] + np.sqrt((kcoord[i+1,0] - kcoord[i+1,0])**2 + (kcoord[i+1,1] - kcoord[i+1,1])**2 + (kcoord[i+1,2] - kcoord[i+1,2])**2)
    else :
        kDis[i+1] = kDis[i] + np.sqrt((kcoord[i+1,0] - kcoord[i,0])**2 + (kcoord[i+1,1] - kcoord[i,1])**2 + (kcoord[i+1,2] - kcoord[i,2])**2)

# Define the special high-symmetry points as numpy arrays.
special_points = {
    'Γ': np.array([0.0, 0.0, 0.0]),
    'X': np.array([0.5, 0.0, 0.0]),
    # 'M': np.array([0.5, 0.5, 0.0]),
    'R': np.array([0.5, 0.5, 0.5]),
    #'K': np.array([0.375, 0.375, 0.75]),
    #'L': np.array([0.5, 0.5, 0.5])
    'S': np.array([0.5, 0.5, 0.0]),
    'Y': np.array([0.0, 0.5, 0.0]),
    'Z': np.array([0.0, 0.0, 0.5]),
    'U': np.array([0.5, 0.0, 0.5]),
    'T': np.array([0.0, 0.5, 0.5]),
}

# Enter here the desired order (even with repeated points)
path_order = ['Γ', 'X', 'S', 'Y', 'Γ', 'Z', 'U', 'R', 'T', 'Z']
tol = 1e-3  # Tolerance for matching coordinates
tick_positions = [] 
tick_labels = []
current_search_index = 0

for sym in path_order:
    found_index = None
    # Search for the next occurrence of the special point (given by the symbol) starting at current_search_index
    for i in range(current_search_index, len(kcoord)):
        if np.linalg.norm(kcoord[i] - special_points[sym]) < tol:
            found_index = i
            break
    if found_index is not None:
        pos = kDis[found_index]
        # Check if this position is essentially the same as the previous one (within tolerance)
        if tick_positions and abs(tick_positions[-1] - pos) < tol:
            # Combine the labels with a vertical bar
            tick_labels[-1] = tick_labels[-1] + '|' + sym
        else:
            tick_positions.append(pos)
            tick_labels.append(sym)
        current_search_index = found_index + 1  # Continue search from here
    else:
        print(f"Warning: Could not find high-symmetry point {sym} starting from index {current_search_index}.")



#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


#*******************************************************************************************
#****************************** PROJECTIONS for Mn only *************************************
#*******************************************************************************************
# Create the arrays to stock the values
ions_list2 = ['Mndown', 'Mnup']
ions_index = [5, 7]

Dosions_projup = {ions: {proj: np.zeros(blocksize_DOS) for proj in projections_list} for ions in ions_list2}
Dosions_projdown = {ions: {proj: np.zeros(blocksize_DOS) for proj in projections_list} for ions in ions_list2}

#Loop for the filling
index=-1
for ions in ions_list2:
    index += 1
    indexproj = 0
    for proj in projections_list:
        indexproj += 1
        for i in range(blocksize_DOS):
            index_block = 6 + (1+ions_index[index]) * (blocksize_DOS+1) + i 
            Dosions_projdown[ions][proj][i]=float(DOSCAR[index_block].strip().split()[indexproj*2-1])
            Dosions_projup[ions][proj][i]=float(DOSCAR[index_block].strip().split()[indexproj*2])

#Extract the PROCAR values
Bandions_projdown = {ions: {proj: np.zeros((nkpoint, nband)) for proj in projections_list} for ions in ions_list2}
Bandions_projup = {ions: {proj: np.zeros((nkpoint, nband)) for proj in projections_list} for ions in ions_list2}

index = -1
for ions in ions_list2:
    index += 1
    proj_index = 0
    for proj in projections_list:  
        proj_index += 1  
        for i in range(nkpoint):
            for j in range(nband):
                Bandions_projdown[ions][proj][i][j] = float(PROCAR[8+i*(blocksize_PRO+3)+j*bandsize+ions_index[index]].strip().split()[proj_index])
                Bandions_projup[ions][proj][i][j] = float(PROCAR[9+nkpoint*(blocksize_PRO+3)+i*(blocksize_PRO+3)+j*bandsize+ions_index[index]].strip().split()[proj_index])


#inverting the contribution for projection dxy and dx2-y2 for the structure Pnma
for ion in ions_list2:
    projBandupdxy = Bandions_projup[ion]['dxy']
    projBandupdx2y2 = Bandions_projup[ion]['x2-y2']
    projBanddowndxy = Bandions_projdown[ion]['dxy']  
    projBanddowndx2y2 = Bandions_projdown[ion]['x2-y2'] 
    Bandions_projup[ion]['dxy'] = projBandupdx2y2
    Bandions_projup[ion]['x2-y2'] = projBandupdxy
    Bandions_projdown[ion]['dxy'] = projBanddowndx2y2
    Bandions_projdown[ion]['x2-y2'] = projBanddowndxy

#inverting the contribution for projection dxy and dx2-y2 for the structure Pnma
for ion in ions_list2:
    projdosupdxy = Dosions_projup[ion]['dxy']
    projdosupdx2y2 = Dosions_projup[ion]['x2-y2']
    projdosdowndxy = Dosions_projdown[ion]['dxy']  
    projdosdowndx2y2 = Dosions_projdown[ion]['x2-y2'] 
    Dosions_projup[ion]['dxy'] = projdosupdx2y2
    Dosions_projup[ion]['x2-y2'] = projdosupdxy
    Dosions_projdown[ion]['dxy'] = projdosdowndx2y2
    Dosions_projdown[ion]['x2-y2'] = projdosdowndxy

for ion in ions_list2:
    Bandions_projup[ion]['dxz+dyz'] = Bandions_projup[ion]['dxz'] + Bandions_projup[ion]['dyz'] 
    Bandions_projdown[ion]['dxz+dyz'] = Bandions_projdown[ion]['dxz'] + Bandions_projdown[ion]['dyz'] 

    Dosions_projup[ion]['dxz+dyz'] = Dosions_projup[ion]['dxz'] + Dosions_projup[ion]['dyz'] 
    Dosions_projdown[ion]['dxz+dyz'] = Dosions_projdown[ion]['dxz'] + Dosions_projdown[ion]['dyz'] 

projections_list = ['s', 'py', 'px', 'pz', 'dxy', 'dxz+dyz', 'dz2', 'x2-y2']

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------- Plotting projections contributions for Mn ions  ---------------------------------
#---------------------------------------------------------------------------------------------------------------------

color_list = ['gray', 'cyan', 'green', 'indianred', 'blue', 'green', 'red', 'orange', 'indianred']
colors = {}

for i, proj in enumerate(projections_list):
    colors[proj] = color_list[i % len(color_list)]
projections_list2 = projections_list
projections_list = projections_list[4:9]
# Marker size scaling factor for the scatter plots
marker_scale = 100


ok = ['positive', 'negative']
sign = {}
for i, ions in enumerate(ions_list2):
    sign[ions] = ok[i % len(ok)]

# For each ion, create a figure with DOS and band structure subplots.
for ion in ions_list2:
    # Create a new figure for each ion
    fig_proj, (ax_band_projdown, ax_dos, ax_band_projup) = plt.subplots(1,3, sharey=True, 
                                            gridspec_kw={'width_ratios': [2, 1, 2]},
                                            figsize=(12, 8))
    # fig_proj.suptitle(f"Mn {sign[ion]} magnetic moment ", fontsize=16)
    #----------------------- Plot the DOS contributions ---------------------------
    #   Fill the area under the DOS curves with grey (assuming these variables are global,
    # i.e. they are the same for each ion; if they are ion specific, you should replace them accordingly)
    ax_dos.fill_betweenx(Energy_dos - E_fermi, 0, Dosp, facecolor='grey', alpha=0.3)
    ax_dos.fill_betweenx(Energy_dos - E_fermi, -Dosm, 0, facecolor='grey', alpha=0.3)
    # Plot the DOS lines on top
    ax_dos.plot(Dosp, Energy_dos - E_fermi, color='black', lw=1, label='DOS (spin maj)')
    ax_dos.plot(-Dosm, Energy_dos - E_fermi, color='black', ls="--", lw=1, label='DOS (spin min)')
    ax_dos.set_xlabel("DOS", fontsize=12)
    ax_dos.set_title("Density of States", fontsize=14)
    ax_dos.axhline(0, color='black', linestyle='--', lw=0.5)
    ax_dos.grid(True, linestyle='--', alpha=0.5)    
    ax_dos.tick_params(axis='both', which='both', labelsize=12)
    

    for proj in projections_list:
        ax_dos.plot(Dosions_projup[ion][proj]*2, Energy_dos - E_fermi, color=colors[proj])
        ax_dos.plot(-Dosions_projdown[ion][proj]*2, Energy_dos - E_fermi, color=colors[proj],ls="--")

    # Create a handle with the color corresponding to the projection.
    legend_handles = []
    up_handle = mlines.Line2D([], [], color='black', linestyle='-',
                          label=f"spin maj")
    down_handle = mlines.Line2D([], [], color='black', linestyle='--',
                            label=f"spin min")
    legend_handles.extend([up_handle, down_handle])
    ax_dos.legend(handles=legend_handles, loc='upper right', fontsize=10)
    ax_dos.set_xlim(-20,20)
    
    ax_dos.axvline(0, color='black', linestyle='-', lw=0.5)
    #----------------------- Plot the band structure ---------------------------
    # Loop over each projection and scatter plot the points where it is dominant.
    for proj in projections_list:
        for i in range(nband):
                ax_band_projup.scatter(kDis, kbandenergyup[:, i] - E_fermi,
                                       s=Bandions_projup[ion][proj][:, i] * marker_scale,
                                       color=colors[proj], marker='o', facecolors = 'none', alpha = 0.5)
                ax_band_projdown.scatter(kDis, kbandenergydown[:, i] - E_fermi,
                                       s=Bandions_projdown[ion][proj][:, i] * marker_scale,
                                       color=colors[proj], marker='o', facecolors = 'none', alpha = 0.5)

    # Plot high-symmetry vertical lines
    for pos in tick_positions:
        ax_band_projup.axvline(x=pos, color='grey', linestyle='--', lw=0.8)
        ax_band_projdown.axvline(x=pos, color='grey', linestyle='--', lw=0.8)

    ax_band_projup.set_xticks(tick_positions)
    ax_band_projup.set_xticklabels(tick_labels, fontsize=20)
    ax_band_projup.axhline(y=0, color='red', lw=1.5, linestyle='--', label=r'$E_F$')
    ax_band_projup.set_xlabel("k-path", fontsize=12)
    ax_band_projup.set_title("Band Structure majority spin", fontsize=14, loc = 'left')
    ax_band_projup.set_xlim(0, kDis[-1])
    ax_band_projup.grid(True, linestyle='--', alpha=0.5)
    ax_band_projup.tick_params(axis='both', which='both', labelsize=12)
    ax_band_projdown.set_xticks(tick_positions)
    ax_band_projdown.set_xticklabels(tick_labels, fontsize=20)
    ax_band_projdown.axhline(y=0, color='red', lw=1.5, linestyle='--', label=r'$E_F$')
    ax_band_projdown.set_xlabel("k-path", fontsize=12)
    ax_band_projdown.set_title("Band Structure minority spin", fontsize=14, loc = 'left')
    ax_band_projdown.set_xlim(0, kDis[-1])
    ax_band_projdown.grid(True, linestyle='--', alpha=0.5)
    ax_band_projdown.tick_params(axis='both', which='both', labelsize=12)

    # Optional: Create a custom legend for the projections.
    legend_handles = []

    # Add projection legend entries
    for proj in projections_list:
        handle = mlines.Line2D([], [], color=colors[proj], marker='o', linestyle='None',
                           markersize=10, label=f"{proj}")
        legend_handles.append(handle)

    # Add black solid and dashed lines for spin majority and minority
    # spin_maj_handle = mlines.Line2D([], [], color='black', linestyle='-', linewidth=2, label='Spin maj')
    # spin_min_handle = mlines.Line2D([], [], color='black', linestyle='--', linewidth=2, label='Spin min')

    # legend_handles.extend([spin_maj_handle, spin_min_handle])

    # Add the legend to the plot
    ax_band_projup.legend(handles=legend_handles, loc='upper right', fontsize=10)

    # Label the shared y-axis (energy) and set axes limits
    ax_band_projdown.set_ylabel("Energy (eV)", fontsize=12)
    # ax_band_projdown.set_ylabel("Energy (eV)", fontsize=12)
    ax_band_projdown.set_ylim(-1, 6)


#Saving and printing the figures
plt.tight_layout()
plt.show()