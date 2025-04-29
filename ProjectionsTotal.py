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
    if (np.sqrt((kcoord[i+1,0] - kcoord[i,0])**2 + (kcoord[i+1,1] - kcoord[i,1])**2 + (kcoord[i+1,2] - kcoord[i,2])**2)>0.2):
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
#****************************** PROJECTIONS by ions *************************************
#*******************************************************************************************
# # Create the arrays to stock the values
# Dosions_projup = {ions: {proj: np.zeros(blocksize_DOS) for proj in projections_list} for ions in ions_list}
# Dosions_projdown = {ions: {proj: np.zeros(blocksize_DOS) for proj in projections_list} for ions in ions_list}

# #Loop for the filling
# index=-1
# for ions in ions_list:
#     index += 1
#     indexproj = 0
#     for proj in projections_list:
#         indexproj += 1
#         for i in range(blocksize_DOS):
#             sumup = 0
#             sumdown = 0
#             for j in range(ions_number[index]):
#                 index_block = 6 + (j+1+ions_numbersum[index]) * (blocksize_DOS+1) + i 
#                 sumup+=float(DOSCAR[index_block].strip().split()[indexproj*2-1])
#                 sumdown+=float(DOSCAR[index_block].strip().split()[indexproj*2])
#             Dosions_projup[ions][proj][i]=sumup
#             Dosions_projdown[ions][proj][i]=sumdown


# #Extract the PROCAR values
# Bandions_projup = {ions: {proj: np.zeros((nkpoint, nband)) for proj in projections_list} for ions in ions_list}
# Bandions_projdown = {ions: {proj: np.zeros((nkpoint, nband)) for proj in projections_list} for ions in ions_list}

# index = -1
# for ions in ions_list:
#     index += 1
#     proj_index = 0
#     for proj in projections_list:  
#         proj_index += 1  
#         for i in range(nkpoint):
#             for j in range(nband):
#                 projup = 0
#                 projdown = 0
#                 for k in range(ions_number[index]):
#                     projup += float(PROCAR[8+i*(blocksize_PRO+3)+j*bandsize+k+ions_numbersum[index]].strip().split()[proj_index])
#                     projdown += float(PROCAR[9+nkpoint*(blocksize_PRO+3)+i*(blocksize_PRO+3)+j*bandsize+k+ions_numbersum[index]].strip().split()[proj_index])
#                 Bandions_projup[ions][proj][i][j]=projup
#                 Bandions_projdown[ions][proj][i][j]=projdown


#---------------------------------------------------------------------------------------------------------------------
#---------------------------------- Plotting projections contributions for each ions ---------------------------------
#---------------------------------------------------------------------------------------------------------------------

# # Define a list of colors to cycle through (extend as needed)
# color_list = ['gray', 'skyblue', 'skyblue', 'skyblue', 'indianred', 'indianred', 'indianred', 'indianred', 'indianred']
# # Create a dictionary for colors for each projection (this will be re-used for each ion)
# colors = {}
# for i, proj in enumerate(projections_list):
#     colors[proj] = color_list[i % len(color_list)]
# projections_list2 = projections_list
# projections_list = projections_list[1:9]
# # Marker size scaling factor for the scatter plots
# marker_scale = 10

# # For each ion, create a figure with DOS and band structure subplots.
# for ion in ions_list:
#     # Create a new figure for each ion
#     fig_proj, (ax_dos_proj, ax_band_proj) = plt.subplots(1, 2, sharey=True,
#                                             gridspec_kw={'width_ratios': [1, 2]},
#                                             figsize=(12, 8))
#     fig_proj.suptitle(f"Ion: {ion}", fontsize=16)

#     #----------------------- Plot the DOS contributions ---------------------------
#     # Fill the area under the DOS curves with grey (assuming these variables are global,
#     # i.e. they are the same for each ion; if they are ion specific, you should replace them accordingly)
#     ax_dos_proj.fill_betweenx(Energy_dos - E_fermi, 0, Dosp, facecolor='grey', alpha=0.3)
#     ax_dos_proj.fill_betweenx(Energy_dos - E_fermi, -Dosm, 0, facecolor='grey', alpha=0.3)
#     # Plot the DOS lines on top
#     ax_dos_proj.plot(Dosp, Energy_dos - E_fermi, color='black', lw=1, label='DOS (spin maj)')
#     ax_dos_proj.plot(-Dosm, Energy_dos - E_fermi, color='black', ls="--", lw=1, label='DOS (spin min)')
#     ax_dos_proj.set_xlabel("DOS", fontsize=12)
#     ax_dos_proj.set_title("Density of States", fontsize=14)
#     ax_dos_proj.axhline(0, color='black', linestyle='--', lw=0.5)
#     ax_dos_proj.grid(True, linestyle='--', alpha=0.5)
#     ax_dos_proj.tick_params(axis='both', which='both', labelsize=10)
    
#     # Plot projection contributions for each projection on the DOS panel:
#     # (Assuming dos_projup and dos_projdown are nested dictionaries: dos_projup[ion][proj])
#     for proj in projections_list:
#         ax_dos_proj.plot(Dosions_projup[ion][proj], Energy_dos - E_fermi, color=colors[proj],
#                          label=f"{proj} spin maj")
#         ax_dos_proj.plot(-Dosions_projdown[ion][proj], Energy_dos - E_fermi, color=colors[proj],
#                          ls="--", label=f"{proj} spin min")
#     ax_dos_proj.legend(loc='upper left', fontsize=10)

#     #----------------------- Plot the band structure ---------------------------
#     # Number of k-points and bands
#     num_k = len(kDis)
#     num_bands = nband
    
#     # Here, proj_Bandup is assumed to be a nested dictionary with structure: proj_Bandup[ion][proj]
#     weights_arrayup = np.stack([Bandions_projup[ion][proj] for proj in projections_list], axis=-1)
#     max_indicesup = np.argmax(weights_arrayup, axis=-1) 
#     weights_arraydown = np.stack([Bandions_projdown[ion][proj] for proj in projections_list], axis=-1)
#     max_indicesdown = np.argmax(weights_arraydown, axis=-1) 
#     # Prepare arrays of k-points and adjusted energies
#     k_points = np.tile(np.array(kDis).reshape(-1, 1), (1, num_bands))
#     energy_pointsup = kbandenergyup - E_fermi   # adjust energies by the Fermi energy
#     energy_pointsdown = kbandenergydown - E_fermi   # adjust energies by the Fermi energy

#     # Loop over each projection and scatter plot the points where it is dominant.
#     for i, proj in enumerate(projections_list):
#         mask = (max_indicesup == i)
#         if np.any(mask):
#             ax_band_proj.scatter(k_points[mask],
#                                  energy_pointsup[mask],
#                                  s=(Bandions_projup[ion][proj][mask] * marker_scale),
#                                  color=colors[proj],
#                                  marker='o')
    
#     # Loop over each projection and scatter plot the points where it is dominant.
#     for i, proj in enumerate(projections_list):
#         mask = (max_indicesdown == i)
#         if np.any(mask):
#             ax_band_proj.scatter(k_points[mask],
#                                  energy_pointsup[mask],
#                                  s=(Bandions_projdown[ion][proj][mask] * marker_scale),
#                                  color=colors[proj],
#                                  marker='s')
    
#     # Plot high-symmetry vertical lines
#     for pos in tick_positions:
#         ax_band_proj.axvline(x=pos, color='grey', linestyle='--', lw=0.8)
#     ax_band_proj.set_xticks(tick_positions)
#     ax_band_proj.set_xticklabels(tick_labels, fontsize=10)
#     ax_band_proj.axhline(y=0, color='red', lw=1.5, linestyle='--', label=r'$E_F$')
#     ax_band_proj.set_xlabel("k-path", fontsize=12)
#     ax_band_proj.set_title("Band Structure majority spin", fontsize=14)
#     ax_band_proj.set_xlim(0, kDis[-1])
#     ax_band_proj.grid(True, linestyle='--', alpha=0.5)
#     ax_band_proj.tick_params(axis='both', which='both', labelsize=10)

#     # Optional: Create a custom legend for the projections.
#     legend_handles = []
#     for proj in projections_list:
#         handle = mlines.Line2D([], [], color=colors[proj], marker='o', linestyle='None',
#                                  markersize=10, label=f"{proj}")
#         legend_handles.append(handle)
#     ax_band_proj.legend(handles=legend_handles, loc='upper right', fontsize=10)

#     # Label the shared y-axis (energy) and set axes limits
#     ax_dos_proj.set_ylabel("Energy (eV)", fontsize=12)
#     ax_dos_proj.set_ylim(-5, 5)
#     ax_dos_proj.set_xlim(-6, 6)




#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# projections_list = projections_list2
#*******************************************************************************************
#****************************** PROJECTIONS DOS + BAND *************************************
#*******************************************************************************************
#Extract the DOSCAR values
# Create the arrays to stock the values
dos_projup={}
dos_projdown={}
for proj in projections_list:
    dos_projup[proj]=np.zeros(blocksize_DOS)
    dos_projdown[proj]=np.zeros(blocksize_DOS)

#Loop for the filling
index=0
for proj in projections_list:
    index += 1
    for i in range(blocksize_DOS):
        sumup = 0
        sumdown = 0
        for j in range(nions-1):
            index_block = 6 + (j+1) * (blocksize_DOS+1) + i 
            sumup+=float(DOSCAR[index_block].strip().split()[index*2-1])
            sumdown+=float(DOSCAR[index_block].strip().split()[index*2])
        dos_projup[proj][i]=sumup
        dos_projdown[proj][i]=sumdown

#Extract the PROCAR values
# Initialize matrices for all projections.
proj_Bandup = {}
proj_Banddown = {}
for proj in projections_list:
    proj_Bandup[proj] = np.zeros((nkpoint, nband))
    proj_Banddown[proj] = np.zeros((nkpoint, nband))

# Fill in the data for each projection.
index = 0
for proj in projections_list:
    index += 1  
    for i in range(nkpoint):
        for j in range(nband):
            proj_Bandup[proj][i, j] = float(PROCAR[3+bandsize+i*(blocksize_PRO+3)+j*bandsize].strip().split()[index])
            proj_Banddown[proj][i, j] = float(PROCAR[4+bandsize+nkpoint*(blocksize_PRO+3)+i*(blocksize_PRO+3)+j*bandsize].strip().split()[index])

projection_listsum = ['psum', 'dsum'] 
dos_projup['psum'] = dos_projup['px'] + dos_projup['py'] + dos_projup['px']
dos_projdown['psum'] = dos_projdown['px'] + dos_projdown['py'] + dos_projdown['px']
dos_projup['dsum'] = dos_projup['dxy'] + dos_projup['dyz'] + dos_projup['dxz'] + dos_projup['x2-y2'] + dos_projup['dz2']
dos_projdown['dsum'] = dos_projdown['dxy'] + dos_projdown['dyz'] + dos_projdown['dxz'] + dos_projdown['x2-y2'] + dos_projdown['dz2']

# #inverting the contribution for projection dxy and dx2-y2 for the structure Pnma
# projBandupdxy = proj_Bandup['dxy']
# projBandupdx2y2 = proj_Bandup['x2-y2']
# projBanddowndxy = proj_Banddown['dxy']  
# projBanddowndx2y2 = proj_Banddown['x2-y2'] 
# proj_Bandup['dxy'] = projBandupdx2y2
# proj_Bandup['x2-y2'] = projBandupdxy
# proj_Banddown['dxy'] = projBanddowndx2y2
# proj_Banddown['x2-y2'] = projBanddowndxy

# #inverting the contribution for projection dxy and dx2-y2 for the structure Pnma
# projdosupdxy = dos_projup['dxy']
# projdosupdx2y2 = dos_projup['x2-y2']
# projdosdowndxy = dos_projdown['dxy']  
# projdosdowndx2y2 = dos_projdown['x2-y2'] 
# dos_projup['dxy'] = projdosupdx2y2
# dos_projup['x2-y2'] = projdosupdxy
# dos_projdown['dxy'] = projdosdowndx2y2
# dos_projdown['x2-y2'] = projdosdowndxy


# proj_Bandup['dxz+dyz'] = proj_Bandup['dxz'] + proj_Bandup['dyz'] 
# proj_Banddown['dxz+dyz'] = proj_Banddown['dxz'] + proj_Banddown['dyz'] 

# dos_projup['dxz+dyz'] = dos_projup['dxz'] + dos_projup['dyz'] 
# dos_projdown['dxz+dyz'] = dos_projdown['dxz'] + dos_projdown['dyz'] 

# projections_list = ['s', 'py', 'px', 'pz', 'dxy', 'dxz+dyz', 'dz2', 'x2-y2']

#---------------------------------------------------------------------------------------------------------------------
#-------------------------------------- Plotting projections contributions + DOS -------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

# fig_proj, (ax_band_down, ax_dos, ax_band_up) = plt.subplots(1, 3, sharey=True,
#                                         # gridspec_kw={'width_ratios': [3, 1, 3]},
#                                         figsize=(12, 8))
# fig_proj.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.1, wspace=0.1)

fig_proj, (ax_dos, ax_band_up) = plt.subplots(1, 2, sharey=True,
                                        gridspec_kw={'width_ratios': [1,2]},
                                        figsize=(12, 8))

# Define a list of colors to cycle through (extend as needed) without extra spaces.
# color_list = ['gray', 'teal', 'cyan', 'skyblue', 'rosybrown', 'lightcoral', 'indianred', 'brown', 'firebrick']
color_list = ['gray', 'cyan', 'cyan', 'cyan', 'indianred', 'indianred', 'indianred', 'indianred', 'indianred']
colors = {}

for i, proj in enumerate(projections_list):
    colors[proj] = color_list[i % len(color_list)]

colorslsitsum = ['cyan', 'red']
colorssum = {}
for i, proj in enumerate(projection_listsum):
    colorssum[proj] = colorslsitsum[i % len(colorslsitsum)]

projections_list2 = projections_list
projections_list = projections_list[1:9]
marker_scale = 100  # Adjust as needed for visual clarity

#----------------------- Plot the DOS contributions ---------------------------
# Fill the area under the DOS curves with grey (assuming these variables are global,
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
    

for proj in projection_listsum:
    ax_dos.plot(dos_projup[proj], Energy_dos - E_fermi, color=colorssum[proj])
    ax_dos.plot(-dos_projdown[proj], Energy_dos - E_fermi, color=colorssum[proj],ls="--")

# Create a handle with the color corresponding to the projection.
legend_handles = []
up_handle = mlines.Line2D([], [], color='black', linestyle='-',
                          label=f"spin maj projection")
down_handle = mlines.Line2D([], [], color='black', linestyle='--',
                            label=f"spin min projection")
legend_handles.extend([up_handle, down_handle])
ax_dos.legend(handles=legend_handles, loc='lower right', fontsize=10)
ax_dos.set_xlim(-6,6)


# Number of k-points and bands (assuming these are defined)
num_k = len(kDis)
num_bands = nband

# Each slice along the third axis corresponds to one orbital’s projection weight for maj spins.
weights_arrayup = np.stack([proj_Bandup[proj] for proj in projections_list], axis=-1) 
max_indicesup = np.argmax(weights_arrayup, axis=-1)  # shape: (num_k, num_bands)
k_points = np.tile(np.array(kDis).reshape(-1, 1), (1, num_bands))
energy_pointsup = kbandenergyup - E_fermi
# Each slice along the third axis corresponds to one orbital’s projection weight for min spins.
weights_arraydown = np.stack([proj_Banddown[proj] for proj in projections_list], axis=-1) 
max_indicesdown = np.argmax(weights_arraydown, axis=-1)  # shape: (num_k, num_bands)
energy_pointsdown = kbandenergydown - E_fermi

# Loop over each projection in the defined list and plot the points where it is dominant
for i, proj in enumerate(projections_list):
    # Create a boolean mask for the points where the current projection is dominant
    mask = (max_indicesup == i)
    # Proceed with plotting only if there are points for this projection
    if np.any(mask):
        # Scatter plot the points in one call for the given projection
        ax_band_up.scatter(k_points[mask],
                             energy_pointsup[mask], ls = "-",
                            #  s=(proj_Bandup[proj][mask] * marker_scale),
                             s=1, 
                             color=colors[proj],
                             marker='o')

# Loop over each projection in the defined list and plot the points where it is dominant
for i, proj in enumerate(projections_list):
    # Create a boolean mask for the points where the current projection is dominant
    mask = (max_indicesdown == i)
    # Proceed with plotting only if there are points for this projection
    if np.any(mask):
        # Scatter plot the points in one call for the given projection
        ax_band_up.scatter(k_points[mask],
                             energy_pointsdown[mask], ls = "--",
                            #  s=(proj_Banddown[proj][mask] * marker_scale), 
                             s=1,
                             color=colors[proj],
                             marker='x')

# Then proceed with plotting high-symmetry lines, labels:
for pos in tick_positions:
    ax_band_up.axvline(x=pos, color='grey', linestyle='--', lw=0.8)
    # ax_band_down.axvline(x=pos, color='grey', linestyle='--', lw=0.8)

ax_band_up.set_xticks(tick_positions)
ax_band_up.set_xticklabels(tick_labels, fontsize=12)
ax_band_up.axhline(y=0, color='red', lw=1.5, linestyle='--', label=r'$E_F$')
ax_band_up.set_xlabel("k-path", fontsize=12)
ax_band_up.set_title("Band Structure majority spin", fontsize=14)
ax_band_up.set_xlim(0, kDis[-1])
ax_band_up.grid(True, linestyle='--', alpha=0.5)
ax_band_up.tick_params(axis='both', which='both', labelsize=12)

# ax_band_down.set_xticks(tick_positions)
# ax_band_down.set_xticklabels(tick_labels, fontsize=10)
# ax_band_down.axhline(y=0, color='red', lw=1.5, linestyle='--', label=r'$E_F$')
# ax_band_down.set_xlabel("k-path", fontsize=12)
# ax_band_down.set_title("Band Structure minority spin", fontsize=14)
# ax_band_down.set_xlim(0, kDis[-1])
# ax_band_down.grid(True, linestyle='--', alpha=0.5)
# ax_band_down.tick_params(axis='both', which='both', labelsize=12)

# Create a custom legend for the projections that ended up being dominant
# (You might need to adapt this if you want a legend entry only if at least one point comes from that projection.)
legend_handlesup = []
legend_handlesdown = []
projlist = ['p maj', 'p min', 'd maj', 'd min']
colorlist2 = ['blue', 'red', 'green','orange']
colors2={}
for i, proj in enumerate(projlist):
    colors[proj] = colorlist2[i % len(colorlist2)]
for proj in projlist:
    # Create a handle with the color corresponding to the projection.
    up_handle = mlines.Line2D([], [], color=colors[proj], marker='o', linestyle='None',
                              markersize=10, label=f"{proj}-orb")
    down_handle = mlines.Line2D([], [], color=colors[proj], marker='o', linestyle='None',
                              markersize=10, label=f"{proj}-orb")
    legend_handlesup.append(up_handle)
    legend_handlesdown.append(down_handle)

ax_band_up.legend(handles=legend_handlesup, loc='upper right', fontsize=10)
# ax_band_down.legend(handles=legend_handlesdown, loc='upper right', fontsize=10)

# Label the shared y-axis (energy)
ax_dos.set_ylabel("Energy (eV)", fontsize=12)
ax_dos.set_ylim(-3, 8)


#---------------------------------------------------------------------------------------------------------------------
#-------------------------------------- Plotting projections d and p separ -------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

# fig_proj, ((ax_band_dup, ax_band_ddown), (ax_band_pup, ax_band_pdown)) = plt.subplots(2,2, sharey=False,
#                                         figsize=(12, 8))

# # Define a list of colors to cycle through (extend as needed) without extra spaces.
# # color_list = ['gray', 'teal', 'cyan', 'skyblue', 'rosybrown', 'lightcoral', 'indianred', 'brown', 'firebrick']
# color_list = ['gray', 'cyan', 'green', 'indianred', 'cyan', 'purple', 'mediumseagreen', 'sandybrown', 'indianred']
# colors = {}

# projections_list = projections_list2

# for i, proj in enumerate(projections_list):
#     colors[proj] = color_list[i % len(color_list)]
# projections_list = projections_list[1:4]
# marker_scale = 100  # Adjust as needed for visual clarity

# for proj in projections_list:
#     for i in range(nband):
#             ax_band_pup.scatter(kDis, kbandenergyup[:, i] - E_fermi,
#                                    s=proj_Bandup[proj][:, i] * marker_scale,
#                                    color=colors[proj], marker='o', alpha = 0.8)
#             ax_band_pdown.scatter(kDis, kbandenergydown[:, i] - E_fermi,
#                                    s=proj_Banddown[proj][:, i] * marker_scale,
#                                    color=colors[proj], marker='o', facecolors = 'none', alpha = 0.5)

# # Then proceed with plotting high-symmetry lines, labels:
# for pos in tick_positions:
#     ax_band_pup.axvline(x=pos, color='grey', linestyle='--', lw=0.8)
#     ax_band_pdown.axvline(x=pos, color='grey', linestyle='--', lw=0.8)

# ax_band_pup.set_xticks(tick_positions)
# ax_band_pup.set_xticklabels(tick_labels, fontsize=10)
# ax_band_pup.axhline(y=0, color='red', lw=1.5, linestyle='--', label=r'$E_F$')
# ax_band_pup.set_xlabel("k-path", fontsize=12)
# ax_band_pup.set_title("p-projection majority spin", fontsize=14)
# ax_band_pup.set_xlim(0, kDis[-1])
# ax_band_pup.grid(True, linestyle='--', alpha=0.5)
# ax_band_pup.tick_params(axis='both', which='both', labelsize=10)

# ax_band_pdown.set_xticks(tick_positions)
# ax_band_pdown.set_xticklabels(tick_labels, fontsize=10)
# ax_band_pdown.axhline(y=0, color='red', lw=1.5, linestyle='--', label=r'$E_F$')
# ax_band_pdown.set_xlabel("k-path", fontsize=12)
# ax_band_pdown.set_title("p-projection minority spin", fontsize=14)
# ax_band_pdown.set_xlim(0, kDis[-1])
# ax_band_pdown.grid(True, linestyle='--', alpha=0.5)
# ax_band_pdown.tick_params(axis='both', which='both', labelsize=10)

# # Create a custom legend for the projections that ended up being dominant
# # (You might need to adapt this if you want a legend entry only if at least one point comes from that projection.)
# legend_handlesup = []
# legend_handlesdown = []


# for proj in projections_list:
#     # Create a handle with the color corresponding to the projection.
#     up_handle = mlines.Line2D([], [], color=colors[proj], marker='o', linestyle='None',
#                               markersize=10, label=f"{proj}")
#     down_handle = mlines.Line2D([], [], color=colors[proj], marker='o', linestyle='None',
#                               markersize=10, label=f"{proj}")
#     legend_handlesup.append(up_handle)
#     legend_handlesdown.append(down_handle)

# ax_band_pup.legend(handles=legend_handlesup, loc='upper right', fontsize=10)
# ax_band_pdown.legend(handles=legend_handlesdown, loc='upper right', fontsize=10)

# # Label the shared y-axis (energy)
# ax_band_pup.set_ylabel("Energy (eV)", fontsize=12)
# ax_band_pup.set_ylim(-6, 6)
# ax_band_pdown.set_ylim(-6, 6)


# projections_list=projections_list2
# marker_scale = 100  # Adjust as needed for visual clarity
# projections_list = projections_list[4:9]


# for proj in projections_list:
#     for i in range(nband):
#             ax_band_dup.scatter(kDis, kbandenergyup[:, i] - E_fermi,
#                                    s=proj_Bandup[proj][:, i] * marker_scale,
#                                    color=colors[proj], marker='o', alpha = 0.8)
#             ax_band_ddown.scatter(kDis, kbandenergydown[:, i] - E_fermi,
#                                    s=proj_Banddown[proj][:, i] * marker_scale,

#                                    color=colors[proj], marker='o', alpha = 0.5)

# # Then proceed with plotting high-symmetry lines, labels:
# for pos in tick_positions:
#     ax_band_dup.axvline(x=pos, color='grey', linestyle='--', lw=0.8)
#     ax_band_ddown.axvline(x=pos, color='grey', linestyle='--', lw=0.8)

# ax_band_dup.set_xticks(tick_positions)
# ax_band_dup.set_xticklabels(tick_labels, fontsize=10)
# ax_band_dup.axhline(y=0, color='red', lw=1.5, linestyle='--', label=r'$E_F$')
# ax_band_dup.set_xlabel("k-path", fontsize=12)
# ax_band_dup.set_title("d-projection", fontsize=14)
# ax_band_dup.set_xlim(0, kDis[-1])
# ax_band_dup.grid(True, linestyle='--', alpha=0.5)
# ax_band_dup.tick_params(axis='both', which='both', labelsize=10)

# ax_band_ddown.set_xticks(tick_positions)
# ax_band_ddown.set_xticklabels(tick_labels, fontsize=10)
# ax_band_ddown.axhline(y=0, color='red', lw=1.5, linestyle='--', label=r'$E_F$')
# ax_band_ddown.set_xlabel("k-path", fontsize=12)
# ax_band_ddown.set_title("d-projection minority spin", fontsize=14)
# ax_band_ddown.set_xlim(0, kDis[-1])
# ax_band_ddown.grid(True, linestyle='--', alpha=0.5)
# ax_band_ddown.tick_params(axis='both', which='both', labelsize=10)

# # Create a custom legend for the projections that ended up being dominant
# # (You might need to adapt this if you want a legend entry only if at least one point comes from that projection.)
# legend_handlesup = []
# legend_handlesdown = []

# for proj in projections_list:
#     # Create a handle with the color corresponding to the projection.
#     up_handle = mlines.Line2D([], [], color=colors[proj], marker='o', linestyle='None',
#                               markersize=10, label=f"{proj}")
#     down_handle = mlines.Line2D([], [], color=colors[proj], marker='o', linestyle='None',
#                               markersize=10, label=f"{proj}")
#     legend_handlesup.append(up_handle)
#     legend_handlesdown.append(down_handle)

# ax_band_dup.legend(handles=legend_handlesup, loc='upper right', fontsize=10)
# ax_band_ddown.legend(handles=legend_handlesdown, loc='upper right', fontsize=10)

# # Label the shared y-axis (energy)
# ax_band_dup.set_ylabel("Energy (eV)", fontsize=12)
# ax_band_dup.set_ylim(-3, 7)
# ax_band_ddown.set_ylim(-1, 7)
#Saving and printing the figures
plt.tight_layout()
plt.show()