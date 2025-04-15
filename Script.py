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
print(blocksize_DOS)
#Extract ions list and ions numbers
ions_list = []
ions_number = []
for i in range(len(POSCAR[5].strip().split())):
    ions_list.append(str(POSCAR[5].strip().split()[i]))
    ions_number.append(int(POSCAR[6].strip().split()[i]))
projections_list = PROCAR[7].strip().split()[1:-1]


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
    'X': np.array([0.0, 0.5, 0.0]),
    'M': np.array([0.5, 0.5, 0.0]),
    'R': np.array([0.0, 0.5, 0.5]),
    #'K': np.array([0.375, 0.375, 0.75]),
    #'L': np.array([0.5, 0.5, 0.5])
}

# Enter here the desired order (even with repeated points)
path_order = ['Γ', 'X', 'M', 'Γ', 'R', 'X', 'R', 'M']
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



#******************************************************************************************
#************************************** DOS + BAND total **********************************
#******************************************************************************************
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

#----------------------------------------------------------------------------------
#------------------------------ Plotting DOS + BAND -------------------------------
#----------------------------------------------------------------------------------
fig_total, (ax_dos, ax_band) = plt.subplots(1, 2, sharey=True,
                                        gridspec_kw={'width_ratios': [1, 2]},
                                        figsize=(12, 8))

#----- DOS Plot (left) -----
# Fill the area under the DOS curves with grey
ax_dos.fill_betweenx(Energy_dos - E_fermi, 0, Dosp, facecolor='grey', alpha=0.3)
ax_dos.fill_betweenx(Energy_dos - E_fermi, -Dosm, 0, facecolor='grey', alpha=0.3)
# Plot the DOS lines on top
ax_dos.plot(Dosp, Energy_dos - E_fermi, color='black', lw=1, label='DOS (spin up)')
ax_dos.plot(-Dosm, Energy_dos - E_fermi, color='red', ls="--", lw=1, label='DOS (spin down)')
ax_dos.set_xlabel("DOS", fontsize=12)
ax_dos.set_title("Density of States", fontsize=14)
ax_dos.axhline(0, color='black', linestyle='--', lw=0.5)
ax_dos.legend(loc='upper left', fontsize=10)
ax_dos.invert_xaxis()  # Optional: Invert x-axis so that positive DOS appears on the right
ax_dos.grid(True, linestyle='--', alpha=0.5)
ax_dos.tick_params(axis='both', which='both', labelsize=10)

#----- Band Structure Plot (right) -----

for i in range(nband):
    if i == 0:
        ax_band.plot(kDis, kbandenergyup[:, i] - E_fermi, lw=1.0, color='black', label="spin up")
        ax_band.plot(kDis, kbandenergydown[:, i] - E_fermi, lw=1.0, ls="--", color='red', label="spin down")
    else:
        ax_band.plot(kDis, kbandenergyup[:, i] - E_fermi, lw=1.0, color='black')
        ax_band.plot(kDis, kbandenergydown[:, i] - E_fermi, lw=1.0, ls="--", color='red')

    
# Draw vertical dashed lines for high-symmetry points
for pos in tick_positions:
    ax_band.axvline(x=pos, color='grey', linestyle='--', lw=0.8)

ax_band.set_xticks(tick_positions)
ax_band.set_xticklabels(tick_labels, fontsize=10)
ax_band.axhline(y=0, color='red', lw=1.5, linestyle='--', label=r'$E_F$')
ax_band.set_xlabel("k-path", fontsize=12)
ax_band.set_title("Band Structure", fontsize=14)
ax_band.set_xlim(0, kDis[-1])
ax_band.legend(loc='upper right', fontsize=10)
ax_band.grid(True, linestyle='--', alpha=0.5)
ax_band.tick_params(axis='both', which='both', labelsize=10)

# Label the shared y-axis (energy)
ax_dos.set_ylabel("Energy (eV)", fontsize=12)
ax_dos.set_ylim(-10,5)
ax_dos.set_xlim(-10,10)



#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


#*******************************************************************************************
#************************** Ions contributions DOS + BAND **********************************
#*******************************************************************************************
#New array to stock the number of ions studied by species 
ions_numbersum = np.zeros_like(ions_number)
for i in range(len(ions_number)-1):
    ions_numbersum[i+1] = ions_numbersum[i]+ions_number[i]
ions_dosup={}
ions_dosdown={}
for ions in ions_list:
    ions_dosup[ions]=np.zeros(blocksize_DOS)
    ions_dosdown[ions]=np.zeros(blocksize_DOS)

index=0
# Loop over ions and energy blocks
for ions in ions_list:
    index += 1
    for i in range(blocksize_DOS):
        # Calculate the index for the current ion block; adjust as required by the file format
        index_block = 6 + (ions_numbersum[index-1]+1) * (blocksize_DOS+1) 
        sumup=0
        sumdown=0
        for j in range(int((linesize-1)/2)-1):
            for k in range(ions_numbersum[index-1]):
                sumup+=float(DOSCAR[index_block+(blocksize_DOS+1)*k+i].strip().split()[j*2+1])
                sumdown+=float(DOSCAR[index_block+(blocksize_DOS+1)*k+i].strip().split()[j*2+2])
        ions_dosup[ions][i]=sumup
        ions_dosdown[ions][i]=sumdown

#Extract the energy for each kpoints by ions of the PROCAR

# Initialize matrices for each ion.
ions_Bandup = {}
ions_Banddown = {}
for ions in ions_list:
    ions_Bandup[ions] = np.zeros((nkpoint, nband))
    ions_Banddown[ions] = np.zeros((nkpoint, nband))

# Fill in the data for each projection.
index = -1
for ions in ions_list:
    index += 1
    for i in range(nkpoint):
        for j in range(nband):
            Bandup = 0
            Banddown = 0
            for k in range(ions_numbersum[index]): 
                Bandup += float(PROCAR[8+i*(blocksize_PRO+3)+j*bandsize+k+ions_numbersum[index]].strip().split()[10])
                Banddown += float(PROCAR[9+nkpoint*(blocksize_PRO+3)+i*(blocksize_PRO+3)+j*bandsize+k+ions_numbersum[index]].strip().split()[10])
            ions_Bandup[ions][i,j]=Bandup
            ions_Banddown[ions][i,j]=Banddown


#--------------------------------------------------------------------------------------------
#---------------------------- Plotting ions contribution ------------------------------------
#--------------------------------------------------------------------------------------------

fig_ions, (ax_ionsdos, ax_ionsband) = plt.subplots(1, 2, sharey=True,
                                        gridspec_kw={'width_ratios': [1, 2]},
                                        figsize=(12, 8))
# Define a list of colors to cycle through (extend as needed) without extra spaces.
color_list = ['indianred', 'teal', 'purple']
colors = {}
for i, ions in enumerate(ions_list):
    colors[ions] = color_list[i % len(color_list)]


#Plotting the DOS for each contribution
# ----- DOS Plot (left) -----
# Fill the area under the DOS curves with grey
ax_ionsdos.fill_betweenx(Energy_dos - E_fermi, 0, Dosp, facecolor='grey', alpha=0.3)
ax_ionsdos.fill_betweenx(Energy_dos - E_fermi, -Dosm, 0, facecolor='grey', alpha=0.3)
# Plot the DOS lines on top
ax_ionsdos.plot(Dosp, Energy_dos - E_fermi, color='black', lw=1, label='DOS (spin up)')
ax_ionsdos.plot(-Dosm, Energy_dos - E_fermi, color='black', ls="--", lw=1, label='DOS (spin down)')
ax_ionsdos.set_xlabel("DOS", fontsize=12)
ax_ionsdos.set_title("Density of States", fontsize=14)
ax_ionsdos.axhline(0, color='black', linestyle='--', lw=0.5)
ax_ionsdos.grid(True, linestyle='--', alpha=0.5)
ax_ionsdos.tick_params(axis='both', which='both', labelsize=10)

# Add the contribution of each ion projection
# Loop over all ions and plot their contributions
for ions in ions_list:
    # Ensure that the length of the projection data matches the energy array length;
    # adjust indices if needed.
    ax_ionsdos.plot(ions_dosup[ions], Energy_dos - E_fermi, color=colors[ions],label=f"{ions} spin maj")
    ax_ionsdos.plot(-ions_dosdown[ions], Energy_dos - E_fermi, color=colors[ions],ls="--", label=f"{ions} spin min")
ax_ionsdos.legend(loc='upper left', fontsize=10)

marker_scale = 100  # Adjust as needed for visual clarity
# Plot the projected band structure for each projection.
for ions in ions_list:
    for i in range(nband):
        if i == 0:
            ax_ionsband.scatter(kDis, kbandenergyup[:, i] - E_fermi,
                                   s=ions_Bandup[ions][:, i] * marker_scale,
                                   color=colors[ions], marker='o', label=f"{ions}")
            # Uncomment the following block if you also want to plot spin down data.
            
            # ax_ionsband.scatter(kDis, kbandenergydown[:, i] - E_fermi,
            #                        s=ions_Banddown[ions][:, i] * marker_scale,
            #                        color=colors[ions], marker='s', label=f"{ions}")
            
        else:
            ax_ionsband.scatter(kDis, kbandenergyup[:, i] - E_fermi,
                                   s=ions_Bandup[ions][:, i] * marker_scale,
                                   color=colors[ions], marker='o')
            
            # ax_ionsband.scatter(kDis, kbandenergydown[:, i] - E_fermi,
            #                        s=ions_Banddown[ions][:, i] * marker_scale,
            #                        color=colors[ions], marker='s')
# Draw vertical dashed lines for high-symmetry points
for pos in tick_positions:
    ax_ionsband.axvline(x=pos, color='grey', linestyle='--', lw=0.8)

#parameters for u
ax_ionsband.set_xticks(tick_positions)
ax_ionsband.set_xticklabels(tick_labels, fontsize=10)
ax_ionsband.axhline(y=0, color='red', lw=1.5, linestyle='--', label=r'$E_F$')
ax_ionsband.set_xlabel("k-path", fontsize=12)
ax_ionsband.set_title("Band Structure majority spin", fontsize=14)
ax_ionsband.set_xlim(0, kDis[-1])
ax_ionsband.grid(True, linestyle='--', alpha=0.5)
ax_ionsband.tick_params(axis='both', which='both', labelsize=10)

# Create custom legend handles with a fixed markersize (e.g., 10)
legend_handles = []
for ions in ions_list:
    # Handle for spin up (circle)
    up_handle = mlines.Line2D([], [], color=colors[ions], marker='o', linestyle='None',
                              markersize=10, label=f"{ions} up")
    # Handle for spin down (square)
    #down_handle = mlines.Line2D([], [], color=colors[ions], marker='s', linestyle='None',
    #                            markersize=10, label=f"{ions} down")
    #legend_handles.extend([up_handle, down_handle])
    legend_handles.append(up_handle)
ax_ionsband.legend(handles=legend_handles, loc='upper right', fontsize=10)

# Label the shared y-axis (energy)
ax_ionsdos.set_ylabel("Energy (eV)", fontsize=12)
ax_ionsdos.set_ylim(-10, 5)
ax_ionsdos.set_xlim(-10,10)



#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



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


#---------------------------------------------------------------------------------------------------------------------
#-------------------------------------- Plotting projections contributions -------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

fig_proj, (ax_dos_proj, ax_band_proj) = plt.subplots(1, 2, sharey=True,
                                        gridspec_kw={'width_ratios': [1, 2]},
                                        figsize=(12, 8))

# Define a list of colors to cycle through (extend as needed) without extra spaces.
color_list = ['gray', 'teal', 'cyan', 'skyblue', 'rosybrown', 'lightcoral', 'indianred', 'brown', 'firebrick']
colors = {}
for i, proj in enumerate(projections_list):
    colors[proj] = color_list[i % len(color_list)]

marker_scale = 100  # Adjust as needed for visual clarity

#Plot the DOS left part
# Fill the area under the DOS curves with grey
ax_dos_proj.fill_betweenx(Energy_dos - E_fermi, 0, Dosp, facecolor='grey', alpha=0.3)
ax_dos_proj.fill_betweenx(Energy_dos - E_fermi, -Dosm, 0, facecolor='grey', alpha=0.3)
# Plot the DOS lines on top
ax_dos_proj.plot(Dosp, Energy_dos - E_fermi, color='black', lw=1, label='DOS (spin up)')
ax_dos_proj.plot(-Dosm, Energy_dos - E_fermi, color='black', ls="--", lw=1, label='DOS (spin down)')
ax_dos_proj.set_xlabel("DOS", fontsize=12)
ax_dos_proj.set_title("Density of States", fontsize=14)
ax_dos_proj.axhline(0, color='black', linestyle='--', lw=0.5)
ax_dos_proj.grid(True, linestyle='--', alpha=0.5)
ax_dos_proj.tick_params(axis='both', which='both', labelsize=10)
for proj in projections_list:
    ax_dos_proj.plot(dos_projup[proj], Energy_dos - E_fermi, color=colors[proj],label=f"{proj} spin maj")
    ax_dos_proj.plot(-dos_projdown[proj], Energy_dos - E_fermi, color=colors[proj],ls="--", label=f"{proj} spin min")
ax_dos_proj.legend(loc='upper left', fontsize=10)

# Plot the projected band structure for each projection.
# Using different marker shapes for spin up (circle) and spin down (square).
for proj in projections_list:
    for i in range(nband):
        if i == 0:
            ax_band_proj.scatter(kDis, kbandenergyup[:, i] - E_fermi,
                                   s=proj_Bandup[proj][:, i] * marker_scale,
                                   color=colors[proj], marker='o', label=f"{proj}")
            # Uncomment the following block if you also want to plot spin down data.
            
            # ax_band_proj.scatter(kDis, kbandenergydown[:, i] - E_fermi,
            #                        s=proj_Banddown[proj][:, i] * marker_scale,
            #                        color=colors[proj], marker='s', label=f"{proj}")
            
        else:
            ax_band_proj.scatter(kDis, kbandenergyup[:, i] - E_fermi,
                                   s=proj_Bandup[proj][:, i] * marker_scale,
                                   color=colors[proj], marker='o')
            
            # ax_band_proj.scatter(kDis, kbandenergydown[:, i] - E_fermi,
            #                        s=proj_Banddown[proj][:, i] * marker_scale,
            #                        color=colors[proj], marker='s')
            

# Draw vertical dashed lines for high-symmetry points
for pos in tick_positions:
    ax_band_proj.axvline(x=pos, color='grey', linestyle='--', lw=0.8)

#parameters for up
ax_band_proj.set_xticks(tick_positions)
ax_band_proj.set_xticklabels(tick_labels, fontsize=10)
ax_band_proj.axhline(y=0, color='red', lw=1.5, linestyle='--', label=r'$E_F$')
ax_band_proj.set_xlabel("k-path", fontsize=12)
ax_band_proj.set_title("Band Structure majority spin", fontsize=14)
ax_band_proj.set_xlim(0, kDis[-1])
ax_band_proj.grid(True, linestyle='--', alpha=0.5)
ax_band_proj.tick_params(axis='both', which='both', labelsize=10)

# Create custom legend handles with a fixed markersize (e.g., 10)
legend_handles = []
for proj in projections_list:
    # Handle for spin up (circle)
    up_handle = mlines.Line2D([], [], color=colors[proj], marker='o', linestyle='None',
                              markersize=10, label=f"{proj} up")
    # Handle for spin down (square)
    # down_handle = mlines.Line2D([], [], color=colors[proj], marker='s', linestyle='None',
    #                             markersize=10, label=f"{proj} down")
    # legend_handles.extend([up_handle, down_handle])
    legend_handles.append(up_handle)

ax_band_proj.legend(handles=legend_handles, loc='upper right', fontsize=10)
# Label the shared y-axis (energy)
ax_dos_proj.set_ylabel("Energy (eV)", fontsize=12)
ax_dos_proj.set_ylim(-10, 5)
ax_dos_proj.set_xlim(-10,10)




#Saving and printing the figures
plt.tight_layout()

plt.show()