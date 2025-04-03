import numpy as np
import matplotlib.pyplot as plt
import re
import itertools
import matplotlib.lines as mlines
import matplotlib.colors as mcolors

#*************************************
#***********Block for DOSCAR**********
#*************************************
# Read the PROCAR file
f = open("PROCAR", "r")
data = f.readlines()
#Extract the numlber of kpt band and ions
nkpoint=int(data[1].strip().split()[3])
nband=int(data[1].strip().split()[7])
nions=int(data[1].strip().split()[11])

# Read the DOSCAR file
f2=open('DOSCAR', 'r')
data_dos=f2.readlines()


# ----------- Block to extract the energy and DOS values from the DOSCAR file -----------
# Extract the Fermi energy from OUTCAR
#E_fermi = float([line.split()[2] for line in open("DOSCAR") if "E-fermi" in line][-1])
E_fermi = float(data_dos[5].strip().split()[3])
nions_dos=int(data_dos[0].strip().split()[1])

Energy=[]
DOSp=[]
DOSm=[]

# Extract the energy and DOS values from the DOSCAR file
for i in range(int((len(data_dos)-5)/(nions_dos+1))-1):
    index=6+i
    eigenval = float(data_dos[index].strip().split()[0])
    doseigenp= float(data_dos[index].strip().split()[1])
    doseigenm= float(data_dos[index].strip().split()[2])
    Energy.append(eigenval)
    DOSp.append(doseigenp)
    DOSm.append(doseigenm)

# Convert lists to numpy arrays
Energy_dos=np.array(Energy)
DOSp=np.array(DOSp)
DOSm=np.array(DOSm)

#-------------- Block to extract the DOS of each ions --------------
#Write the ions in order
ions_list =  ['Ni', 'O']
blocksize = int((len(data_dos) - 5) / (nions_dos + 1))

# Get the number of columns in one of the blocks (adjust the index if needed)
linesize = len(data_dos[blocksize + 8].strip().split())
ions_dosup={}
ions_dosdown={}
for ions in ions_list:
    ions_dosup[ions]=np.zeros(blocksize-1)
    ions_dosdown[ions]=np.zeros(blocksize-1)

index=0
# Loop over ions and energy blocks
for ions in ions_list:
    index += 2
    for i in range(blocksize - 1):
        # Calculate the index for the current ion block; adjust as required by the file format
        index_block = 6 + (index) * blocksize + i   
        sumup=0
        sumdown=0
        for j in range(int((linesize-1)/2)-1):
            sumup+=float(data_dos[index_block].strip().split()[j*2+1])
            sumdown+=float(data_dos[index_block].strip().split()[j*2+2])
        ions_dosup[ions][i]=sumup
        ions_dosdown[ions][i]=sumdown

#********************* block for each projection contributions *************************

# Extract the projection name
projection_line = data[7].strip().split()
projections_list = projection_line[1:-1]
blocksize = int((len(data_dos) - 5) / (nions_dos + 1))

# Get the number of columns in one of the blocks (adjust the index if needed)
linesize = len(data_dos[blocksize + 8].strip().split())
dos_projup={}
dos_projdown={}
for proj in projections_list:
    dos_projup[proj]=np.zeros(blocksize-1)
    dos_projdown[proj]=np.zeros(blocksize-1)

index=0
for proj in projections_list:
    index += 1
    for i in range(blocksize -1):
        sumup = 0
        sumdown = 0
        index_block=0
        for j in range(nions-1):
            index_block = 6 + (j+1) * blocksize + i 
            sumup+=float(data_dos[index_block].strip().split()[index*2-1])
            sumdown+=float(data_dos[index_block].strip().split()[index*2])
        dos_projup[proj][i]=sumup
        dos_projdown[proj][i]=sumdown


#*************************************
#**********Block for PROCAR***********
#*************************************

# Read the PROCAR file
f = open("PROCAR", "r")
data = f.readlines()

#Extract the numlber of kpt band and ions
nkpoint=int(data[1].strip().split()[3])
nband=int(data[1].strip().split()[7])
nions=int(data[1].strip().split()[11])

#extract the energy of each band for each k points
#For spin up
bandsize = 5+nions
blocksize = nband*bandsize
kbandenergy=np.zeros([nkpoint, nband])
if (nions==1):
    for i in range(nkpoint):
        for j in range(nband):
            kbandenergy[i,j]=float(data[5+i*((nband*5)+3)+j*5].strip().split()[4])
else :
    for i in range(nkpoint):
        for j in range(nband):
            kbandenergy[i,j]=float(data[5+i*(blocksize+3)+j*bandsize].strip().split()[4])
kbandenergyup=np.array(kbandenergy)

#For spin down
kbandenergy=np.zeros([nkpoint, nband])
if (nions==1):
    for i in range(nkpoint):
        for j in range(nband):
            kbandenergy[i,j]=float(data[6+nkpoint*(nband*5+3)+i*((nband*5)+3)+j*(5+nions)].strip().split()[4])
else :
    for i in range(nkpoint):
        for j in range(nband):
            kbandenergy[i,j]=float(data[6+nkpoint*(blocksize+3)+i*(blocksize+3)+j*bandsize].strip().split()[4])
kbandenergydown = np.array(kbandenergy)


#**************Block to extract for each projection the total contributions********************
# Extract the projection name
projection_line = data[7].strip().split()
projections = projection_line[1:-1]

#------------------For spin up----------------------
# Initialize matrices for all projections.
proj_matricesup = {}
for proj in projections:
    proj_matricesup[proj] = np.zeros((nkpoint, nband))

# Fill in the data for each projection.
index = 0
for proj in projections:
    index += 1  
    if nions == 1:
        for i in range(nkpoint):
            for j in range(nband):
                # Modify the index expression based on the projection, here as an example using the same index.
                proj_matricesup[proj][i, j] = float(data[13 + nkpoint*(blocksize+3) + i * ((nband * 5) + 3) + j * (5 + nions)].strip().split()[index])
    else:
        for i in range(nkpoint):
            for j in range(nband):
                proj_matricesup[proj][i, j] = float(data[12+i*(blocksize+3)+j*bandsize].strip().split()[index])

#-----------------For spin down--------------------
# Initialize matrices for all projections.
proj_matricesdown = {}
for proj in projections:
    proj_matricesdown[proj] = np.zeros((nkpoint, nband))

# Fill in the data for each projection.
index = 0
for proj in projections:
    index += 1  
    if nions == 1:
        for i in range(nkpoint):
            for j in range(nband):
                # Modify the index expression based on the projection, here as an example using the same index.
                proj_matricesdown[proj][i, j] = float(data[13+nkpoint*(blocksize+3)+i*((nband * 5)+3)+j*(5+nions)].strip().split()[index])
    else:
        for i in range(nkpoint):
            for j in range(nband):
                proj_matricesdown[proj][i, j] = float(data[13+nkpoint*(blocksize+3)+i*(blocksize+3)+j*bandsize].strip().split()[index])


#******************* Block to extract the contribution of each ions **********************
#------------------For spin up----------------------
#Write the ions in order
ions_list =  ['Ni', 'O']
# Initialize matrices for each ion.
ions_contribup = {}
for ions in ions_list:
    ions_contribup[ions] = np.zeros((nkpoint, nband))

# Fill in the data for each projection.
index = -1
for ions in ions_list:
    index += 2 
    for i in range(nkpoint):
        for j in range(nband):
            ions_contribup[ions][i, j] = float(data[8+i*(blocksize+3)+j*bandsize+index].strip().split()[10])

#-------------------For spin down-----------------------------
# Initialize matrices for each ion.
ions_contribdown = {}
for ions in ions_list:
    ions_contribdown[ions] = np.zeros((nkpoint, nband))

# Fill in the data for each projection.
index = -1
for ions in ions_list:
    index += 2
    for i in range(nkpoint):
        for j in range(nband):
            ions_contribdown[ions][i, j] = float(data[9+nkpoint*(blocksize+3)+i*(blocksize+3)+j*bandsize+index].strip().split()[10])


#******************Extract the k coordinate for each kpoints****************************
kcoord = []
for i in range(nkpoint):
    if (nions == 1):
        line = data[i*((nband*5)+3)+3].strip()
    else :
        line = data[i*(blocksize+3)+3].strip()
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
    # Calculate Euclidean distance between consecutive k-points
    kDis[i+1] = kDis[i] + np.sqrt((kcoord[i+1,0] - kcoord[i,0])**2 + (kcoord[i+1,1] - kcoord[i,1])**2 + (kcoord[i+1,2] - kcoord[i,2])**2)

# Define the special high-symmetry points as numpy arrays.
special_points = {
    'Γ': np.array([0.0, 0.0, 0.0]),
    'X': np.array([0.5, 0.0, 0.5]),
    'W': np.array([0.5, 0.25, 0.75]),
    'U': np.array([0.625, 0.25, 0.625]),
    'K': np.array([0.375, 0.375, 0.75]),
    'L': np.array([0.5, 0.5, 0.5])
}

# Enter here The desired order (even with repeated points) for high symmetry points in the 1st BZ:
path_order = ['Γ', 'X', 'W', 'K', 'Γ', 'L', 'U', 'W', 'L', 'K', 'U', 'X']
tol = 1e-3  # Tolerance for matching coordinates
tick_positions = [] 
tick_labels = []
current_search_index = 0

for sym in path_order:
    found_index = None
    # Search for the next occurrence of the special point starting at current_search_index
    for i in range(current_search_index, len(kcoord)):
        if np.linalg.norm(kcoord[i] - special_points[sym]) < tol:
            found_index = i
            break
    if found_index is not None:
        tick_positions.append(kDis[found_index])
        tick_labels.append(sym)
        current_search_index = found_index + 1  # Continue search from here
    else:
        print(f"Warning: Could not find high-symmetry point {sym} starting from index {current_search_index}.")


#******************************************************
#******* Create a figure with DOS and Band subplots ****​
#******************************************************
fig, (ax_dos, ax_band) = plt.subplots(1, 2, sharey=True,
                                        gridspec_kw={'width_ratios': [1, 2]},
                                        figsize=(12, 8))

#----- DOS Plot (left) -----
# Fill the area under the DOS curves with grey
ax_dos.fill_betweenx(Energy_dos - E_fermi, 0, DOSp, facecolor='grey', alpha=0.3)
ax_dos.fill_betweenx(Energy_dos - E_fermi, -DOSm, 0, facecolor='grey', alpha=0.3)
# Plot the DOS lines on top
ax_dos.plot(DOSp, Energy_dos - E_fermi, color='black', lw=1, label='DOS (spin up)')
ax_dos.plot(-DOSm, Energy_dos - E_fermi, color='red', ls="--", lw=1, label='DOS (spin down)')
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


#*******************************************************************************************************************
#***************************PLotting the band structure for the contribution of each projection**********************
#*******************************************************************************************************************
fig, (ax_band_projup, ax_band_projdown) = plt.subplots(1, 2, sharey=True,
                                        gridspec_kw={'width_ratios': [1, 1]},
                                        figsize=(12, 8))

# Define a list of colors to cycle through (extend as needed) without extra spaces.
color_list = ['gray', 'teal', 'cyan', 'skyblue', 'rosybrown', 'lightcoral', 'indianred', 'brown', 'firebrick']
colors = {}
for i, proj in enumerate(projections):
    colors[proj] = color_list[i % len(color_list)]

marker_scale = 100  # Adjust as needed for visual clarity
# Plot the projected band structure for each projection.
# Using different marker shapes for spin up (circle) and spin down (square).
for proj in projections:
    for i in range(nband):
        if i == 0:
            ax_band_projup.scatter(kDis, kbandenergyup[:, i] - E_fermi,
                                   s=proj_matricesup[proj][:, i] * marker_scale,
                                   color=colors[proj], marker='o', label=f"{proj}")
            # Uncomment the following block if you also want to plot spin down data.
            
            ax_band_projdown.scatter(kDis, kbandenergydown[:, i] - E_fermi,
                                   s=proj_matricesdown[proj][:, i] * marker_scale,
                                   color=colors[proj], marker='o', label=f"{proj}")
            
        else:
            ax_band_projup.scatter(kDis, kbandenergyup[:, i] - E_fermi,
                                   s=proj_matricesup[proj][:, i] * marker_scale,
                                   color=colors[proj], marker='o')
            
            ax_band_projdown.scatter(kDis, kbandenergydown[:, i] - E_fermi,
                                   s=proj_matricesdown[proj][:, i] * marker_scale,
                                   color=colors[proj], marker='o')
            

# Draw vertical dashed lines for high-symmetry points
for pos in tick_positions:
    ax_band_projup.axvline(x=pos, color='grey', linestyle='--', lw=0.8)
    ax_band_projdown.axvline(x=pos, color='grey', linestyle='--', lw=0.8)

#parameters for up
ax_band_projup.set_xticks(tick_positions)
ax_band_projup.set_xticklabels(tick_labels, fontsize=10)
ax_band_projup.axhline(y=0, color='red', lw=1.5, linestyle='--', label=r'$E_F$')
ax_band_projup.set_xlabel("k-path", fontsize=12)
ax_band_projup.set_title("Band Structure majority spin", fontsize=14)
ax_band_projup.set_xlim(0, kDis[-1])
ax_band_projup.grid(True, linestyle='--', alpha=0.5)
ax_band_projup.tick_params(axis='both', which='both', labelsize=10)
#parameters for down
ax_band_projdown.set_xticks(tick_positions)
ax_band_projdown.set_xticklabels(tick_labels, fontsize=10)
ax_band_projdown.axhline(y=0, color='red', lw=1.5, linestyle='--', label=r'$E_F$')
ax_band_projdown.set_xlabel("k-path", fontsize=12)
ax_band_projdown.set_title("Band Structure minority spin", fontsize=14)
ax_band_projdown.set_xlim(0, kDis[-1])
ax_band_projdown.grid(True, linestyle='--', alpha=0.5)
ax_band_projdown.tick_params(axis='both', which='both', labelsize=10)

# Create custom legend handles with a fixed markersize (e.g., 10)
legend_handlesup = []
legend_handlesdown = []
for proj in projections:
    # Handle for spin up (circle)
    up_handle = mlines.Line2D([], [], color=colors[proj], marker='o', linestyle='None',
                              markersize=10, label=f"{proj} up")
    # Handle for spin down (square)
    down_handle = mlines.Line2D([], [], color=colors[proj], marker='s', linestyle='None',
                                markersize=10, label=f"{proj} down")
    #legend_handles.extend([up_handle, down_handle])
    legend_handlesup.extend([up_handle])
    legend_handlesdown.extend([down_handle])

ax_band_projup.legend(handles=legend_handlesup, loc='upper right', fontsize=10)
ax_band_projdown.legend(handles=legend_handlesdown, loc='upper right', fontsize=10)

# Label the shared y-axis (energy)
ax_band_projup.set_ylabel("Energy (eV)", fontsize=12)
ax_band_projup.set_ylim(-10, 5)



#*******************************************************************************************
#****************************** Plotting each ions contributions ***************************
#*******************************************************************************************
fig, (ax_ionsdos, ax_ionsband) = plt.subplots(1, 2, sharey=True,
                                        gridspec_kw={'width_ratios': [1, 2]},
                                        figsize=(12, 8))
# Define a list of colors to cycle through (extend as needed) without extra spaces.
color_list = ['indianred', 'teal']
colors = {}
for i, ions in enumerate(ions_list):
    colors[ions] = color_list[i % len(color_list)]


#Plotting the DOS for each contribution
# ----- DOS Plot (left) -----
# Fill the area under the DOS curves with grey
ax_ionsdos.fill_betweenx(Energy_dos - E_fermi, 0, DOSp, facecolor='grey', alpha=0.3)
ax_ionsdos.fill_betweenx(Energy_dos - E_fermi, -DOSm, 0, facecolor='grey', alpha=0.3)
# Plot the DOS lines on top
ax_ionsdos.plot(DOSp, Energy_dos - E_fermi, color='black', lw=1, label='DOS (spin up)')
ax_ionsdos.plot(-DOSm, Energy_dos - E_fermi, color='black', ls="--", lw=1, label='DOS (spin down)')
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
                                   s=ions_contribup[ions][:, i] * marker_scale,
                                   color=colors[ions], marker='o', label=f"{ions}")
            # Uncomment the following block if you also want to plot spin down data.
            
            ax_ionsband.scatter(kDis, kbandenergydown[:, i] - E_fermi,
                                   s=ions_contribdown[ions][:, i] * marker_scale,
                                   color=colors[ions], marker='s', label=f"{ions}")
            
        else:
            ax_ionsband.scatter(kDis, kbandenergyup[:, i] - E_fermi,
                                   s=ions_contribup[ions][:, i] * marker_scale,
                                   color=colors[ions], marker='o')
            
            ax_ionsband.scatter(kDis, kbandenergydown[:, i] - E_fermi,
                                   s=ions_contribdown[ions][:, i] * marker_scale,
                                   color=colors[ions], marker='s')
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
    down_handle = mlines.Line2D([], [], color=colors[ions], marker='s', linestyle='None',
                                markersize=10, label=f"{ions} down")
    #legend_handles.extend([up_handle, down_handle])
    legend_handles.extend([up_handle, down_handle])

ax_ionsband.legend(handles=legend_handles, loc='upper right', fontsize=10)

# Label the shared y-axis (energy)
ax_ionsdos.set_ylabel("Energy (eV)", fontsize=12)
ax_ionsdos.set_ylim(-10, 5)


#********************************************************************************************
#***************************** PROJECTION contribution DOS and BAND **************************
#********************************************************************************************
fig, (ax_dos_proj, ax_band_proj) = plt.subplots(1, 2, sharey=True,
                                        gridspec_kw={'width_ratios': [1, 2]},
                                        figsize=(12, 8))

# Define a list of colors to cycle through (extend as needed) without extra spaces.
color_list = ['gray', 'teal', 'cyan', 'skyblue', 'rosybrown', 'lightcoral', 'indianred', 'brown', 'firebrick']
colors = {}
for i, proj in enumerate(projections):
    colors[proj] = color_list[i % len(color_list)]

marker_scale = 100  # Adjust as needed for visual clarity

#Plot the DOS left part
# Fill the area under the DOS curves with grey
ax_dos_proj.fill_betweenx(Energy_dos - E_fermi, 0, DOSp, facecolor='grey', alpha=0.3)
ax_dos_proj.fill_betweenx(Energy_dos - E_fermi, -DOSm, 0, facecolor='grey', alpha=0.3)
# Plot the DOS lines on top
ax_dos_proj.plot(DOSp, Energy_dos - E_fermi, color='black', lw=1, label='DOS (spin up)')
ax_dos_proj.plot(-DOSm, Energy_dos - E_fermi, color='black', ls="--", lw=1, label='DOS (spin down)')
ax_dos_proj.set_xlabel("DOS", fontsize=12)
ax_dos_proj.set_title("Density of States", fontsize=14)
ax_dos_proj.axhline(0, color='black', linestyle='--', lw=0.5)
ax_dos_proj.grid(True, linestyle='--', alpha=0.5)
ax_dos_proj.tick_params(axis='both', which='both', labelsize=10)
for proj in projections_list:
    ax_dos_proj.plot(dos_projup[proj], Energy_dos - E_fermi, color=colors[proj],label=f"{proj} spin maj")
    ax_dos_proj.plot(-dos_projdown[proj], Energy_dos - E_fermi, color=colors[proj],ls="--", label=f"{proj} spin min")
ax_dos_proj.legend(loc='upper left', fontsize=10)


#Plot the band right part
# Plot the projected band structure for each projection.
# Using different marker shapes for spin up (circle) and spin down (square).
for proj in projections:
    for i in range(nband):
        if i == 0:
            ax_band_proj.scatter(kDis, kbandenergyup[:, i] - E_fermi,
                                   s=proj_matricesup[proj][:, i] * marker_scale,
                                   color=colors[proj], marker='o', label=f"{proj}")
            # Uncomment the following block if you also want to plot spin down data.
            
            ax_band_proj.scatter(kDis, kbandenergydown[:, i] - E_fermi,
                                   s=proj_matricesdown[proj][:, i] * marker_scale,
                                   color=colors[proj], marker='s', label=f"{proj}")
            
        else:
            ax_band_proj.scatter(kDis, kbandenergyup[:, i] - E_fermi,
                                   s=proj_matricesup[proj][:, i] * marker_scale,
                                   color=colors[proj], marker='o')
            
            ax_band_proj.scatter(kDis, kbandenergydown[:, i] - E_fermi,
                                   s=proj_matricesdown[proj][:, i] * marker_scale,
                                   color=colors[proj], marker='s')
            

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
for proj in projections:
    # Handle for spin up (circle)
    up_handle = mlines.Line2D([], [], color=colors[proj], marker='o', linestyle='None',
                              markersize=10, label=f"{proj} up")
    # Handle for spin down (square)
    down_handle = mlines.Line2D([], [], color=colors[proj], marker='s', linestyle='None',
                                markersize=10, label=f"{proj} down")
    legend_handles.extend([up_handle, down_handle])
    #legend_handlesup.extend([up_handle, down_handle])

ax_band_proj.legend(handles=legend_handles, loc='upper right', fontsize=10)

# Label the shared y-axis (energy)
ax_dos_proj.set_ylabel("Energy (eV)", fontsize=12)
ax_dos_proj.set_ylim(-10, 5)


plt.tight_layout()
plt.show()
plt.close()

