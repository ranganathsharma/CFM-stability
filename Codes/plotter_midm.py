# This file has the code to merge the results of stability analysis for different densities for a given model

# Imports

import os, sys, yaml, glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as patches


# Merge the results.csv into a single csv file



# Define the pattern for matching CSV files (modify 'model' accordingly)
model = "midm"  # Change this to the actual model name
file_path = os.path.dirname(os.path.abspath(__file__))
csv_files = glob.glob(os.path.join(file_path, f"{model}_*_stability_status.csv"))  # Finds all matching files

# List to store individual DataFrames
dfs = []


file_path = os.path.join(file_path, f'{model}_stability_status.csv')

if not os.path.exists(file_path):
    dfs = []  # Ensure dfs list is initialized

    # Read each CSV file and append it to the list
    for file in csv_files:
        df = pd.read_csv(file, index_col=False)  # Avoid index shift
        dfs.append(df)

    # Merge all DataFrames
    merged_df = pd.concat(dfs, ignore_index=True)  # Reset index

    # Drop extra index column if it appears
    if "Unnamed: 0" in merged_df.columns:
        merged_df = merged_df.drop(columns=["Unnamed: 0"])

    merged_df.to_csv(file_path, index=False)
    data = merged_df
    del merged_df, dfs, df

    for file in csv_files:
        os.remove(file)

    print('The merged file has been created')
else:
    data = pd.read_csv(file_path, index_col=False)  # Avoid index issues
    data = data[(data['dens']>6) & (data['dens'] <142)]
    print('The data is loaded')

# Plotting the types of instabilities observed at different densities and perturbations depending on the type of instability

def data_correction():
    # df['column_name'] = df['column_name'].fillna(0)
    data['convective_analytical'] = data['convective_analytical'].fillna(0)
    data.to_csv(f'{model}_stability_status.csv')

# data_correction()

def patching(ax, points, color):
    x = np.linspace(0.1, 4, 10)
    y = np.linspace(0.1, 4, 10)
    dx = x[1] - x[0]
    dy = y[1] - x[0]

    for count, coordinates in enumerate(points):
        vertices = [(x[coordinates[0]-1] - dx/2, y[coordinates[1]-1] - dy/2),
                    (x[coordinates[0]-1] + dx/2, y[coordinates[1]-1] - dy/2),
                    (x[coordinates[0]-1] + dx/2, y[coordinates[1]-1] + dy/2),
                    (x[coordinates[0]-1] - dx/2, y[coordinates[1]-1] + dy/2)]
        polygon = patches.Polygon(vertices, closed = True, linewidth = 2, edgecolor = 'none', facecolor = color, alpha = 0.1)
        ax.add_patch(polygon)

    return ax


# Plotting the phase plots for varying perturbation strength and density conditions
def phase_plot(data: pd.DataFrame,
                instability: str,
                density: float,
                da: float,
                dtau: float):
    # This function plots different instabilities over the parameter phase space
    
    dv = da*dtau
    data['dv'] = data['da']*data['dtau']
    data = data[data['dens'] == density]
    data = data[(data['dens'] == density) & (data['dv'] == dv)]   
    
    fig, ax = plt.subplots(figsize=(3, 3), dpi= 600)

    for i in range(len(data)):
        temp_data = data.iloc[i]
        if temp_data[instability] == 1:
            ax.scatter(temp_data['max_acc'], temp_data['safe_head'],
                        marker = 's', alpha = 1, s = 20, color = 'b')
        elif temp_data[instability] == 2:
            ax.scatter(temp_data['max_acc'], temp_data['safe_head'],
                        marker = 's', alpha = 1, s = 20, color = 'r')
        elif temp_data[instability] == 3:
            ax.scatter(temp_data['max_acc'], temp_data['safe_head'],
                        marker = 's', alpha = 1, s = 20, color = 'g')
        elif temp_data[instability] == 0:
            ax.scatter(temp_data['max_acc'], temp_data['safe_head'],
                        marker = 'o', alpha = 1, s = 20, color = 'm')
        else:
            ax.scatter(temp_data['max_acc'], temp_data['safe_head'],
                       marker = 'o', alpha  = 1, s = 20, color = 'y')
            
    if np.any(data[instability] == 1):    
        ax.scatter(-1, 0, marker = 's', alpha = 1, s = 20, color = 'b', label = 'Downstream')

    if np.any(data[instability] == 2):
        ax.scatter(-1, 0, marker = 's', alpha = 1, s = 20, color = 'r', label = 'Absolute')

    if np.any(data[instability] == 3):
        ax.scatter(-1, 0, marker = 's', alpha = 1, s = 20, color = 'g', label = 'Upstream')
    
    if np.any(data[instability] == 0):
        ax.scatter(-1, 0, marker = 'o', alpha = 1, s = 20, color = 'm', label = 'Stable')
    
    if np.any(data[instability] == 0):
        ax.scatter(-1, 0, marker = 'o', alpha = 1, s = 20, color = 'y', label = 'Accident')

    handles, labels = plt.gca().get_legend_handles_labels()

    if len(labels) <= 2:
        ax.scatter(-1, 0, marker = 'o', alpha = 0, s = 20, color = 'm', label = ' ')
    else: 
        pass

    if density == 22:

        # comparison with analytical

        # coordinates = [(2, 5), (2, 6), (2, 7), (2, 8), (2, 9) ]
        # ax = patching(ax, coordinates, 'g')


        # comparison with large perturbation

        coordinates = [(1, 1), (2, 1), (1, 7), (1, 8), (1, 9), (1, 10), (2, 9)]
        ax = patching(ax, coordinates, 'b')


        pass

    if density == 130:

        coordinates = [(1, 5)]
        ax = patching(ax, coordinates, 'r')

        coordinates = [(4, 1), (7, 2)]
        ax = patching(ax, coordinates, 'g')

        # comparison with large perturbation



        pass
    

    plt.xlabel('$a_{max}$', fontsize = 15)
    plt.ylabel('T', fontsize = 15)
    plt.tick_params(axis = 'both', which = 'major', labelsize = 12)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
    plt.xlim(-0.1, 4.1)
    plt.ylim(-0.1, 4.1)
    plt.savefig(f'{model}_phase_space_{instability}_{density}_{da}_{dtau}.png', bbox_inches = 'tight')
    # plt.show()
    plt.close()
    
# phase_plot(data= data, instability = 'convective_analytical', density =  130, da = 0.1, dtau = 1)
# phase_plot(data= data, instability = 'convective_numerical', density =  22 , da = 5, dtau = 1)

def error_proportions(data):

    data = data[(data['dens']>6) & (data['dens'] <142)]
    data['dv'] = data['da'] * data['dtau']

    dens_array = np.unique(data['dens'])
    dv_list = np.unique(data['dv'])
    dv_list = dv_list[dv_list < 10]
    # dv_list = [0.1]

    fig, ax = plt.subplots(figsize=(5, 4), dpi = 600)

    colors = ['r', 'g']  # Colors for the cases
    labels = ["Case 1", "Case 2"]  # Labels for the legend

    for count, dens in enumerate(dens_array):
        temp_data = data[data['dens'] == dens]

        cases = [
            # Merging Case 1 and Case 2 with OR (|)
            temp_data.loc[
                ((temp_data['convective_analytical'] != 0) & (temp_data['convective_numerical'] == 0)) | 
                ((temp_data['convective_numerical'].isin([1, 3])) & (temp_data['convective_analytical'] == 2))
            ],  
            # Case 3 remains unchanged
            temp_data.loc[
                ((temp_data['convective_analytical'].isin([1, 3])) & (temp_data['convective_numerical'] == 2) |
                ((temp_data['convective_analytical'] == 0) & (temp_data['convective_numerical'] != 0)))
            ]]

        for case_idx, stability_data in enumerate(cases):
            counts = [len(stability_data[stability_data['dv'] == dv]) / 100 for dv in dv_list]

            if counts:
                mean_val = np.mean(counts)
                min_val = np.min(counts)
                max_val = np.max(counts)

                color = colors[case_idx]
                label = labels[case_idx]

                # Plot min-max range as a bar
                ax.vlines(dens, min_val, max_val, color=color, linewidth=4, alpha=0.5)

                # Scatter for mean value with border
                ax.scatter(dens, mean_val, facecolors=color, label=label if count == 0 else "", s=30)

    ax.set_xlabel('$\\rho_e$', fontsize=15)
    ax.set_ylabel('Proportion', fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=12)
    ax.set_ylim(-0.05, 0.5)
    
    if len(dv_list) == 1:
        plt.savefig(f'{model}_error_proportions_small_pert.png', bbox_inches='tight')
    else:
        plt.savefig(f'{model}_error_proportions.png', bbox_inches='tight')

    plt.close()

# error_proportions(data)

def nonlinear(data):

    data['dv'] = data['da']*data['dtau']

    reference_data = data[data['dv'] == 0.1]

    dens_list = np.unique(data['dens'])
    dv_list = np.unique(data['dv'])
    dv_list = dv_list[dv_list < 20]

    NL1 = np.zeros((len(dv_list), len(dens_list)))
    NL2 = np.zeros_like(NL1)

    for i, dv in enumerate(dv_list):
        for j, dens in enumerate(dens_list):

            target_data = data[(data['dens'] == dens) & (data['dv'] == dv)]
            temper_data = reference_data[reference_data['dens'] == dens]

            if len(target_data) != len(temper_data):
                print('Something wrong', dens, dv, len(target_data), len(temper_data))

            else:
                
                merged = temper_data.merge(target_data[['safe_head', 'max_acc', 'convective_numerical']],  # Only the column of interest from df2
                                    on=['safe_head', 'max_acc'],
                                    how='inner',
                                    suffixes=('_ref', '_tar'))

                NL1_count = ((merged['convective_numerical_ref'] == 0) & (merged['convective_numerical_tar'] != 0)).sum()
                NL2_count = ((merged['convective_numerical_ref'].isin([1, 3])) & (merged['convective_numerical_tar'] == 2)).sum()

                NL1[i, j] = NL1_count/100
                NL2[i, j] = NL2_count/100

    DENS, DV = np.meshgrid(dens_list, dv_list)

    fig = plt.figure(figsize=(7, 6), dpi=300)

    # NL1 surface
    ax1 = fig.add_subplot(121, projection='3d')
    surf1 = ax1.plot_surface(DENS, DV, NL1, cmap='Blues', edgecolor='k', linewidth=0.5, antialiased=True)
    ax1.set_xlabel('$\\rho_e$', labelpad=15)
    ax1.set_ylabel('$\\delta v$', labelpad=15)
    ax1.set_zlabel('Proportion', labelpad=15)
    ax1.set_title('Nonlinear Instability: Type 1')
    ax1.view_init(elev=30, azim=-80)
    ax1.tick_params(axis = 'both', which = 'major', labelsize = 10)
    fig.colorbar(surf1, ax=ax1, orientation='horizontal', shrink=0.6, pad=0.1)

    # NL2 surface
    ax2 = fig.add_subplot(122, projection='3d')
    surf2 = ax2.plot_surface(DENS, DV, NL2, cmap='Reds', edgecolor='k', linewidth=0.5, antialiased=True)
    ax2.set_xlabel('$\\rho_e$', labelpad=15)
    ax2.set_title('Nonlinear Instability: Type 2')
    ax2.view_init(elev=30, azim=-80)
    fig.colorbar(surf2, ax=ax2, orientation='horizontal', shrink=0.6, pad=0.1)

    # Save
    fig.subplots_adjust(left=0.05, right=0.95, top=0.92, bottom=0.08)
    plt.savefig(f'{model}_nonlinear_instability_surface.png', bbox_inches='tight')
    plt.close()       
        
nonlinear(data)

# Plotting the proportion of different instabilities 
def proportions(data, quantity):

    data = data[(data['dens']>6) & (data['dens'] <142)]

    dens_array = np.unique(data['dens'])
    data['dv'] = data['da'] * data['dtau']
    data = data[data['dv'] < 10]

    dv_list = np.unique(data['dv'])

    stability_cases = {
        0: ('stable', 'm'),
        1: ('downstream', 'b'),
        2: ('absolute', 'r'),
        3: ('upstream', 'g')
    }

    fig, ax = plt.subplots(figsize=(5, 4))

    for count, dens in enumerate(dens_array):
        temp_data = data[data['dens'] == dens]

        for case, (label, color) in stability_cases.items():
            stability_data = temp_data[temp_data[quantity] == case]

            counts = [len(stability_data[stability_data['dv'] == dv]) / 100 for dv in dv_list]

            mean_val = np.mean(counts)
            min_val = np.min(counts)
            max_val = np.max(counts)

            # Plot min-max range as a bar
            ax.vlines(dens, min_val, max_val, color=color, linewidth=4, alpha=0.5)

            # Scatter for mean value
            ax.scatter(dens, mean_val, color=color, label=label if count == 0 else "", s=30)

    ax.set_xlabel('$\\rho_e$', fontsize=15)
    ax.set_ylabel('Proportion', fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=12)
    
    plt.savefig(f'{model}_proportions_{quantity}.png', bbox_inches='tight')
    plt.close()

# proportions(data, 'convective_numerical')
# proportions(data, 'convective_analytical')  

def parameter_perturbation(data: pd.DataFrame, density: float, instability: str):
    data = data[data['dens'] == density].copy()  # Filter by density
    
    data['dv'] = data['da'] * data['dtau']  # Compute dv

    # Extract unique da and dtau values
    unique_da = sorted(data['da'].unique())  # 3 unique da values (rows)
    unique_dtau = sorted(data['dtau'].unique())  # 3 unique dtau values (columns)

    fig, axes = plt.subplots(3, 3, figsize=(8,8), dpi=600, sharex=True, sharey=True)  # 3x3 grid
    colors = {1: 'b', 2: 'r', 3: 'g', 0: 'm'}  # Color mapping for instabilities
    labels = {1: 'Downstream unstable', 2: 'Absolute unstable', 3: 'Upstream unstable', 0: 'Stable'}  # Legends
    markers = {1: 's', 2: 's', 3: 's', 0: 'o'}

    scatter_plots = []  # Store scatter plots for legend

    for row, da_value in enumerate(unique_da):  # Loop over da (rows)
        for col, dtau_value in enumerate(unique_dtau):  # Loop over dtau (columns)
            ax = axes[row, col]  # Select correct subplot

            # Filter data for the given da and dtau
            subset = data[(data['da'] == da_value) & (data['dtau'] == dtau_value)]

            for instability_type, color in colors.items():
                sub_subset = subset[subset[instability] == instability_type]
                scatter = ax.scatter(sub_subset['max_acc'], sub_subset['safe_head'], color=color, alpha=0.5, marker=markers[instability_type])
                if row == 0 and col == 0:  # Only add to legend once
                    scatter_plots.append(scatter)

            ax.set_title(f"$\delta a$ = {da_value:.2f}, $\\delta \\tau$ = {dtau_value:.2f}", fontsize=10)
            ax.set_xlim(-0.1, 4.1)
            ax.set_ylim(-0.1, 4.1)

            # Set x-labels only for bottom row
            if row != 2:  
                ax.tick_params(axis='x', labelbottom=False)
            else:
                ax.set_xlabel("$a_{max}$", fontsize=15)

            # Set y-labels only for leftmost column
            if col != 0:  
                ax.tick_params(axis='y', labelleft=False)
            else:
                ax.set_ylabel("T", fontsize=15)

    # Adjust layout to make space for the legend
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)  # Shift plots up to make room for legend

    # Add a single legend at the bottom
    fig.legend(scatter_plots, labels.values(), loc='lower center', ncol=3, fontsize=12)
    plt.savefig(f'{model}_param_perturb_{density}.png')
    plt.close()
# parameter_perturbation(data = data, density = 22, instability = 'convective_numerical')

def density_perturbation_space(data: pd.DataFrame, instability: str, max_acc: float, safe_head: float):
    max_acc = np.round(max_acc, 3)
    safe_head = np.round(safe_head, 3)

    data[['max_acc', 'safe_head']] = data[['max_acc', 'safe_head']].round(3)
    
    data = data[(data['max_acc'] == max_acc) & (data['safe_head'] == safe_head)].copy()
    data['dv'] = data['da'] * data['dtau']
    unique_dv = sorted(data['dv'].unique())  # Get unique dv values
    dv_mapping = {dv: i+1 for i, dv in enumerate(unique_dv)}  # Map dv values to evenly spaced integers

    plt.figure(figsize=(4, 3), dpi = 100)

    for i in range(len(data)):
        dv_mapped = dv_mapping[data.iloc[i]['dv']]  # Use mapped values for even spacing
        if data.iloc[i][instability] == 1:
            plt.scatter(data.iloc[i]['dens'], dv_mapped, marker='s', color='b', alpha=0.5, s = 20)
        elif data.iloc[i][instability] == 2:
            plt.scatter(data.iloc[i]['dens'], dv_mapped, marker='s', color='r', alpha=0.5, s = 20)
        elif data.iloc[i][instability] == 3:
            plt.scatter(data.iloc[i]['dens'], dv_mapped, marker='s', color='g', alpha=0.5, s = 20)
        elif data.iloc[i][instability] == 0:
            plt.scatter(data.iloc[i]['dens'], dv_mapped, marker='o', color='m', alpha=0.5, s = 20)

    # Legend for instability types
    legend_handles = [
        plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='g', markersize=6.5, alpha = 0.5, label='Upstream instability'),
        plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='b', markersize=6.5, alpha = 0.5, label='Downstream instability'),
        plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='r', markersize=6.5, alpha = 0.5, label='Absolute instability'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='m', markersize=6.5, alpha = 0.5, label='Stability')
    ]

    # Set labels and format y-ticks
    plt.xlabel('$\\rho_e$', fontsize=15)
    plt.ylabel('$\\delta v$', fontsize=15)
    plt.xticks(fontsize=12)
    plt.xlim(0, 150)
    plt.yticks(ticks=list(dv_mapping.values()), labels=[f"{dv:.2f}" for dv in unique_dv], fontsize=12)  # Replace y-axis ticks with dv values
    plt.legend(handles=legend_handles, loc='lower center', bbox_to_anchor=(0.5, -0.5), ncol=2, fontsize=10)
    # plt.savefig(f'density_perturbation_{max_acc}_{safe_head}.png', bbox_inches = 'tight')
    plt.show()
    plt.close()

# max_acc_array = safe_head_array = np.linspace(0.1, 4, 10)
# for max_acc in max_acc_array:
#     for safe_head in safe_head_array:

#         density_perturbation_space(data = data, instability = 'convective_analytical', max_acc = max_acc, safe_head = safe_head)

# density_perturbation_space(data = data, instability = 'convective_analytical', max_acc = 0.533, safe_head = 4)

def stability_class(quantity):
    data = pd.read_csv(f'stability_diagram_{model}.csv')
    column_name = f'stability_{quantity}'

    print(np.unique(data[column_name]))
    fig, ax = plt.subplots(figsize=(3, 3), dpi=600)
    for i in range(len(data)):
        if data.iloc[i][column_name] == 's':
            ax.scatter(data.iloc[i]['max_acc'], data.iloc[i]['safe_head'], marker = 'p', alpha = 1, s = 20, color = 'g')
        elif data.iloc[i][column_name] in ['sna', 'snau', 'snaua', 'snua', 'sa', 'sau', 'sua']:
            ax.scatter(data.iloc[i]['max_acc'], data.iloc[i]['safe_head'], marker = '^', alpha = 1, s = 20, color = 'r')
        elif data.iloc[i][column_name] in ['snu', 'su', 'snu ']:
            ax.scatter(data.iloc[i]['max_acc'], data.iloc[i]['safe_head'], marker = 'v', alpha = 1, s = 20, color = 'b')
        elif data.iloc[i][column_name] in ['snus', 'sds', 'sns', 'snds', 'sus', 'sdns']:
            ax.scatter(data.iloc[i]['max_acc'], data.iloc[i]['safe_head'], marker = 'v', alpha = 1, s = 20, color = 'c')
        elif data.iloc[i][column_name] in ['snaus', 'sdaus', 'sndaus', 'sdaus', 'snas', 'sdnaus', 'sndaus', 'sndas', 'sas']:
            ax.scatter(data.iloc[i]['max_acc'], data.iloc[i]['safe_head'], marker = '^', alpha = 1, s = 20, color = 'y')

    # analytical

    # if quantity == 'analytical':

    #     coordinates = [(8, 1), (9, 1), (10, 1), (10, 2), (3, 2), (4, 3), (5, 3), (3, 3), (3, 4), (3, 5), (2, 5), (2, 6), (2, 7), (2, 8), (2, 9), (2, 10), (1, 9), (1, 10)]
    #     ax = patching(ax, coordinates, 'b')

    # elif quantity == 'numerical':

    #     coordinates = [(1, 10), (1, 9), (2, 10), (2, 9), (2, 8), (2, 4), (2, 5), (3, 7), (3, 6), (3, 5), (3, 4), (4, 4), (4, 3), (5, 3), (9, 2), (10, 2)]
    #     ax = patching(ax, coordinates, 'b')


    plt.scatter(-1, -1, marker = 'v', s = 20, color = 'c', alpha = 1, label = '2b')
    plt.scatter(-1, -1, marker = 'v', s = 20, color = 'b', alpha = 1, label = '2a')
    plt.scatter(-1, -1, marker = '^', s = 20, color = 'r', alpha = 1, label = '1a')
    plt.scatter(-1, -1, marker = '^', s = 20, color = 'y', alpha = 1, label = '1b')
    plt.scatter(-1, -1, marker = 'p', s = 20, color = 'g', alpha = 1, label = '3')  

    plt.xlabel('$a_{max}$', fontsize = 15)
    plt.ylabel('T', fontsize = 15)
    plt.tick_params(axis = 'both', which = 'major', labelsize = 12)
    plt.xlim(-0.1, 4.1)
    plt.ylim(-0.1, 4.1)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=12)
    plt.savefig(f'{model}_stability_classes_{quantity}.png', bbox_inches = 'tight')
    plt.close()

# stability_class('numerical')
# stability_class('analytical')

def proportions_compare(quantity):
    model = 'idm'
    data = pd.read_csv(f'{model}_stability_status.csv')
    data = data[data['dens'] != 2]
    data = data[data['dens'] != 6]

    dens_array = np.unique(data['dens'])
    data['dv'] = data['da']*data['dtau']

    
    dv_list = np.unique(data['dv'])

    upstream_mean = np.zeros_like(dens_array)
    upstream_error = np.zeros_like(dens_array)

    downstream_mean = np.zeros_like(dens_array)
    downstream_error = np.zeros_like(dens_array)

    absolute_mean = np.zeros_like(dens_array)
    absolute_error = np.zeros_like(dens_array)

    stable_mean = np.zeros_like(dens_array)
    stable_error = np.zeros_like(dens_array)

    stability_cases = {0: 'stable',
                       1: 'downstream',
                       2: 'absolute',
                       3: 'upstream'}

    for count, dens in enumerate(dens_array):
        temp_data = data[data['dens']== dens]

        for case, value in stability_cases.items():
            stability_data = temp_data[temp_data[quantity] == case]
            
            counts = []
            for dv in dv_list:
                counts.append(len(stability_data[stability_data['dv'] == dv])/100)
            locals()[f"{value}_mean"][count] = np.mean(counts)
            locals()[f"{value}_error"][count] = np.std(counts)

    plt.errorbar(dens_array, upstream_mean, yerr = upstream_error, label = 'upstream instability', linewidth = 0.8, color = 'g', fmt = '-o', capsize = 2, markersize = 3)
    plt.errorbar(dens_array, downstream_mean, yerr = downstream_error, label = 'downstream instability', linewidth = 0.8, color = 'b', fmt = '-o', capsize = 2, markersize = 3)
    plt.errorbar(dens_array, absolute_mean, yerr = absolute_error, label = 'absolute instability - idm', linewidth = 0.8, color = 'r', fmt = '-o', capsize = 2, markersize = 3)
    plt.errorbar(dens_array, stable_mean, yerr = stable_error, label = 'stability - idm', linewidth = 0.8, color = 'm', fmt = '-o', capsize = 2, markersize = 3)
    
    model = 'midm'
    data = pd.read_csv(f'{model}_stability_status.csv')

    data = data[data['dens'] != 2]
    data = data[data['dens'] != 6]

    dens_array = np.unique(data['dens'])
    data['dv'] = data['da']*data['dtau']

    
    dv_list = np.unique(data['dv'])

    upstream_mean = np.zeros_like(dens_array)
    upstream_error = np.zeros_like(dens_array)

    downstream_mean = np.zeros_like(dens_array)
    downstream_error = np.zeros_like(dens_array)

    absolute_mean = np.zeros_like(dens_array)
    absolute_error = np.zeros_like(dens_array)

    stable_mean = np.zeros_like(dens_array)
    stable_error = np.zeros_like(dens_array)

    stability_cases = {0: 'stable',
                       1: 'downstream',
                       2: 'absolute',
                       3: 'upstream'}

    for count, dens in enumerate(dens_array):
        temp_data = data[data['dens']== dens]

        for case, value in stability_cases.items():
            stability_data = temp_data[temp_data[quantity] == case]
            
            counts = []
            for dv in dv_list:
                counts.append(len(stability_data[stability_data['dv'] == dv])/100)
            locals()[f"{value}_mean"][count] = np.mean(counts)
            locals()[f"{value}_error"][count] = np.std(counts)

    plt.errorbar(dens_array, upstream_mean, yerr = upstream_error, label = 'upstream instability', alpha = 0.5, linewidth = 0.8, color = 'g', fmt = '--', capsize = 2, markersize = 3)
    plt.errorbar(dens_array, downstream_mean, yerr = downstream_error, label = 'downstream instability',alpha = 0.5, linewidth = 0.8, color = 'b', fmt = '--', capsize = 2, markersize = 3)
    plt.errorbar(dens_array, absolute_mean, yerr = absolute_error, label = 'absolute instability - midm',alpha = 0.5, linewidth = 0.8, color = 'r', fmt = '--', capsize = 2, markersize = 3)
    plt.errorbar(dens_array, stable_mean, yerr = stable_error, label = 'stability - midm',alpha = 0.5, linewidth = 0.8, color = 'm', fmt = '--', capsize = 2, markersize = 3)
    
    
    plt.xlabel('$\\rho_e$', fontsize = 15)
    plt.ylabel('Proportion', fontsize = 15)
    plt.tick_params(axis = 'both', which = 'major', labelsize = 12)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=12)
    plt.savefig(f'{model}_proportions_compare_{quantity}.png', bbox_inches = 'tight')
    plt.close()

# proportions_compare('convective_numerical')




