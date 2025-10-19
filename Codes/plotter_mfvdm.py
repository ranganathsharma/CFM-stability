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
model = "mfvdm"  # Change this to the actual model name
file_path = os.path.dirname(os.path.abspath(__file__))
csv_files = glob.glob(os.path.join(file_path, f"{model}_*_stability_status_analytical.csv"))  # Finds all matching files

# List to store individual DataFrames
dfs = []


# file_path = os.path.join(file_path, f'{model}_stability_status_analytical.csv') # changed
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

    # Save the merged DataFrame without an index
    # merged_df = merged_df[merged_df['string_numerical2'].isin([0, 1])]
    # Save the merged DataFrame without an index
    merged_df.to_csv(file_path, index=False)
    data = merged_df
    del merged_df, dfs, df

    for file in csv_files:
        os.remove(file)

    data = data[(data['dens']>6) & (data['dens'] <142)]
    print('The merged file has been created, and intermediary files deleted')

else:
    data = pd.read_csv(file_path, index_col=False)  # Avoid index issues
    data = data[(data['dens']>6) & (data['dens'] <142)]
    print('The data is loaded')

def density_perturbation_space(data: pd.DataFrame, instability: str, sens_relvel: float, sens_vel: float):
    sens_relvel = round(sens_relvel, 3)
    sens_vel = round(sens_vel, 3)
    data[['sens_relvel', 'sens_vel']] = data[['sens_relvel', 'sens_vel']].round(3)
    
    data = data[(data['sens_relvel'] == sens_relvel) & (data['sens_vel'] == sens_vel)].copy()
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
    plt.xlim(-1, 201)
    plt.yticks(ticks=list(dv_mapping.values()), labels=[f"{dv:.2f}" for dv in unique_dv], fontsize=12)  # Replace y-axis ticks with dv values
    plt.legend(handles=legend_handles, loc='lower center', bbox_to_anchor=(0.5, -0.5), ncol=2, fontsize=10)
    # plt.savefig(f'density_perturbation_{sens_relvel}_{sens_vel}.png', bbox_inches = 'tight')
    plt.show()
    plt.close()

# sens_vel_array = np.linspace(0.2, 1, 10)
# sens_relvel_array = np.linspace(0.2, 1,10)
# for sens_vel in sens_vel_array:
#     for sens_relvel in sens_relvel_array:
#         density_perturbation_space(data = data, instability = 'convective_analytical', sens_relvel = sens_relvel, sens_vel = sens_vel)

# density_perturbation_space(data = data, instability = 'convective_numerical', sens_relvel = 0.467, sens_vel = 5.556)

# Plotting the phase plots for varying perturbation strength and density conditions

def patching(ax, points, color):
    x = np.linspace(0.2, 1, 10)
    y = np.linspace(0.2, 1, 10)
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

def phase_plot(data: pd.DataFrame,
                instability: str,
                density: float,
                da: float,
                dtau: float):
    # This function plots different instabilities over the parameter phase space
    
    dv = da*dtau
    data['dv'] = data['da']*data['dtau']
    data = data[data['dens'] == density]
    # print(data.head())
    data = data[(data['dens'] == density) & (data['dv'] == dv)]
    
    fig, ax = plt.subplots(figsize=(3, 3), dpi=600)

    for i in range(len(data)):
        temp_data = data.iloc[i]
        if temp_data[instability] == 1:
            ax.scatter(temp_data['sens_relvel'], temp_data['sens_vel'],
                        marker = 's', alpha = 1, s = 20, color = 'b')
        elif temp_data[instability] == 2:
            ax.scatter(temp_data['sens_relvel'], temp_data['sens_vel'],
                        marker = 's', alpha = 1, s = 20, color = 'r')
        elif temp_data[instability] == 3:
            ax.scatter(temp_data['sens_relvel'], temp_data['sens_vel'],
                        marker = 's', alpha = 1, s = 20, color = 'g')
        elif temp_data[instability] == 0:
            ax.scatter(temp_data['sens_relvel'], temp_data['sens_vel'],
                        marker = 'o', alpha = 1, s = 20, color = 'm')
        else:
            ax.scatter(temp_data['sens_relvel'], temp_data['sens_vel'],
                        marker = 'o', alpha = 1, s = 20, color = 'y')
            
    if np.any(data[instability] == 1):    
        ax.scatter(-1, 0, marker = 's', alpha = 1, s = 20, color = 'b', label = 'Downstream')

    if np.any(data[instability] == 2):
        ax.scatter(-1, 0, marker = 's', alpha = 1, s = 20, color = 'r', label = 'Absolute')

    if np.any(data[instability] == 3):
        ax.scatter(-1, 0, marker = 's', alpha = 1, s = 20, color = 'g', label = 'Upstream')
    
    if np.any(data[instability] == 0):
        ax.scatter(-1, 0, marker = 'o', alpha = 1, s = 20, color = 'm', label = 'Stable')
    
    if np.any(np.isnan(data[instability])):
        ax.scatter(-1, 0, marker = 'o', alpha = 1, s = 20, color = 'y', label = 'Accident')

    handles, labels = plt.gca().get_legend_handles_labels()

    if len(labels) <= 2:
        ax.scatter(-1, 0, marker = 'o', alpha = 0, s = 20, color = 'm', label = ' ')
    else: 
        pass
    

    if density == 70:

        # comparison with analytical results

        # coordinates = [(2, 1), (2, 2), (1, 3), (4, 3), (3, 4), (1, 7)]
        # ax = patching(ax, coordinates, 'r')

        # comparison with larger perturbation

        coordinates = [(1, 2)]
        ax = patching(ax, coordinates, 'b')

        pass
    
    plt.xlabel('$\\kappa_2$', fontsize = 15)
    plt.ylabel('$\\kappa_1$', fontsize = 15)
    plt.tick_params(axis = 'both', which = 'major', labelsize = 12)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
    plt.xlim(0.15, 1.05)
    plt.ylim(0.15, 1.05)
    plt.savefig(f'{model}_phase_space_{instability}_{density}_{da}_{dtau}.png', bbox_inches = 'tight')
    # plt.show()
    plt.close()
    
# phase_plot(data= data, instability = 'convective_analytical', density =  70, da = 0.1, dtau = 1)
# phase_plot(data= data, instability = 'convective_numerical', density =  70, da = 0.1, dtau = 1)

def parameter_perturbation(data: pd.DataFrame, density: float, instability: str):

    data = data.loc[np.isclose(data['dens'], density)]

    data['dv'] = data['da'] * data['dtau']
    

    unique_da = sorted(data['da'].unique())
    unique_dtau = sorted(data['dtau'].unique())

    fig, axes = plt.subplots(3, 3, figsize=(10, 8), dpi=600, sharex=True, sharey=True)

    # Define categories with desired order
    categories = [0, 1, 2, 3, 'accident']
    colors = {0: 'm', 1: 'b', 2: 'r', 3: 'g', 'accident': 'y'}
    labels = {0: 'Stable', 1: 'Downstream unstable', 2: 'Absolute unstable', 3: 'Upstream unstable', 'accident': 'Accident'}
    markers = {0: 'o', 1: 's', 2: 's', 3: 's', 'accident': 'X'}

    scatter_handles = {}

    for row, da_value in enumerate(unique_da):
        for col, dtau_value in enumerate(unique_dtau):
            ax = axes[row, col]
            subset = data[(data['da'] == da_value) & (data['dtau'] == dtau_value)]

            for cat in categories:
                if cat == 'accident':
                    sub_subset = subset[subset[instability].isna()]
                else:
                    sub_subset = subset[subset[instability] == cat]

                if not sub_subset.empty:
                    scatter = ax.scatter(
                        sub_subset['sens_relvel'], sub_subset['sens_vel'],
                        color=colors[cat], alpha=0.5, marker=markers[cat], label=labels[cat]
                    )
                    # Store only one handle per category
                    if cat not in scatter_handles:
                        scatter_handles[cat] = scatter

            ax.set_title(f"$\delta a$ = {da_value:.2f}, $\\delta \\tau$ = {dtau_value:.2f}", fontsize=10)
            ax.set_xlim(0.15, 1.05)
            ax.set_ylim(0.15, 1.05)

            if row != 2:
                ax.tick_params(axis='x', labelbottom=False)
            else:
                ax.set_xlabel("$\\kappa_2$", fontsize=15)

            if col != 0:
                ax.tick_params(axis='y', labelleft=False)
            else:
                ax.set_ylabel("$\kappa_1$", fontsize=15)

    # Build legend from handles in correct order
    ordered_handles = [scatter_handles[cat] for cat in categories if cat in scatter_handles]
    ordered_labels = [labels[cat] for cat in categories if cat in scatter_handles]

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)
    fig.legend(ordered_handles, ordered_labels, loc='lower center', ncol=3, fontsize=12)
    plt.savefig(f'{model}_param_perturb_{density}.png')
    plt.close()

# parameter_perturbation(data = data, density = 58, instability = 'convective_numerical')

# Plotting the proportion of different instabilities at different densities 
def proportions(data, quantity):

    data = data[~data['dens'].isin([2, 6])]
    data = data[data['dens'] < 138]
    dens_array = np.unique(data['dens'])
    data['dv'] = data['da'] * data['dtau']
    data = data[data['dv'] < 10]
    dv_list = np.unique(data['dv'])

    stability_cases = {
        0: ('stable', 'm'),
        1: ('downstream', 'b'),
        2: ('absolute', 'r'),
        3: ('upstream', 'g'),
        'nan': ('accident', 'y')  # Adding NaN case
    }

    fig, ax = plt.subplots(figsize=(5, 4), dpi = 600)

    for count, dens in enumerate(dens_array):
        temp_data = data[data['dens'] == dens]

        for case, (label, color) in stability_cases.items():
            if case == 'nan':
                stability_data = temp_data[np.isnan(temp_data[quantity])]
            else:
                stability_data = temp_data[temp_data[quantity] == case]

            counts = [len(stability_data[stability_data['dv'] == dv]) / 100 for dv in dv_list]
            
            if counts:  # Avoid errors if counts is empty
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

def error_proportions(data):
    data['dv'] = data['da'] * data['dtau']
    data = data[~data['dens'].isin([2, 6])]
    data = data[data['dens'] < 138]

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
            
            temp_data.loc[
                ((temp_data['convective_analytical'] != 0) & (temp_data['convective_numerical'] == 0)) | 
                ((temp_data['convective_numerical'].isin([1, 3])) & (temp_data['convective_analytical'] == 2))
            ],  
            
            temp_data.loc[
                ((temp_data['convective_analytical'].isin([1, 3])) & (temp_data['convective_numerical'] == 2)) |
                ((temp_data['convective_analytical'].isin([1, 3])) & (np.isnan(temp_data['convective_numerical'])))
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
                
                merged = temper_data.merge(target_data[['sens_vel', 'sens_relvel', 'convective_numerical']],  # Only the column of interest from df2
                                    on=['sens_vel', 'sens_relvel'],
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

def stability_class(quantity):
    data = pd.read_csv(f'stability_diagram_{model}.csv')
    column_name = f'stability_{quantity}'
    print(data[column_name].unique())
    fig, ax = plt.subplots(figsize=(3, 3), dpi=600)
    for i in range(len(data)):
        if data.iloc[i][column_name] == 's':
            ax.scatter(data.iloc[i]['sens_relvel'], data.iloc[i]['sens_vel'], marker = 'p', alpha = 1, s = 20, color = 'g')
        elif data.iloc[i][column_name] in ['sna']:
            ax.scatter(data.iloc[i]['sens_relvel'], data.iloc[i]['sens_vel'], marker = '^', alpha = 1, s = 20, color = 'r')
        elif data.iloc[i][column_name] in ['snu', 'su']:
            ax.scatter(data.iloc[i]['sens_relvel'], data.iloc[i]['sens_vel'], marker = 'v', alpha = 1, s = 20, color = 'b')
        elif data.iloc[i][column_name] in ['snus', 'sds', 'sns', 'snds', 'sus', 'sdns']:
            ax.scatter(data.iloc[i]['sens_relvel'], data.iloc[i]['sens_vel'], marker = 'v', alpha = 1, s = 20, color = 'c')
        elif data.iloc[i][column_name] in ['snaus', 'saus', 'sus', 'suaus']:
            ax.scatter(data.iloc[i]['sens_relvel'], data.iloc[i]['sens_vel'], marker = '^', alpha = 1, s = 20, color = 'y')

    # analytical

    # if quantity == 'analytical':

    #     coordinates = [(8, 1), (9, 1), (10, 1), (10, 2), (3, 2), (4, 3), (5, 3), (3, 3), (3, 4), (3, 5), (2, 5), (2, 6), (2, 7), (2, 8), (2, 9), (2, 10), (1, 9), (1, 10)]
    #     ax = patching(ax, coordinates, 'b')

    # elif quantity == 'numerical':

    #     coordinates = [(1, 10), (1, 9), (2, 10), (2, 9), (2, 8), (2, 4), (2, 5), (3, 7), (3, 6), (3, 5), (3, 4), (4, 4), (4, 3), (5, 3), (9, 2), (10, 2)]
    #     ax = patching(ax, coordinates, 'b')

    # if quantity == 'analytical':

    #     coordinates = [(1, 6), (1, 7), (1, 8), (1,9), (1, 10), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8), (2, 9), (2, 10), (3, 3), (3, 4), (3, 5) ,(3, 6), (3, 7), (3, 8), (3, 9), (3, 10), (4, 2), (4, 3), (4, 4), (4, 5), (4, 6), (4, 7), (4, 8), (4, 9), (4, 10), (5, 1), (5, 2), (5, 4), (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6), (7, 1), (7, 2), (7, 3), (7, 4), (7, 5), (8, 1), (9, 1)]
    #     ax = patching(ax, coordinates, 'b')

    # elif quantity == 'numerical':

    #     coordinates = [(1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1,9), (1, 10), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8), (2, 9), (2, 10), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (3, 10), (4, 1), (4, 2), (4, 3), (4, 4), (4, 5), (4, 6), (4, 7), (4, 8), (4, 9), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (5, 6), (5, 7), (6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (7, 1), (7, 2), (7, 3), (8, 1)]
    #     ax = patching(ax, coordinates, 'b')


    plt.scatter(-1, -1, marker = 'v', s = 20, color = 'c', alpha = 1, label = '2b')
    # plt.scatter(-1, -1, marker = 'v', s = 20, color = 'b', alpha = 1, label = '2a')
    # plt.scatter(-1, -1, marker = '^', s = 20, color = 'r', alpha = 1, label = '1a')
    plt.scatter(-1, -1, marker = '^', s = 20, color = 'y', alpha = 1, label = '1b')
    plt.scatter(-1, -1, marker = 'p', s = 20, color = 'g', alpha = 1, label = '3')  

    plt.xlabel('$\\kappa_2$', fontsize = 15)
    plt.ylabel('$\\kappa_1$', fontsize = 15)
    plt.tick_params(axis = 'both', which = 'major', labelsize = 12)
    plt.xlim(0.15, 1.05)
    plt.ylim(0.15, 1.05)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=12)
    plt.savefig(f'{model}_stability_classes_{quantity}.png', bbox_inches = 'tight')
    plt.close()

# stability_class('numerical')
# stability_class('analytical')