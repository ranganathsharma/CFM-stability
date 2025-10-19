""" 
This file has the codes to simulate mixed flow based on the chosen set of models to plot the trajectories.
"""

# Imports

import os, sys, yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, root_scalar
from scipy import optimize
import logging
from matplotlib.collections import LineCollection
from matplotlib.ticker import ScalarFormatter

# Import the parameters

logging.basicConfig(filename= f'{__file__}.log', 
                    level=logging.WARNING, 
                    filemode = 'w')

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)

yaml_file_path = os.path.join(parent_dir, 'parameters.yaml')
with open(yaml_file_path, 'r') as file:
    model_param_dict = yaml.safe_load(file)

from model_equations import Model_equations
from solvers import solvers

print('The processing is complete')

class Mixed_flow():

    def __init__(self, 
                 model: str,
                 params: dict,
                 autonomous_percent: float, 
                 density: float,
                 da: float = 0,
                 dtau: float = 1):
        
        self.model = model
        self.percentage = autonomous_percent
        self.init_dens = density
        self.da = da
        self.dtau = dtau
        self.dt = 0.1
        self.num_veh = 600
        self.param_list = params
        self.param_list['m'] = 3
        self.param_list['spacing_weight'] = np.array([0.7, 0.2, 0.1])
        self.param_list['relvel_weight'] = np.array([0.7, 0.2, 0.1])
        self.threshold = 1e-5
        self.unpack_params()

    def unpack_params(self):
        """
        Method to input the parameter values
        """
        for key in self.param_list.keys():
            setattr(self, key, self.param_list[key]) 

    def equilibrium_vel(self, 
                        dx: float) -> None:
        """
        Method to calculate the equilibrium velocity for the given initial gap

        Parameters
        ----------
        dx : float
            Initial spacing between vehicles including the length of the vehicle
        """
        self.equations = Model_equations(self.param_list)
        init_vel_function = getattr(self.equations, 'init_vel_' + str(self.model), None) 
        
        if callable(init_vel_function):
            init_vel = init_vel_function(dx)
        else:
            logging.error(f'The initial velocity cannot be computed at {dx} for the model {self.model}')

        self.init_vel = init_vel

    def setup_arrays(self):

        """
        Based on the list of vehicles, set the corresponding equations as an array
        Solve all the equations as a vector function in a time loop instead of solving each of them in a loop
        """
        sim_dur = 600
        self.tim_array = np.arange(0, sim_dur, self.dt)
        self.pos_array = np.zeros((self.num_veh, int(sim_dur/self.dt)))
        self.vel_array = np.zeros_like(self.pos_array)

        # Generate a list of vehicle types where the first three are human driven (0) and the rest are picked according to the autonomous percentage (1 for autonomous)
        vehicle_types = [0] * 3  # First three are human driven
        num_autonomous = int((self.num_veh - 3) * self.percentage)
        num_human = (self.num_veh - 3) - num_autonomous
        rest_types = [1] * num_autonomous + [0] * num_human
        np.random.shuffle(rest_types)
        vehicle_types += rest_types
        self.vehicle_types = vehicle_types

        init_gap = 1000/self.init_dens

        if init_gap < model_param_dict['Min_head'][self.model]:
            logging.warning(f'The initial spacing is lower than the threshold for the model {self.model}')
            raise ValueError() # Value error when the system cannot be processed
        else:
            pass

        try:
            self.equilibrium_vel(init_gap)
        except:
            logging.error(f'The equilibrium velocity could not be calculated for init_gap = {init_gap}')
            raise ValueError() # Value error when the system cannot be processed

        self.init_gap = init_gap

        if self.init_vel < 0:
            logging.warning(f'The equilibrium velocity is predicted to be negative for init_gap = {init_gap}')
            raise ValueError() # Value error when the system cannot be processed

        self.vel_array[:,0] = self.init_vel # Initial condition for all vehicles
        self.vel_array[0,:] = self.init_vel # Boundary condition for the leader vehicle. Has to be changed based on the perturbation kind
        
        for i in range(self.num_veh):
            self.pos_array[i,0] = -i*self.init_gap

        if self.da + self.dtau != 0:
            da = self.da
            dur = self.dtau
        else:
            logging.warning('The input for the kind of perturbation is invalid')

        temp_time = np.arange(0, 2*dur, self.dt)
        temp_vel = self.init_vel*np.ones((temp_time.shape[0]))

        for i in range(1, temp_time.shape[0]):
            if temp_time[i] < dur + self.dt:
                temp_vel[i] = temp_vel[i-1] - da*self.dt
            else:
                temp_vel[i] = temp_vel[i-1] + da*self.dt

        temp_vel = np.where(temp_vel < 0, 0, temp_vel)

        self.vel_array[0, int(60/self.dt): int(60/self.dt) + temp_time.shape[0]] = temp_vel

    def run_simulation(self):

        self.setup_arrays() 
        solver_class = solvers()
        numerical_solver = getattr(solver_class, 'euler_second_TK')

        for i in range(1, self.num_veh):
            print(i)

            if self.vehicle_types[i] == 0:
                model_name = self.model
            else:
                model_name = 'm' + self.model
                self.m = self.param_list['m']

            self.equation = Model_equations(self.param_list)
            v_equation = getattr(self.equation, 'u_dot_' + model_name + '_open')
            x_equation = getattr(self.equation, 'x_dot_' + model_name + '_open')

            if model_name[0] == 'm':
                a_old = np.zeros((self.m+1))
                for count, time in enumerate(self.tim_array[:-1]):
                    x_new, v_new, a_new = numerical_solver(self.pos_array[i-self.m: i+1, count],
                                                           self.vel_array[i-self.m: i+1, count],
                                                           a_old,
                                                           self.dt,
                                                           v_equation,
                                                           x_equation)
                    a_old = a_new
                    a_old[:self.m] = 0

                    if np.min(x_new[:-1] - x_new[1:]) < 0:
                        logging.warning(f'The vehicles have crashed for the vehicle {np.where(x_new[:-1] - x_new[1:] == np.min(x_new[:-1] - x_new[1:]))}')
                        raise OverflowError()
                    
                    if np.min(v_new[1:]) < 0:
                        logging.warning(f'The vehicles are moving back for the vehicle {np.where(v_new[1:] == np.min(v_new[1:]))}')
                        raise OverflowError()
                        
                    # Update the new positions and velocities
                    self.pos_array[i, count + 1], self.vel_array[i, count + 1] = x_new[-1], np.maximum(v_new[-1], 0)

            else:
                a_old = 0
                for count, time in enumerate(self.tim_array[:-1]):
                    x_new, v_new, a_new = numerical_solver(self.pos_array[i-1:i+1, count],
                                                            self.vel_array[i-1:i+1, count],
                                                            a_old, 
                                                            self.dt,
                                                            v_equation,
                                                            x_equation)
                        
                    a_old = a_new
                    a_old[0] = 0

                    # print(self.pos_array[i-1:i+1, count])

                    if np.min(x_new[:-1] - x_new[1:]) < 0:
                        logging.warning(f'The vehicles have crashed for the vehicle {np.where(x_new[:-1] - x_new[1:] == np.min(x_new[:-1] - x_new[1:]))}')
                        logging.warning(f'The vehicle that has crashed is {i}')
                        logging.warning(f'The type of the vehicle is {model_name}')
                        raise OverflowError()
                    
                    if np.min(v_new[1:]) < 0:
                        logging.warning(f'The vehicles are moving back for the vehicle {np.where(v_new[1:] == np.min(v_new[1:]))}')
                        raise OverflowError()
                        
                    # Update the new positions and velocities
                    if i!=1:
                        self.pos_array[i, count + 1], self.vel_array[i, count + 1] = x_new[1], np.maximum(v_new[1], 0)
                    else:
                        self.pos_array[i-1:i+1, count + 1], self.vel_array[i, count + 1] = x_new, np.maximum(v_new[-1], 0)
        
        plt.figure(dpi=600)

        lc = None  

        vmin = np.min(self.vel_array)
        vmax = np.max(self.vel_array)

        for i in range(1, self.num_veh, 5):

            x = self.tim_array/60
            y = self.pos_array[i]/1000     

            values = self.vel_array[i]

            points = np.array([x, y]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

            lc = LineCollection(segments, cmap="rainbow", linewidth=0.8,
                                norm=plt.Normalize(vmin=vmin, vmax=vmax))
            lc.set_array(values)
            plt.gca().add_collection(lc)

        cbar = plt.colorbar(lc)
        cbar.set_label(f"velocity (m/s)", rotation=270, labelpad=20, fontsize=15)
        cbar.ax.tick_params(labelsize=15)

        fmt = ScalarFormatter(useMathText=True)
        fmt.set_scientific(True)
        fmt.set_powerlimits((-2, 2))   # controls when sci-notation kicks in
        cbar.ax.yaxis.set_major_formatter(fmt)

        plt.tick_params(axis = 'both', which = 'major', labelsize = 15)
        plt.ylim(-4, 0.2)
        plt.xlim(0.5, 7)
        plt.xlabel("Time (min)" , fontsize = 18)
        plt.ylabel("Position (km)", fontsize = 18)
        plt.tight_layout()
        # plt.show()
        plt.savefig(f'{self.model}_{self.percentage}_trajectory.png', bbox_inches = 'tight')
        plt.close()

if __name__ == '__main__':

    dens = 130
    
    pct_list = [0, 0.25, 0.5, 1]

    for pct in pct_list:

        params = model_param_dict['Models']['idm'].copy()
        params['safe_head'] = 0.96 
        params['max_acc'] = 1.5
        sim = Mixed_flow(model='idm', 
                        params=params,
                        autonomous_percent = pct, 
                        density=dens, 
                        da=0.1, 
                        dtau=1)
        sim.run_simulation()
