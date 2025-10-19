"""
This file has the codes to calculate stability conditions for a given model at different traffic states and perturbation strengths.
The size of the platoon and duration of simulation are the hyperparameters of the system
"""

# Imports

import os, sys, yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, root_scalar
from scipy import optimize
import logging
import multiprocessing as mp

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

class Stability_analysis():

    def __init__(self, 
                 model: str,
                 dens: float,
                 pert_kind: str,
                 dt: float = 0.1, 
                 da: float = 0,
                 dtau: float = 0) -> None:
        """
        Class to catalogue the types of instabilities observed in the CFM

        Parameters
        ----------
        model : str
            name of the model
        dens : float
            density of the initial positions in vehs/m
        pert_kind : str
            strength of the perturbation: small, medium, large, custom
        dt : float, optional
            simulation time step, by default 0.1
        da: float, optional
            perturbation strength by acceleration, by default 0
        dtau: float, optional
            perturbation strength by duration, by default 0
        """
        self.model = model
        self.init_dens = dens
        self.pert_kind = pert_kind
        self.dt = dt
        self.num_veh = 1000
        self.da = da
        self.dtau = dtau
        self.param_list = model_param_dict['Models'][model]
        self.unpack_params()
        self.model_order = model_param_dict['Model_order'][model]
        self.threshold = 1e-5

        # Create the file with the header
        file_path = os.path.dirname(os.path.abspath(__file__))
        file_path = os.path.join(file_path, f'{model}_{dens}_stability_status.csv')
        if not os.path.exists(file_path):

            temp_dict = {'n': [],
                        'dens': [],
                        'max_acc': [],
                        'safe_head': [],
                        'pert': [], 
                        'string_analytical': [], 
                        'string_numerical1': [],
                        'string_numerical2': [],
                        'string_numerical_value': [],
                        'convective_analytical': [],
                        'convective_numerical': [],
                        'downstream': [],
                        'upstream': [],
                        'sigma_0': [],
                        'v_g':[],
                        'upper_lim':[],
                        'da': [],
                        'dtau': []}
            
            pd.DataFrame(temp_dict).to_csv(file_path, mode = 'w', header = True, index = False)
        else:
            pass
    
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

    def string_stability_numerical(self) -> int:
        """
        Method to verify the string stability numerically

        Returns
        -------
        int
            values indicate different types of string instabilities observed. 
            0: stable
            1: velocity range is bounded by leader
            2: velocity range is not bounded by the leader, velocity drop is bounded by the leader
            3: velocity drop is not bounded by the leader
        """

        vel_max = np.max(self.vel_array[:,int(50/self.dt):int(1000/self.dt)], axis = 1)
        vel_min = np.min(self.vel_array[:,int(50/self.dt):int(1000/self.dt)], axis = 1)
        vel_drop = vel_max - vel_min

        vel_dev = self.init_vel - vel_min

        delta = vel_drop[1:] - vel_drop[:-1] 

        if np.all(delta < -self.threshold):
            stab_status = 0 # Stable
        else:
            if vel_drop[0] >= np.max(vel_drop[1:]):
                stab_status = 1 # Velocity drop increases
            else:
                if vel_dev[0] >= np.max(vel_dev[1:]):
                    stab_status = 2 # Perturbation capped
                else:
                    stab_status = 3 # Uncapped range

        vel_dev = np.max(np.abs(self.vel_array[:,int(50/self.dt):int(1000/self.dt)] - self.init_vel), axis = 1)
        delta = vel_dev[1:] - vel_dev[:-1]

        self.string_numerical_value = np.max(delta)
        
        if np.all(delta < -self.threshold):
            base_status = 0 # Stable
        else:
            base_status = 1 # Unstable

        return stab_status, base_status
                
    def string_stability_analytical(self, 
                                    model: str,
                                    eq_gap: float, 
                                    eq_vel: float) -> int:
        """
        Analytical expression to calculate the string stability of a model. 
        This function is particular to the model

        Parameters
        ----------
        model: str
            name of the model for verification
        eq_gap : float
            Gap between vehicles in equilibrium
        eq_vel : float
            Velocity of all the vehicles in equilibrium

        Returns
        -------
        int
            different values indicating analytical string stability status
            0: string stable
            1: string unstable
        """
        # The expressions for string stability is particular to the model. This method contains the expressions for the stability equation and gives the string stability status
        
        if model == 'idm':
            eq_gap -= self.len_veh

            f_s = 2* self.max_acc/eq_gap**3 *(self.jam_dis_0 + self.safe_head*eq_vel)**2
            f_v = - self.max_acc *(4/self.des_vel**4*eq_vel**3 + 2*self.safe_head*(self.jam_dis_0 + self.safe_head*eq_vel)/eq_gap**2)
            f_dv = np.sqrt(self.max_acc/self.des_dec)*eq_vel/eq_gap**2 *(self.jam_dis_0 + self.safe_head*eq_vel)

            if 0.5 - f_dv/f_v - f_s/f_v**2 > 0:
                return 0 # Stable
            else:
                return 1 # Unstable
        else:
            logging.warning(f'The string stability conditions are not applied for the right model')

    def convective_string_stability_numerical(self) -> int:
        """
        Method to calculate the convective string stability in the model

        Returns
        -------
        int
            values indicate different types of convective instabilities
            0: convective stable
            1: downstream convective instability
            2: absolute convective instability
            3: upstream convective instability
        """
        self.pert_veh = 0

        dx, dt = 100, 10
        dv = np.abs(self.vel_array - self.init_vel)
        aggregate_df = pd.DataFrame()
        aggregate_df['x'] = pd.DataFrame(self.pos_array[1:].flatten())
        aggregate_df['dv'] = pd.DataFrame(dv[1:].flatten())
        big_t = np.ones((len(self.pos_array), len(self.pos_array[0])))
        big_t*= self.tim_array
        aggregate_df['t'] = pd.DataFrame(big_t.flatten())
        aggregate_df['x_index'] = aggregate_df['x']//dx
        aggregate_df['t_index'] = aggregate_df['t']//dt

        start_time, end_time = self.tim_array[0], self.tim_array[-1]
        time_range = end_time - start_time
        start_pos, end_pos = np.min(self.pos_array), np.max(self.pos_array)
        pos_range = end_pos - start_pos
        del big_t
    
        grouped = aggregate_df.groupby(['x_index', 't_index'])
        del aggregate_df

        space_0, time_0 = start_pos//dx, start_time
        Velocity_matrix = np.zeros((int(np.ceil(pos_range/dx) + 1), int(np.ceil(time_range/dt) + 1)))

        for (space, time), group_df in grouped:
            Velocity_matrix[int(np.floor(space-space_0)), int(np.floor(time-time_0))] += group_df.dv.max()

        pert_start_index = int((self.pos_array[self.pert_veh, int(60/self.dt)]-start_pos)/dx)
        pert_end_index = int((self.pos_array[self.pert_veh, int((60 + 2*self.dtau)/self.dt)]-start_pos)/dx)

        Velocity_matrix = Velocity_matrix[:-2, :-2]
        self.grouped = grouped
        self.space_0 = space_0
        self.time_0 = time_0
        self.start_pos, self.end_pos = start_pos, end_pos
        self.start_time, self.end_time = start_time, end_time

        self.upstream = 0
        self.downstream = 0
        self.upstream_value = 0
        self.downstream_value = 0

        upstream = np.max(Velocity_matrix[:pert_end_index], axis=1)
        upstream[-1] = self.da * self.dtau
        upstream_diff = upstream[:-1] - upstream[1:] # pert upstream - pert downstream; if > 0: unstable

        if np.any(upstream_diff > self.threshold):
            self.upstream_value = np.max(upstream_diff)
            self.upstream = 1

        elif np.all(upstream_diff[-10:] == 0) and (np.all(upstream[-10:]) > self.threshold):
            self.upstream = 1

        downstream = np.max(Velocity_matrix[pert_start_index:], axis = 1)
        downstream[0] = self.da * self.dtau
        downstream_diff = downstream[:-1] - downstream[1:] # pert upstream - pert downstream; if  < 0: unstable

        del Velocity_matrix

        if np.any(downstream_diff < -self.threshold):
            self.downstream_value = np.min(downstream_diff)
            self.downstream = 1

        elif np.all(downstream_diff[:10] == 0) and (np.all(downstream[:10]) > self.threshold):
            self.downstream = 1

        del Velocity_matrix

        if self.upstream + self.downstream == 2:
            return 2
        if self.upstream == 1:
            return 3
        if self.downstream == 1:
            return 1
        else:
            return 0
   
    def convective_string_stability_analytical(self, 
                                                eq_gap: float,
                                                eq_vel: float) -> int:
        """
        Method to calculate the convective string stability of a model analytically

        Parameters
        ----------
        eq_gap : float
            Gap between vehicles in equilibrium
        eq_vel : float
            Velocity of all vehicles in equilibrium

        Returns
        -------
        int
            values indicating different types of convective string instabilities
            0: Stable
            1: downstream convective string instability
            2: absolute convective string instability
            3: upstream convective string instability
        """
        eq_gap -= self.len_veh

        f_s = 2* self.max_acc/eq_gap**3 *(self.jam_dis_0 + self.safe_head*eq_vel)**2
        f_v = - self.max_acc *(4/self.des_vel**4*eq_vel**3 + 2*self.safe_head*(self.jam_dis_0 + self.safe_head*eq_vel)/eq_gap**2)
        f_dv = np.sqrt(self.max_acc/self.des_dec)*eq_vel/eq_gap**2 *(self.jam_dis_0 + self.safe_head*eq_vel)

        def lambda_func(kappa):
            p = -f_v + f_dv*(1 - np.exp(-1j*kappa))
            q = f_s*(1 - np.exp(-1j*kappa))

            return (-p + np.sqrt(p**2 - 4*q))*0.5 # Assuming only the positive branch is the unstable root
        
        def lambda_prime(kappa):
            
            exp_term = np.exp(-1j*kappa)
            return -0.5*1j*f_dv*exp_term + 0.5*(1j*f_dv*(f_dv*(1 - exp_term) - f_v)*exp_term - 2*1j*f_s*exp_term)/np.sqrt(-4*f_s*(1 - exp_term) + (f_dv*(1 - exp_term) - f_v)**2)
            
        def lambda_pprime(kappa):
            exp_term = np.exp(-1j*kappa)
            return -0.5*f_dv*exp_term + 0.5*(-f_dv**2*exp_term**2 + f_dv*(f_dv*(1 - exp_term) - f_v)*exp_term - 2*f_s*exp_term)/np.sqrt(-4*f_s*(1 - exp_term) + (f_dv*(1 - exp_term) - f_v)**2) + 0.5*(-1j*f_dv*(f_dv*(1 - exp_term) - f_v)*exp_term + 2*1j*f_s*exp_term)*(1j*f_dv*(f_dv*(1 - exp_term) - f_v)*exp_term - 2*1j*f_s*exp_term)/(-4*f_s*(1 - exp_term) + (f_dv*(1 - exp_term) - f_v)**2)**(3/2)
        
        kappa_list = np.linspace(0, np.pi, 100)
        lambda_list = np.zeros_like(kappa_list)
        for count, kappa in enumerate(kappa_list):
            lambda_list[count] = lambda_func(kappa).real

        count = np.where(lambda_list == np.max(lambda_list))[0][-1]

        try:
            result = optimize.minimize_scalar(lambda x: -lambda_func(x).real, bounds=(kappa_list[count - 1], kappa_list[count + 1]), method='bounded')
        except: 
            # optimize.minimize_scalar raises an error when the first value is the largest making count - 1 = -1 i.e. last element of kappa_list
            if count == 0:
                result = optimize.minimize_scalar(lambda x: -lambda_func(x).real, bounds=(0, kappa_list[count + 1]), method='bounded')
            else:
                result = optimize.minimize_scalar(lambda x: -lambda_func(x).real, bounds=(kappa_list[count-1], np.pi), method='bounded')
        
        kappa_0 = result.x
        lambda_0 = lambda_func(kappa_0)
        sigma_0 = lambda_0.real
        rho_e = 1/(eq_gap + self.len_veh)

        v_g = eq_vel + lambda_prime(kappa_0).imag/rho_e
        sigma_kk = lambda_pprime(kappa_0).real/rho_e**2
        omega_kk = lambda_pprime(kappa_0).imag/rho_e**2

        D_2 = -sigma_kk*(1 + omega_kk**2/sigma_kk**2)

        self.sigma_0 = sigma_0
        self.v_g = v_g
        self.upper_lim = v_g**2/2/D_2

        try:
            cs_1 = v_g + np.sqrt(2*D_2*sigma_0)
            cs_2 = v_g - np.sqrt(2*D_2*sigma_0)
        except:
            return 0

        if (0 < sigma_0) & (sigma_0 <= self.upper_lim):

            if (cs_1 < 0) & (cs_2 < 0):
                return 3 #upstream
            elif (cs_1 > 0) & (cs_2 > 0):
                return 1 #downstream
            else:
                logging.warning(f'The convective instability could not be calculated analytically')
        else:
            if cs_1*cs_2 < 0:
                return 2 #absolute
            else:
                return 0

    def setup_arrays(self,
                     init_gap: float) -> None:
        """
        Method to set up the arrays 

        Parameters
        ----------
        init_gap : float
            initial distance between vehicles including the length of the vehicle
        """
        
        sim_dur = 3600

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

        self.tim_array = np.arange(0, sim_dur, self.dt)
        self.pos_array = np.zeros((self.num_veh, self.tim_array.shape[0]))
        self.vel_array = np.zeros((self.num_veh, self.tim_array.shape[0]))

        self.vel_array[:,0] = self.init_vel # Initial condition for all vehicles
        self.vel_array[0,:] = self.init_vel # Boundary condition for the leader vehicle. Has to be changed based on the perturbation kind
        
        for i in range(self.num_veh):
            self.pos_array[i,0] = -i*self.init_gap

        if self.pert_kind == 'small':
            da = 0.1
            dur = 1
        elif self.pert_kind == 'medium':
            da = 1
            dur = 3
        elif self.pert_kind == 'large':
            da = 2
            dur = 5
        elif self.pert_kind == 'custom':
            if self.da + self.dtau != 0:
                da = self.da
                dur = self.dtau
            else:
                logging.warning('The input for the kind of perturbation is invalid')
        else:
            logging.warning('Invalid input for perturbation')

        temp_time = np.arange(0, 2*dur, self.dt)
        temp_vel = self.init_vel*np.ones((temp_time.shape[0]))

        for i in range(1, temp_time.shape[0]):
            if temp_time[i] < dur + self.dt:
                temp_vel[i] = temp_vel[i-1] - da*self.dt
            else:
                temp_vel[i] = temp_vel[i-1] + da*self.dt

        temp_vel = np.where(temp_vel < 0, 0, temp_vel)

        self.vel_array[0, int(60/self.dt): int(60/self.dt) + temp_time.shape[0]] = temp_vel

    def evolve(self, 
               params: dict) -> None:

        init_gap = 1000/self.init_dens
        self.param_list = params
        self.unpack_params()
        self.setup_arrays(init_gap = init_gap)
        self.equation = Model_equations(self.param_list)

        if self.model_order == 'second':
            v_equation = getattr(self.equations, 'u_dot_' + self.model + '_open')
            x_equation = getattr(self.equations, 'x_dot_' + self.model + '_open')
        else:
            v_equation = getattr(self.equations, 'x_dot_' + self.model + '_open')
            x_equation = getattr(self.equations, 'x_dot_' + self.model + '_open')

        solver_class = solvers()

        numerical_solver = getattr(solver_class, 'euler_second_TK')

        a_old = np.zeros((self.pos_array.shape[0]))
        for count, time in enumerate(self.tim_array[:-1]):
            x_new, v_new, a_new = numerical_solver(self.pos_array[:,count],
                                                   self.vel_array[:,count],
                                                   a_old,
                                                   self.dt,
                                                   v_equation,
                                                   x_equation)
            a_old = a_new
            a_old[0] = 0 # Leader dynamics is always controlled

            if np.max(v_new) > 50:
                logging.warning(f'The velocity has blown up for {np.where(v_new) == np.max(v_new)} with a value {np.max(v_new)} at the count {count}')
                raise OverflowError()

            if np.min(x_new[:-1] - x_new[1:]) < 0:
                logging.warning(f'The vehicles have crashed for the vehicle {np.where(x_new[:-1] - x_new[1:] == np.min(x_new[:-1] - x_new[1:]))}')
                raise OverflowError()
            
            if np.min(v_new[1:]) < 0:
                logging.warning(f'The vehicles are moving back for the vehicle {np.where(v_new[1:] == np.min(v_new[1:]))}')
                raise OverflowError()
                
            # Update the new positions and velocities
            self.pos_array[:, count + 1], self.vel_array[1:, count + 1] = x_new, np.maximum(v_new[1:], v_new[1:]*0)
        
    def looper(self):
        
        num_params = 10
        safe_head_array = np.linspace(0.1, 4, num_params)
        max_acc_array = np.linspace(0.1, 4, num_params)

        temp_params = self.param_list

        for i, safe_head in enumerate(safe_head_array):
            temp_params['safe_head'] = safe_head
            for j, max_acc in enumerate(max_acc_array):
                temp_params['max_acc'] = max_acc
                
                self.processor(params = temp_params, save = True)
        
    def processor(self, params: dict,
                  save: bool = False):
        try:
            self.evolve(params = params)
            str_numerical_status1, str_numerical_status2 = self.string_stability_numerical()
        except (OverflowError, ValueError):
            logging.warning(f'There was a crash at parameter values {params}, da, dtau, dens =  {self.da, self.dtau, self.init_dens}')
            str_numerical_status1 = np.nan
            str_numerical_status2 = np.nan 
            self.string_numerical_value = np.nan

        str_analytical_status = self.string_stability_analytical(model= self.model,
                                                                eq_gap = self.init_gap, 
                                                                eq_vel = self.init_vel)  

        if str_analytical_status != 0:
            con_analytical_status = self.convective_string_stability_analytical(eq_gap = self.init_gap,
                                                                                eq_vel = self.init_vel)                    
        else:
            con_analytical_status = 0
            self.v_g = np.nan
            self.sigma_0 = np.nan
            self.upper_lim = np.nan

        if not np.any(np.isnan([str_numerical_status1, str_numerical_status2])):
            con_numerical_status = self.convective_string_stability_numerical()
        else:
            self.upstream_value = np.nan
            self.downstream_value = np.nan
            con_numerical_status = np.nan


        temp_dict = {'n': [self.num_veh],
                    'dens': [1000/self.init_gap],
                    'max_acc': [params['max_acc']],
                    'safe_head': [params['safe_head']],
                    'pert': [self.pert_kind], 
                    'string_analytical': [str_analytical_status], 
                    'string_numerical1': [str_numerical_status1],
                    'string_numerical2': [str_numerical_status2],
                    'string_numerical_value': [self.string_numerical_value],
                    'convective_analytical': [con_analytical_status],
                    'convective_numerical': [con_numerical_status],
                    'downstream':[self.downstream_value],
                    'upstream': [self.upstream_value],
                    'sigma_0' : [self.sigma_0],
                    'v_g': [self.v_g], 
                    'upper_lim': [self.upper_lim]}
        
        if self.pert_kind == 'small':
            temp_dict['da'] = [0.1]
            temp_dict['dtau'] = [1]
        elif self.pert_kind == 'medium':
            temp_dict['da'] =  [1]
            temp_dict['dtau'] = [3]
        elif self.pert_kind == 'large':
            temp_dict['da'] = [2]
            temp_dict['dtau'] = [5]
        elif self.pert_kind == 'custom':
            temp_dict['da'] = [self.da]
            temp_dict['dtau'] = [self.dtau]
        else:
            logging.warning('The input for the perturbation is invalid at ')
        
        if save:
            pd.DataFrame(temp_dict).to_csv(f'{self.model}_{self.init_dens}_stability_status.csv', mode = 'a', header = False, index=False)
        else:
            pass