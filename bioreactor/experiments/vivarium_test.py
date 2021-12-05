import os
import copy
import pylab as plt
import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
import matplotlib

# Process, Deriver, and Composer base classes
from vivarium.core.process import Process, Deriver
from vivarium.core.composer import Composer
from vivarium.core.registry import process_registry

# helper functions for composition
from vivarium.core.composition import (
    simulate_process,
    process_in_experiment,
    composer_in_experiment,
    composite_in_experiment)

# other vivarium imports
from vivarium.core.engine import Engine, pp
from vivarium.library.units import units

# plotting functions
from vivarium.plots.simulation_output import (
    plot_simulation_output, plot_variables)
from vivarium.plots.simulation_output import _save_fig_to_dir as save_fig_to_dir
from vivarium.plots.agents_multigen import plot_agents_multigen
from vivarium.plots.topology import plot_topology


### setup

OUTDIR = 'out/vivarium_test'
os.makedirs(OUTDIR, exist_ok=True)

AVOGADRO = constants.N_A * 1 / units.mol

store_cmap = matplotlib.cm.get_cmap('Dark2')
dna_color = matplotlib.colors.to_rgba(store_cmap(0))
rna_color = matplotlib.colors.to_rgba(store_cmap(1))
protein_color = matplotlib.colors.to_rgba(store_cmap(2)) 
global_color = matplotlib.colors.to_rgba(store_cmap(7))                             
store_colors = {
    'DNA': dna_color, 
    'DNA\n(counts)': dna_color, 
    'DNA\n(mg/mL)': dna_color, 
    'mRNA': rna_color,
    'mRNA\n(counts)': rna_color,
    'mRNA\n(mg/mL)': rna_color,
    'Protein': protein_color,
    'Protein\n(mg/mL)': protein_color,
    'boundary': global_color}
                       
# plotting configurations
topology_plot_config = {
    'settings': {
        'coordinates': {
            'Tl': (-1,0),
            'Tx': (-1,-1),
            'Protein': (1,0),
            'mRNA': (1,-1),
            'DNA': (1,-2),
        },
        'node_distance': 3.0,
        'process_color': 'k',
        'store_colors': store_colors,
        'dashed_edges': True,
        'graph_format': 'vertical',
        'color_edges': False}, 
    'out_dir': OUTDIR}

plot_var_config = {
    'row_height': 2,
    'row_padding': 0.2,
    'column_width': 10,
    'out_dir': OUTDIR}

stoch_exp_settings = {
    'settings': {
        'experiment_id': 'stochastic_txtl'},
    'initial_state': {
            'DNA\n(counts)': {
                'G': 1.0
            },
            'mRNA\n(counts)': {
                'C': 0.0
            },
            'Protein\n(mg/mL)': {
                'X': 0.0 * units.mg / units.mL
            }}}

stoch_plot_config = {
    'variables':[
        {
            'variable': ('Protein\n(mg/mL)', ('X', 'milligram / milliliter')), 
            'color': store_colors['Protein'],
            'display': 'Protein: X (mg/mL)'}, 
        {
            'variable': ('mRNA\n(mg/mL)', ('C', 'milligram / milliliter')), 
            'color': store_colors['mRNA'],
            'display': 'mRNA: C (mg/mL)'}, 
        {
            'variable': ('DNA\n(counts)', 'G'), 
            'color': store_colors['DNA'],
            'display': 'DNA: G (counts)'}, 
    ],
    'filename': 'stochastic_txtl_output.pdf',
    **plot_var_config}

# plot config
hierarchy_plot_settings = {
    'include_paths': [
        ('RNA_counts', 'C'),
        ('RNA', 'C'),
        ('Protein', 'X'),
        ('DNA', 'G'),
    ],
    'store_order': ('Protein', 'RNA_counts', 'RNA', 'DNA'),
    'titles_map': {
        ('Protein', 'X',): 'Protein', 
        ('RNA_counts', 'C'): 'RNA\n(counts)', 
        ('DNA', 'G',): 'DNA', 
        ('RNA', 'C'): 'RNA\n(concentration)',
    },
    'column_width': 10,
    'row_height': 1.5,
    'stack_column': True,
    'tick_label_size': 10,
    'linewidth': 1.5,
    'title_size': 10}

colony_plot_config = {
    'variables': [('boundary', ('volume', 'femtoliter'))],
    'filename': 'colony_growth.pdf',
    **plot_var_config}


# hierarchy topology plot
agent_0_string = 'agents\n0'
agent_1_string = 'agents\n00'
agent_2_string = 'agents\n01'
row_1 = 0
row_2 = -1
row_3 = -2
row_4 = -3
node_space = 0.75
vertical_space=0.9
bump = 0.1
process_column = -0.2
agent_row = -3.2
agent_column = bump/2 #0.5

hierarchy_topology_plot_config = {
    'settings': {
        'graph_format': 'hierarchy',
        'node_size': 6000,
        'process_color': 'k',
        'store_color': global_color,
        'store_colors': {
            f'{agent_0_string}\nDNA': dna_color,
            f'{agent_0_string}\nRNA': rna_color,
            f'{agent_0_string}\nRNA_counts': rna_color,
            f'{agent_0_string}\nProtein': protein_color,
        },
        'dashed_edges': True,
        'show_ports': False,
        'coordinates': {
            # Processes
            'ColonyVolume': (2.5, 0),
            'agents\n0\nstochastic_Tx': (agent_column, agent_row*vertical_space),
            'agents\n0\nTl': (agent_column+node_space, agent_row*vertical_space),
            'agents\n0\nconcs': (agent_column+2*node_space, agent_row*vertical_space),
            'agents\n0\ndivision': (agent_column+3*node_space, agent_row*vertical_space),
            # Stores
            'agents': (1.5*node_space, row_1*vertical_space),
            'agents\n0': (1.5*node_space, row_2*vertical_space),
            'agents\n0\nagents': (1.5*node_space, row_1*vertical_space),
            'agents\n0\nDNA': (0, row_3*vertical_space),
            'agents\n0\nRNA_counts': (node_space+bump, row_3*vertical_space),
            'agents\n0\nRNA': (node_space, (row_3-bump)*vertical_space),
            'agents\n0\nProtein': (2*node_space+bump, row_3*vertical_space),
            'agents\n0\nboundary': (3*node_space+bump, row_3*vertical_space),
        },
        'node_labels': {
            # Processes
            'ColonyVolume': 'Colony\nVolume',
            'agents\n0\nstochastic_Tx': 'stochastic\nTx',
            'agents\n0\nTl': 'Tl',
            'agents\n0\nconcs': 'counts\nto\nmg/mL',
            'agents\n0\ngrowth': 'growth',
            'agents\n0\ndivision': 'division',
            # Stores
            # third
            'agents\n0': '0',
            'agents\n0\nDNA': 'DNA',
            'agents\n0\nRNA': 'RNA',
            'agents\n0\nrates': 'rates',
            'agents\n0\nRNA_counts': '',
            'agents\n0\nboundary': 'boundary',
            'agents\n0\nProtein': 'Protein',
            # fourth
            'agents\n0\nrates\ngrowth_rate': 'growth_rate',
            'agents\n0\nrates\ngrowth_noise': 'growth_noise',
        },
        'remove_nodes': [
            'agents\n0\nshape',
            'agents\n0\nrates\ngrowth_rate',
            'agents\n0\nrates\ngrowth_noise',
            'agents\n0\nrates',
            'agents\n0\ngrowth',
            'agents\n0\ndivide_condition',
            'agents\n0\nboundary\ndivide',
            'agents\n0\nboundary\nvolume',
        ]
    },
    'out_dir': OUTDIR
}

# topology plot config for after division
agent_2_dist = 3.5
hierarchy_topology_plot_config2 = copy.deepcopy(hierarchy_topology_plot_config)

# redo coordinates, labels, store_colors, and removal
hierarchy_topology_plot_config2['settings']['node_distance'] = 2.5
hierarchy_topology_plot_config2['settings']['coordinates'] = {}
hierarchy_topology_plot_config2['settings']['node_labels'] = {}
hierarchy_topology_plot_config2['settings']['store_colors'] = {}
# hierarchy_topology_plot_config2['settings']['remove_nodes'] = []
for node_id, coord in hierarchy_topology_plot_config['settings']['coordinates'].items():
    if agent_0_string in node_id:
        new_id1 = node_id.replace(agent_0_string, agent_1_string)
        new_id2 = node_id.replace(agent_0_string, agent_2_string)
        hierarchy_topology_plot_config2['settings']['coordinates'][new_id1] = coord
        hierarchy_topology_plot_config2['settings']['coordinates'][new_id2] = (coord[0]+agent_2_dist, coord[1])
    else:
        hierarchy_topology_plot_config2['settings']['coordinates'][node_id] = (coord[0]+agent_2_dist/2, coord[1])
hierarchy_topology_plot_config2['settings']['coordinates']['ColonyVolume'] = (5.5, 0)

for node_id, label in hierarchy_topology_plot_config['settings']['node_labels'].items():
    if agent_0_string in node_id:
        new_id1 = node_id.replace(agent_0_string, agent_1_string)
        new_id2 = node_id.replace(agent_0_string, agent_2_string)
        hierarchy_topology_plot_config2['settings']['node_labels'][new_id1] = label
        hierarchy_topology_plot_config2['settings']['node_labels'][new_id2] = label
    else:
        hierarchy_topology_plot_config2['settings']['node_labels'][node_id] = label
hierarchy_topology_plot_config2['settings']['node_labels']['agents\n00'] = '1'
hierarchy_topology_plot_config2['settings']['node_labels']['agents\n01'] = '2'

for node_id, color in hierarchy_topology_plot_config['settings']['store_colors'].items():
    if agent_0_string in node_id:
        new_id1 = node_id.replace(agent_0_string, agent_1_string)
        new_id2 = node_id.replace(agent_0_string, agent_2_string)
        hierarchy_topology_plot_config2['settings']['store_colors'][new_id1] = color
        hierarchy_topology_plot_config2['settings']['store_colors'][new_id2] = color
    else:
        hierarchy_topology_plot_config2['settings']['store_colors'][node_id] = color
        
for node_id in hierarchy_topology_plot_config['settings']['remove_nodes']:
    if agent_0_string in node_id:
        new_id1 = node_id.replace(agent_0_string, agent_1_string)
        new_id2 = node_id.replace(agent_0_string, agent_2_string)
        hierarchy_topology_plot_config2['settings']['remove_nodes'].extend([new_id1, new_id2])



### Tx: a deterministic transcription process

class Tx(Process): # class for Tx process
    
    defaults = {
        'ktsc': 5e-2,
        'kdeg': 5e-3}

    def __init__(self, parameters=None): # superclass constructor fpr Tx process
        super().__init__(parameters)
        
    def ports_schema(self): # ports for Tx process
        return {
            'DNA': {  
                'G': {
                    '_default': 10 * units.mg / units.mL,         
                    '_updater': 'accumulate',
                    '_emit': True}},
            'mRNA': {
                'C': {
                    '_default': 100 * units.mg / units.mL,
                    '_updater': 'accumulate',
                    '_emit': True}}}
    
    def next_update(self, timestep, states): # update method for Tx process
        G = states['DNA']['G']
        C = states['mRNA']['C']
        dC = (self.parameters['ktsc'] * G - self.parameters['kdeg'] * C) * timestep
        return {
            'mRNA': {
                'C': dC}}

# plot topology of Tx process
fig = plot_topology(Tx(), filename='tx_topology.pdf', **topology_plot_config)

# tsc configuration
tx_config = {'time_step': 10}
tx_sim_settings = {
    'experiment_id': 'TX'}
tx_initial_state = {
    'DNA': {'G': 10.0 * units.mg/units.mL},
    'mRNA': {'C': 0.0 * units.mg/units.mL}}
tx_plot_config = {
    'variables': [
        {
            'variable': ('mRNA', ('C', 'milligram / milliliter')),
            'color': store_colors['mRNA']
        }, 
        {
            'variable': ('DNA', ('G', 'milligram / milliliter')),
            'color': store_colors['DNA']
        }],
    'filename': 'tx_output.pdf',
    **plot_var_config}

# initialize
tx_process = Tx(tx_config)

# make the experiment
tx_exp = process_in_experiment(
    tx_process, tx_sim_settings, tx_initial_state)

# run
tx_exp.update(10000)

# retrieve the data as a timeseries
tx_output = tx_exp.emitter.get_timeseries()

# plot
fig = plot_variables(tx_output,  **tx_plot_config)


### Tl: a deterministic translation process

class Tl(Process):

    defaults = {
        'ktrl': 5e-4,
        'kdeg': 5e-5}

    def ports_schema(self):
        return {
            'mRNA': {
                'C': {
                    '_default': 100 * units.mg / units.mL,
                    '_divider': 'split',
                    '_emit': True}},
            'Protein': {
                'X': {
                    '_default': 200 * units.mg / units.mL,
                    '_divider': 'split',
                    '_emit': True}}}

    def next_update(self, timestep, states):
        C = states['mRNA']['C']
        X = states['Protein']['X']
        dX = (self.parameters['ktrl'] * C - self.parameters['kdeg'] * X) * timestep
        return {
            'Protein': {
                'X': dX}}

# trl configuration
tl_config = {'time_step': 10}
tl_sim_settings = {'experiment_id': 'TL'}
tl_initial_state = {
    'Protein': {'X': 0.0 * units.mg / units.mL}}
tl_plot_config = {
    'variables': [
        {
            'variable': ('Protein', ('X', 'milligram / milliliter')),
            'color': store_colors['Protein']
        },
        {
            'variable': ('mRNA', ('C', 'milligram / milliliter')),
            'color': store_colors['mRNA']
        },  
        ],
    'filename': 'tl_output.pdf',
    **plot_var_config}

# initialize
tl_process = Tl(tl_config)

# make the experiment
tl_exp = process_in_experiment(
    tl_process, tl_sim_settings, tl_initial_state)

# run
tl_exp.update(10000)

# retrieve the data as a timeseries
tl_output = tl_exp.emitter.get_timeseries()

# plot
fig = plot_variables(tl_output,  **tl_plot_config)

