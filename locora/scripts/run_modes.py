from collections import OrderedDict

from run_grid_solvent import grid_solvent, grid_solvent_out
from run_process_data import process_data, process_data_out

### The following will register the run modes (i.e. the functions defined above)
### for run_locora executable.
### Each run mode key has a corresponding list, that contains [processing method, 
### writing output method, method can be run in parallel yes/no]
run_modes = OrderedDict()
run_modes['grid_solvent'] = [grid_solvent, grid_solvent_out, 1]
run_modes['process_data'] = [process_data, process_data_out, 0]