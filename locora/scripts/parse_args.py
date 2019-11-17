import argparse
from locora.scripts.run_modes import run_modes

def parse_args():

    parser = argparse.ArgumentParser(
            description="Executable for solute-solvent local correlation analysis.")

    ### Required Arguments
    ### ~~~~~~~~~~~~~~~~~~

    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input', 
                          required=True, 
                          type=str, 
                          default=None,
                          help="Input file.")

    parser.add_argument('-m', '--mode',
                        required=True,
                        type=str,
                        default=None,
                        help="Run mode.",
                        choices=run_modes.keys())

    ### Optional Arguments ###
    ### ~~~~~~~~~~~~~~~~~~ ###

    parser._action_groups.append(parser._action_groups.pop(1))

    process_data_info = "Only used in mode=process_data."

    parser.add_argument('-noim', '--noimage',
                        action='store_true',
                        help='Swtich off re-imaging during calculations.')

    parser.add_argument('-c', '--cutoff',
                        required=False,
                        type=float,
                        default=3.5,
                        help='Use only water molecules that within cutoff (in Ang.),\
                              around the center of the fractional coordinate system. \
                              Default is 3.5 Ang. %s' %process_data_info)

    parser.add_argument('-ts', '--timestep',
                        required=False,
                        type=float,
                        default=1.,
                        help='Timestep per MD frame in ps. Default is 1 ps. %s' %process_data_info)

    parser.add_argument('-ms', '--minstep',
                        required=False,
                        type=int,
                        default=1,
                        help='Minium number of occurancies per window, in order to analyse lifetime \
                        distribution. Default is 1. %s' %process_data_info)

    parser.add_argument('-tr', '--transient',
                        required=False,
                        type=int,
                        nargs=2,
                        default=[0, 1],
                        help='Scan transient re-entering time of water molecules in multiples \
                              of timestep. Start with -tr[0] and end with -tr[1]. \
                              Default is [0, 1]. %s' %process_data_info)

    parser.add_argument('-np', '--nproc',
                        required=False,
                        type=int,
                        default=1,
                        help='Total number of procesess. Default is 1(=no multiprocessing).')

    parser.add_argument('-pre', '--prefix',
                        required=False,
                        type=str,
                        default="",
                        help='Output prefix. Default is \'\'. %s' %process_data_info)

    parser.add_argument('-w', '--window',
                        required=False,
                        type=int,
                        default=1000,
                        help='window size used for calculation of time correlation.\
                              Default is 1000. %s' %process_data_info)

    parser.add_argument('-b', '--bootstrap',
                        required=False,
                        type=int,
                        default=0,
                        help='Number of bootstrapping resample steps.\
                              If bootstrap=0, no bootstrapping will be carried out and \
                              averages/standard deviation are calculated from subsequent \
                              windows of length --window. If bootstrap>0, --boostrap bootstrapping \
                              iterations are performed for windows of legnth --window. \
                              Default is 0. %s' %process_data_info)

    parser.add_argument('-pl', '--plot',
                        action='store_true',
                        help='Turn on plotting. Default is off (=no plotting). %s' %process_data_info)

    parser.add_argument('-ai', '--anisotropic',
                        action='store_true',
                        help='Perform anisotropic analysis. Default is off (=no anisotropic analysis). %s' %process_data_info)

    parser.add_argument('-lp', '--legendre',
                        required=False,
                        type=int,
                        default=1,
                        help='Legendre polynomial used for calculation of orientational lifetimes. \
                              Can be either 1 or 2. Default is 1(=First order Legendre polynomial). %s' %process_data_info)

    parser.add_argument('-de', '--double-exp',
                        action='store_true',
                        help='Perform double exponential fitting. %s' %process_data_info)

    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='Verbosity output. Default is off.')

    return parser.parse_args()