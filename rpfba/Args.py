from argparse  import ArgumentParser as argparse_ArgumentParser


def build_args_parser():
    parser = argparse_ArgumentParser(prog='rpfba', description='Calculate FBA to generate rpFBA collection')
    parser = _add_arguments(parser)
    return parser


def _add_arguments(parser):
    parser.add_argument('input_sbml',
                        type=str,
                        help='input SBML file')
    parser.add_argument('gem_sbml',
                        type=str,
                        help='GEM file')
    parser.add_argument('outfile',
                        type=str,
                        help='output file')
    parser.add_argument('--pathway_id',
                        type=str,
                        default='rp_pathway',
                        help='id of the heterologous pathway (default: rp_pathway)')
    parser.add_argument('--sink_species_group_id',
                        type=str,
                        default='rp_sink_species',
                        help='id of the central species (default: central_species)')
    parser.add_argument('--species_group_id',
                        type=str,
                        default='central_species',
                        help='id of the sink species (default: rp_sink_species)')
    parser.add_argument('--objective_id',
                        type=str,
                        default='None',
                        help='overwrite the auto-generated id of the results (default: None)')
    parser.add_argument('--compartment_id',
                        type=str,
                        default='MNXC3',
                        help='SBML compartment id (default: MNXC3)')
    parser.add_argument('--sim_type',
                        type=str,
                        default='fraction',
                        help='type of simulation to use. Available simulation types include: \'fraction_reaction\', \'fba\', \'rpfba\'')
    parser.add_argument('--source_reaction',
                        type=str,
                        default='biomass',
                        help='reaction id of the source reaction')
    parser.add_argument('--target_reaction',
                        type=str,
                        default='RP1_sink',
                        help='reaction id of the target reaction. Note that if \'fba\' or \'rpfba\' options are used, then these are ignored')
    parser.add_argument('--source_coefficient',
                        type=float,
                        default=1.0,
                        help='source coefficient')
    parser.add_argument('--target_coefficient',
                        type=float,
                        default=1.0,
                        help='target coefficient')
    # parser.add_argument('--num_workers',
    #                     type=int,
    #                     default=10,
    #                     help='number of workers (multi-threads)')
    parser.add_argument('--is_max',
                        action='store_true',
                        help='maximise the objective (default)')
    parser.add_argument('--fraction_of',
                        type=float,
                        default=0.75,
                        help='fraction of the optimum. Note that this value is ignored is \'fba\' is used')
    parser.add_argument('--dont_merge',
                        action='store_true',
                        help='output the merged model (default)')
    parser.add_argument('--log', metavar='ARG',
                        type=str, choices=['debug', 'info', 'warning', 'error', 'critical'],
                        default='error',
                        help='Adds a console logger for the specified level (default: error)')
    return parser
