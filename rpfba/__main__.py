#!/usr/bin/env python

import logging
from rpfba import runFBA, build_args_parser


def _cli():
    parser = build_args_parser()
    args  = parser.parse_args()

    # Create logger
    logger = logging.getLogger(__name__)
    logger.setLevel(getattr(logging, args.log.upper()))
    logger.formatter = logging.Formatter('%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s')
    # logging.basicConfig(
    #     level=getattr(logging, args.log.upper()),
    #     format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    #     datefmt='%d-%m-%Y %H:%M:%S',
    # )

    result = runFBA(args.input_sbml, args.gem_sbml, args.outfile,
                    args.sim_type,
                    args.source_reaction, args.target_reaction,
                    args.source_coefficient, args.target_coefficient,
                    args.is_max,
                    args.fraction_of,
                    args.dont_merge,
                    args.pathway_id,
                    args.objective_id,
                    args.compartment_id,
                    args.species_group_id,
                    args.sink_species_group_id,
                    logger)

    return result



if __name__ == '__main__':
    _cli()
