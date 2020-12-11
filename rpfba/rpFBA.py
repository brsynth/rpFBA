import logging
from cobra.flux_analysis import pfba
from brs_libs            import rpSBML

# logger = logger.getLogger(__name__)
# logger.setLevel(logger.INFO)
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)


# TODO: add the pareto frontier optimisation as an automatic way to calculate the optimal fluxes


# def runMultiObjective(rpsbml,
#                       reactions,
#                       coefficients,
#                       is_max=True,
#                       pathway_id='rp_pathway',
#                       objective_id=None):
#     """Run FBA using multiple objectives
#
#     :param reactions: The ids of the reactions involved in the objective
#     :param coefficients: The coefficients associated with the reactions id
#     :param is_max: Maximise or minimise the objective (Default: True)
#     :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
#     :param objective_id: Overwrite the default id (Default: None)
#
#     :type reactions: list
#     :type coefficients: list
#     :type is_max: bool
#     :type pathway_id: str
#     :type objective_id: str
#
#     :return: Success or failure of the function
#     :rtype: bool
#     """
#     fbc_plugin = rpsbml.getModel().getPlugin('fbc')
#     rpsbml.checklibSBML(fbc_plugin, 'Getting FBC package')
#     objective_id = rpsbml.findCreateObjective(reactions, coefficients, is_max)
#     rpsbml.checklibSBML(fbc_plugin.setActiveObjectiveId(objective_id),
#                         'Setting active objective '+str(objective_id))
#     cobraModel = rpsbml.convertToCobra()
#     if not cobraModel:
#         return False
#     cobra_results = cobraModel.optimize()
#     rpsbml.addAnalysisResults(objective_id, cobra_results, pathway_id)
#     return rpsbml


def rp_fba(rpsbml,
           reaction_id,
           coefficient=1.0,
           is_max=True,
           pathway_id='rp_pathway',
           objective_id=None,
           logger=None):
    """Run FBA using a single objective

    :param reaction_id: The id of the reactions involved in the objective
    :param coefficient: The coefficient associated with the reactions id (Default: 1.0)
    :param is_max: Maximise or minimise the objective (Default: True)
    :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
    :param objective_id: Overwrite the default id (Default: None)

    :type reaction_id: str
    :type coefficient: float
    :type is_max: bool
    :type pathway_id: str
    :type objective_id: str

    :return: Tuple with the results of the FBA and boolean indicating the success or failure of the function
    :rtype: tuple
    """
    logger = logger or logging.getLogger(__name__)
    fbc_plugin = rpsbml.getModel().getPlugin('fbc')
    rpsbml.checklibSBML(fbc_plugin, 'Getting FBC package')
    objective_id = rpsbml.findCreateObjective([reaction_id], [coefficient], is_max, objective_id)
    # run the FBA
    rpsbml.checklibSBML(fbc_plugin.setActiveObjectiveId(objective_id),
                        'Setting active objective '+str(objective_id))
    cobraModel = rpsbml.convertToCobra()
    if not cobraModel:
        return 0.0, False
    cobra_results = cobraModel.optimize()
    rpsbml.addAnalysisResults(objective_id, cobra_results, pathway_id)
    return cobra_results.objective_value, rpsbml


def rp_pfba(rpsbml,
            reaction_id,
            coefficient=1.0,
            fraction_of_optimum=0.95,
            is_max=True,
            pathway_id='rp_pathway',
            objective_id=None,
            logger=None):
    """Run parsimonious FBA using a single objective

    :param reaction_id: The id of the reactions involved in the objective
    :param coefficient: The coefficient associated with the reactions id (Default: 1.0)
    :param fraction_of_optimum: Between 0.0 and 1.0 determining the fraction of optimum (Default: 0.95)
    :param is_max: Maximise or minimise the objective (Default: True)
    :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
    :param objective_id: Overwrite the default id (Default: None)

    :type reaction_id: str
    :type coefficient: float
    :type fraction_of_optimum: float
    :type is_max: bool
    :type pathway_id: str
    :type objective_id: str

    :return: Tuple with the results of the FBA and boolean indicating the success or failure of the function
    :rtype: tuple
    """
    logger = logger or logging.getLogger(__name__)
    fbc_plugin = rpsbml.getModel().getPlugin('fbc')
    rpsbml.checklibSBML(fbc_plugin, 'Getting FBC package')
    objective_id = rpsbml.findCreateObjective([reaction_id], [coefficient], is_max, objective_id)
    # run the FBA
    rpsbml.checklibSBML(fbc_plugin.setActiveObjectiveId(objective_id),
                        'Setting active objective '+str(objective_id))
    cobraModel = rpsbml.convertToCobra()
    if not cobraModel:
        return 0.0, False
    cobra_results = pfba(cobraModel, fraction_of_optimum)
    rpsbml.addAnalysisResults(objective_id, cobra_results, pathway_id)
    return cobra_results.objective_value, rpsbml


def rp_fraction(rpsbml,
                source_reaction,
                source_coefficient,
                target_reaction,
                target_coefficient,
                fraction_of_source=0.75,
                is_max=True,
                pathway_id='rp_pathway',
                objective_id=None,
                logger=None):
    """Optimise for a target reaction while fixing a source reaction to the fraction of its optimum

    :param source_reaction: The id of the source reaction
    :param source_coefficient: The source coefficient associated with the source reaction id
    :param target_reaction: The id of the target reaction
    :param target_coefficient: The source coefficient associated with the target reaction id
    :param fraction_of_optimum: Between 0.0 and 1.0 determining the fraction of optimum (Default: 0.75)
    :param is_max: Maximise or minimise the objective (Default: True)
    :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
    :param objective_id: Overwrite the default id (Default: None)

    :type source_reaction: str
    :type source_coefficient: float
    :type target_reaction: str
    :type target_coefficient: float
    :type fraction_of_optimum: float
    :type is_max: bool
    :type pathway_id: str
    :type objective_id: str

    :return: Tuple with the results of the FBA and boolean indicating the success or failure of the function
    :rtype: tuple
    """
    logger = logger or logging.getLogger(__name__)
    # retreive the biomass objective and flux results and set as maxima
    fbc_plugin = rpsbml.getModel().getPlugin('fbc')
    rpsbml.checklibSBML(fbc_plugin, 'Getting FBC package')
    logger.debug('findCreateObjective: '+str(source_reaction))
    source_obj_id = rpsbml.findCreateObjective([source_reaction], [source_coefficient], is_max)
    # TODO: use the rpSBML BRSynth annotation parser
    source_flux = None
    try:
        fbc_obj = fbc_plugin.getObjective(source_obj_id)
        # TODO: if this is None need to set it up
        fbc_obj_annot = fbc_obj.getAnnotation()
        if not fbc_obj_annot:
            raise ValueError
        source_flux = float(fbc_obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(0).getAttrValue('value'))
        logger.debug('Already calculated flux for '+str(source_obj_id))
    except (AttributeError, ValueError) as e:
        logger.debug(e)
        logger.debug('Performing FBA to calculate the source reaction')
        ### FBA ###
        # self.runFBA(source_reaction, source_coefficient, is_max, pathway_id)
        rpsbml.checklibSBML(fbc_plugin.setActiveObjectiveId(source_obj_id),
                            'Setting active objective '+str(source_obj_id))
        cobraModel = rpsbml.convertToCobra()
        if not cobraModel:
            logger.error('Converting libSBML to CobraPy returned False')
            rpsbml.addAnalysisResults(source_obj_id, 0.0, pathway_id)
            return 0.0, False
        cobra_results = cobraModel.optimize()
        rpsbml.addAnalysisResults(source_obj_id, cobra_results, pathway_id)
        # cobra_results.objective_value
        fbc_obj = fbc_plugin.getObjective(source_obj_id)
        fbc_obj_annot = fbc_obj.getAnnotation()
        if fbc_obj_annot is None:
            logger.error('No annotation available for: '+str(source_obj_id))
        source_flux = float(fbc_obj_annot.getChild('RDF').getChild('BRSynth').getChild('brsynth').getChild(0).getAttrValue('value'))
    # TODO: add another to check if the objective id exists
    logger.debug('FBA source flux ('+str(source_reaction)+') is: '+str(source_flux))
    if not objective_id:
        objective_id = 'obj_'+str(target_reaction)+'__restricted_'+str(source_reaction)
    logger.debug('findCreateObjective() for '+str(objective_id))
    objective_id = rpsbml.findCreateObjective([target_reaction], [target_coefficient], is_max, objective_id)
    logger.debug('Optimising the objective: '+str(objective_id))
    logger.debug('Setting upper bound: '+str(source_flux*fraction_of_source))
    logger.debug('Setting loer bound: '+str(source_flux*fraction_of_source))
    old_upper_bound, old_lower_bound = rpsbml.setReactionConstraints(source_reaction,
                                                                     source_flux*fraction_of_source,
                                                                     source_flux*fraction_of_source)
    rpsbml.checklibSBML(fbc_plugin.setActiveObjectiveId(objective_id),
                        'Setting active objective '+str(objective_id))
    cobraModel = rpsbml.convertToCobra()
    if not cobraModel:
        # although this may not be the greatest idea, set flux to 0.0 when cobrapy error
        rpsbml.addAnalysisResults(objective_id, 0.0, pathway_id)
        return 0.0, False
    cobra_results = cobraModel.optimize()
    rpsbml.addAnalysisResults(objective_id, cobra_results, pathway_id)
    ##### print the biomass results ######
    logger.debug('Biomass: '+str(cobra_results.fluxes.biomass))
    logger.debug('Target: '+str(cobra_results.fluxes.RP1_sink))
    # reset the bounds to the original values for the target
    old_upper_bound, old_lower_bound = rpsbml.setReactionConstraints(source_reaction,
                                                                     old_upper_bound,
                                                                     old_lower_bound)
    logger.debug('The objective '+str(objective_id)+' results '+str(cobra_results.objective_value))
    return cobra_results.objective_value, rpsbml


# TODO: do not use the species_group_id and the sink_species_group_id. Loop through all the groups (and if the same) and overwrite the annotation instead
def runFBA(sbml_path, gem_sbml, outFile,
           sim_type,
           source_reaction, target_reaction,
           source_coefficient, target_coefficient,
           is_max,
           fraction_of,
           dont_merge=True,
           pathway_id='rp_pathway',
           objective_id=None,
           compartment_id='MNXC3',
           # fill_orphan_species=False,
           species_group_id='central_species',
           sink_species_group_id='rp_sink_species',
           logger=None):
    """Single rpSBML simulation

    :param file_name: The name of the model
    :param sbml_path: Path to the rpSBML file
    :param gem_sbml: Path to the GEM file
    :param sim_type: The type of simulation to use. Available simulation types include: fraction, fba, rpfba
    :param source_reaction: The reaction id of the source reaction.
    :param target_reaction: The reaction id of the target reaction. Note that if fba or rpfba options are used, then these are ignored
    :param source_coefficient: The source coefficient
    :param target_coefficient: The target coefficient
    :param is_max: Maximise or minimise the objective
    :param fraction_of: The fraction of the optimum. Note that this value is ignored is fba is used
    :param tmpOutputFolder: The path to the output document
    :param dont_merge: Output the merged model (Default: True)
    :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
    :param objective_id: Overwrite the auto-generated id of the results (Default: None)
    :param compartment_id: The SBML compartment id (Default: MNXC3)
    :param fill_orphan_species: Add pseudo reactions that consume/produce single parent species. Note in development
    :param species_group_id: The id of the central species (Default: central_species)
    :param sink_species_group_id: The id of the sink species (Default: rp_sink_species)

    :type inputTar: str
    :type gem_sbml: str
    :type sim_type: str
    :type source_reaction: str
    :type target_reaction: str
    :type source_coefficient: float
    :type target_coefficient: float
    :type is_max: bool
    :type fraction_of: float
    :type tmpOutputFolder: str
    :type dont_merge: bool
    :type num_workers: int
    :type pathway_id: str
    :type objective_id: str
    :type compartment_id: str
    :type fill_orphan_species: bool
    :type species_group_id: str
    :type sink_species_group_id: str

    :return: Succcess or failure of the function
    :rtype: bool
    """

    logger = logger or logging.getLogger(__name__)

    logger.info('--------- '+str(outFile)+' ------------')

    rpsbml = rpSBML(sbml_path, logger=logger)
    # Save the central species
    groups = rpsbml.getModel().getPlugin('groups')
    central = groups.getGroup(species_group_id)
    sink_group = groups.getGroup(sink_species_group_id)
    rp_group = groups.getGroup(pathway_id)
    cent_spe = [str(i.getIdRef()) for i in central.getListOfMembers()]
    sink_spe = [str(i.getIdRef()) for i in sink_group.getListOfMembers()]
    rp_reac = [str(i.getIdRef()) for i in rp_group.getListOfMembers()]
    logger.debug('old central species: '+str(cent_spe))
    logger.debug('old sink species: '+str(sink_spe))
    logger.debug('old rp reactions: '+str(rp_reac))

    rpsbml_gem = rpSBML(gem_sbml, logger=logger)
    species_source_target, reactions_convert = rpSBML.mergeModels(rpsbml, rpsbml_gem, logger)
    # NOTE: reactions_convert is organised with key being the rpsbml reaction and value being the rpsbml_gem value`
    # BUG: when merging the RP1_sink (very rare cases) can be recognised if another reaction contains the same species as a reactant
    ## under such as scenario the algorithm will consider that they are the same -- TODO: overwrite it
    if target_reaction in reactions_convert:
        logger.warning('The target_reaction ('+str(target_reaction)+') has been detected in model '+str(outFile)+', ignoring this model...')
        return False
    rev_reactions_convert = {v: k for k, v in reactions_convert.items()}
    logger.debug('species_source_target: '+str(species_source_target))
    logger.debug('reactions_convert: '+str(reactions_convert))
    logger.debug('rev_reactions_convert: '+str(rev_reactions_convert))

    # TO TEST MERGE: TO REMOVE
    # rpsbml_gem.modelName = 'test'
    # rpsbml_gem.writeSBML('/home/mdulac/workspace/Galaxy-SynBioCAD/rpFBA/rpFBA_image/tmp_out/')
    ####### fraction of reaction ######
    # Choose the good function according to 'sim_type' argument
    sim_func = globals()['rp_'+sim_type]
    obj_val,rpsbml_gem = sim_func(rpsbml_gem,
                                  source_reaction, source_coefficient,
                                  target_reaction, target_coefficient,
                                  fraction_of,
                                  is_max,
                                  pathway_id, objective_id,
                                  logger)
    # ####### FBA ########
    # elif sim_type=='fba':
    #     obj_val,rpsbml_gem = fba(rpsbml_gem,
    #                                   source_reaction, source_coefficient,
    #                                   target_reaction, target_coefficient,
    #                                   fraction_of,
    #                                   is_max,
    #                                   pathway_id, objective_id)
    # ####### pFBA #######
    # elif sim_type=='pfba':
    #     obj_val,rpsbml_gem = pfba(rpsbml_gem,
    #                                   source_reaction, source_coefficient,
    #                                   target_reaction, target_coefficient,
    #                                   fraction_of,
    #                                   is_max,
    #                                   pathway_id, objective_id)
    # else:
    #     logger.error('Cannot recognise sim_type: '+str(sim_type))
    '''
    ###### multi objective #####
    elif sim_type=='multi_fba':
        rpfba.runMultiObjective(reactions, coefficients, is_max, pathway_id)
    '''
    if dont_merge:
        logger.info('Returning model with heterologous pathway only')
        groups = rpsbml_gem.getModel().getPlugin('groups')
        rp_pathway = groups.getGroup(pathway_id)
        logger.debug('---- Reactions ----')
        for member in rp_pathway.getListOfMembers():
            #### reaction annotation
            logger.debug(member.getIdRef())
            reacFBA = rpsbml_gem.getModel().getReaction(member.getIdRef())
            logger.debug(reacFBA)
            try:
                #reacIN = rpsbml.model.getReaction(reactions_convert[member.getIdRef()])
                reacIN = rpsbml.getModel().getReaction(rev_reactions_convert[member.getIdRef()])
            except KeyError:
                reacIN = rpsbml.getModel().getReaction(member.getIdRef())
            logger.debug(reacIN)
            logger.debug(reacFBA.getAnnotation())
            reacIN.setAnnotation(reacFBA.getAnnotation())
            #### species TODO: only for shadow price
        #### add groups ####
        source_groups = rpsbml_gem.getModel().getPlugin('groups')
        target_groups = rpsbml.getModel().getPlugin('groups')
        target_groupsID = [i.getId() for i in target_groups.getListOfGroups()]
        for source_group in source_groups.getListOfGroups():
            logger.info('Replacing group id: '+str(source_group.getId()))
            if source_group.getId()==species_group_id:
                target_group = target_groups.getGroup(source_group.getId())
                # TODO: #### replace the new potentially incorect central species with the normal ones #####
                # delete all the previous members
                logger.info('Removing central_species')
                for i in range(target_group.getNumMembers()):
                    logger.debug('Deleting group member: '+str(target_group.getMember(0).getIdRef()))
                    target_group.removeMember(0)
                # add the new ones
                for cs in cent_spe:
                    logger.info('Creating new member: '+str(cs))
                    newM = target_group.createMember()
                    newM.setIdRef(cs)
            elif source_group.getId()==sink_species_group_id:
                target_group = target_groups.getGroup(source_group.getId())
                logger.info('Removing sink species')
                for i in range(target_group.getNumMembers()):
                    logger.info('Deleting group member: '+str(target_group.getMember(0).getIdRef()))
                    target_group.removeMember(0)
                #add the new ones
                for cs in sink_spe:
                    logger.info('Creating new member: '+str(cs))
                    newM = target_group.createMember()
                    newM.setIdRef(cs)
            elif source_group.getId() in target_groupsID:
                target_group = target_groups.getGroup(source_group.getId())
                target_group.setAnnotation(source_group.getAnnotation())
            '''
            elif source_group.getId()==pathway_id:
                target_group = target_groups.getGroup(source_group.getId())
                logger.debug('Removing rp ractions')
                for i in range(target_group.getNumMembers()):
                    logger.debug('Deleting group member: '+str(target_group.getMember(0).getIdRef()))
                    target_group.removeMember(0)
                #add the new ones
                for cs in rp_reac:
                    logger.debug('Creating new member: '+str(cs))
                    newM = target_group.createMember()
                    newM.setIdRef(cs)
            '''
        #### add objectives ####
        source_fbc = rpsbml_gem.getModel().getPlugin('fbc')
        target_fbc = rpsbml.getModel().getPlugin('fbc')
        target_objID = [i.getId() for i in target_fbc.getListOfObjectives()]
        for source_obj in source_fbc.getListOfObjectives():
            source_obj_id = source_obj.getId()
            if source_obj.getId() in target_objID:
                target_obj = target_fbc.getObjective(source_obj.getId())
                target_obj.setAnnotation(source_obj.getAnnotation())
                for target_fluxObj in target_obj.getListOfFluxObjectives():
                    for source_fluxObj in source_obj.getListOfFluxObjectives():
                        if target_fluxObj.getReaction()==source_fluxObj.getReaction():
                            target_fluxObj.setAnnotation(source_fluxObj.getAnnotation())
            else:
                target_fbc.addObjective(source_obj)
        #rpsbml.createMultiFluxObj('obj_RP1_sink', ['RP1_sink'], [1])
        target_fbc.setActiveObjectiveId(source_obj_id) #tmp random assigenement of objective
        rpsbml.writeSBML(outFile)
    else:
        logger.info('Returning the full model')
        rpsbml_gem.writeSBML(outFile)
    return True


########################################################################
############################### FBA pathway ranking ####################
########################################################################

# 1) Number of interventions
# need to calculate the number of steps that are not native to know the number of interventions

# 2) Maximal growth rate

# 3) Minimum product yeild at maximal growth rate

# 4) Minimum product yeild

# 5) Anaerobic condition

# 6) Number of potentially disruptive products

    # Toxicity?

# 7) Number of accessible metabolites (avoid intermediate accumulation)

# 8) Thermodynamics (MDF)

# 9) The overlap of the same changes --> might not be applicable in our case

# 10) Reduced model

# 11) ECM
