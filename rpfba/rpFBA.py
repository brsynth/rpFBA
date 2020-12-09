from cobra.flux_analysis import pfba
from logging             import getLogger

logger = getLogger(__name__)


# TODO: add the pareto frontier optimisation as an automatic way to calculate the optimal fluxes


def runMultiObjective(rpsbml,
                      reactions,
                      coefficients,
                      is_max=True,
                      pathway_id='rp_pathway',
                      objective_id=None):
    """Run FBA using multiple objectives

    :param reactions: The ids of the reactions involved in the objective
    :param coefficients: The coefficients associated with the reactions id
    :param is_max: Maximise or minimise the objective (Default: True)
    :param pathway_id: The id of the heterologous pathway (Default: rp_pathway)
    :param objective_id: Overwrite the default id (Default: None)

    :type reactions: list
    :type coefficients: list
    :type is_max: bool
    :type pathway_id: str
    :type objective_id: str

    :return: Success or failure of the function
    :rtype: bool
    """
    fbc_plugin = rpsbml.getModel().getPlugin('fbc')
    rpsbml.checklibSBML(fbc_plugin, 'Getting FBC package')
    objective_id = rpsbml.findCreateObjective(reactions, coefficients, is_max)
    rpsbml.checklibSBML(fbc_plugin.setActiveObjectiveId(objective_id),
                        'Setting active objective '+str(objective_id))
    cobraModel = rpsbml.convertToCobra()
    if not cobraModel:
        return False
    cobra_results = cobraModel.optimize()
    rpsbml.addAnalysisResults(objective_id, cobra_results, pathway_id)
    return rpsbml


def runFBA(rpsbml,
           reaction_id,
           coefficient=1.0,
           is_max=True,
           pathway_id='rp_pathway',
           objective_id=None):
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


def runParsimoniousFBA(rpsbml,
                       reaction_id,
                       coefficient=1.0,
                       fraction_of_optimum=0.95,
                       is_max=True,
                       pathway_id='rp_pathway',
                       objective_id=None):
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


def runFractionReaction(rpsbml,
                        source_reaction,
                        source_coefficient,
                        target_reaction,
                        target_coefficient,
                        fraction_of_source=0.75,
                        is_max=True,
                        pathway_id='rp_pathway',
                        objective_id=None):
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
