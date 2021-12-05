import pandas as pd

from scipy import constants

import cobra
import cobra.test
from cobra.medium import minimal_medium

from cobra import Model, Reaction, Metabolite, Configuration
from vivarium.library.units import units

AVOGADRO = constants.N_A * 1 / units.mol

EXTERNAL_PREFIX = 'EX_'
DEFAULT_UPPER_BOUND = 100.0
DEFAULT_TOLERANCE = [0.95, 1]


def build_model(
        stoichiometry,
        reversible,
        objective,
        external_molecules,
        default_upper_bound=DEFAULT_UPPER_BOUND
):
    model = Model('fba')
    model.compartments = {'c': 'cytoplasm'}

    metabolite_keys = {}
    for reaction_key, chemistry in stoichiometry.items():
        metabolite_keys.update(chemistry)

    metabolites = {
        metabolite: Metabolite(metabolite, name=metabolite, compartment='c')
        for metabolite in list(metabolite_keys.keys())}

    model.add_metabolites(metabolites.values())

    # make reactions
    reactions = {}
    for reaction_key, chemistry in stoichiometry.items():
        reaction = Reaction(reaction_key, name=reaction_key)

        # set reaction bounds
        reaction.upper_bound = default_upper_bound
        if reaction_key in reversible:
            reaction.lower_bound = -reaction.upper_bound

        # make stoichiometry
        reaction_model = {
            metabolites[metabolite]: value
            for metabolite, value in chemistry.items()}
        reaction.add_metabolites(reaction_model)

        reactions[reaction_key] = reaction

    # make exchange reactions for all external_molecules
    for external in external_molecules:
        external_key = EXTERNAL_PREFIX + external
        reaction = Reaction(external_key, name=external_key)

        # set reaction bounds
        reaction.upper_bound = default_upper_bound
        reaction.lower_bound = -default_upper_bound  # TODO -- should exchanges have symmetric bounds by default?

        # make stoichiometry
        reaction_model = {metabolites[external]: -1}
        reaction.add_metabolites(reaction_model)

        reactions[external_key] = reaction

    model.add_reactions(reactions.values())

    model.objective = {
        reactions[reaction_key]: value
        for reaction_key, value in objective.items()}

    return model

def extract_model(model):
    """
    TODO -- where do demands and sinks go?
    demands = model.demands
    sinks = model.sinks

    # boundary reactions include exchanges, demands, sinks
    boundary = model.boundary
    boundary_reactions = [reaction.id for reaction in boundary]
    """

    reactions = model.reactions
    metabolites = model.metabolites
    exchanges = model.exchanges
    objective_expression = model.objective.expression.args

    # get stoichiometry and flux bounds
    stoichiometry = {}
    flux_bounds = {}
    reversible = []
    for reaction in reactions:
        reaction_id = reaction.id
        reaction_metabolites = reaction.metabolites
        bounds = list(reaction.bounds)
        stoichiometry[reaction.id] = {
            metabolite.id: coeff for metabolite, coeff in reaction_metabolites.items()}
        flux_bounds[reaction_id] = bounds
        if not any(b == 0.0 for b in bounds):
            reversible.append(reaction_id)

    # get external molecules and exchange bounds from exchanges
    external_molecules = []
    exchange_bounds = {}
    for reaction in exchanges:
        reaction_metabolites = list(reaction.metabolites.keys())
        assert len(reaction_metabolites) == 1  # only 1 molecule in the exchange reaction
        metabolite_id = reaction_metabolites[0].id
        external_molecules.append(metabolite_id)
        bounds = list(reaction.bounds)
        exchange_bounds[metabolite_id] = bounds

    # get molecular weights
    molecular_weights = {}
    nonexisting_elements = ['R', 'X']
    for metabolite in metabolites:
        if any([e in nonexisting_elements for e in metabolite.elements]):
            continue
        molecular_weights[metabolite.id] = metabolite.formula_weight

    # get objective
    objective = {}
    for expression in objective_expression:
        exp_str = str(expression)
        coeff, reaction_id = exp_str.split('*')
        try:
            reactions.get_by_id(reaction_id)
            objective[reaction_id] = float(coeff)
        except:
            pass

    return {
        'stoichiometry': stoichiometry,
        'reversible': reversible,
        'external_molecules': external_molecules,
        'objective': objective,
        'flux_bounds': flux_bounds,
        'exchange_bounds': exchange_bounds,
        'molecular_weights': molecular_weights,
    }


class FBA(object):
    """
    This class provides an interface to cobra FBA.
    It can load in BiGG models (http://bigg.ucsd.edu/models) if provided a model_path to a saved JSON BiGG model,
    or load in a novel model specified by stoichiometry, reversibility, and objective.

    TODO (Eran) -- MOMA (Minimization Of Metabolic Adjustment) option is provided, but has not yet been tested.
    """

    cobra_configuration = Configuration()

    def __init__(self, config={}):
        model_path = config.get('model_path')

        # get tolerances
        self.default_tolerance = config.get('default_tolerance', DEFAULT_TOLERANCE)
        self.tolerance = config.get('tolerance', {})

        # set MOMA
        self.moma = config.get('moma', False)

        # initialize the model
        if model_path:
            # load a BiGG model
            self.model = cobra.io.load_json_model(model_path)
            extract = extract_model(self.model)

            self.stoichiometry = extract['stoichiometry']
            self.reversible = extract['reversible']
            self.external_molecules = extract['external_molecules']
            self.objective = extract['objective']
            self.flux_bounds = extract['flux_bounds']
            self.molecular_weights = extract['molecular_weights']
            self.exchange_bounds = extract['exchange_bounds']
            self.default_upper_bound = DEFAULT_UPPER_BOUND  # TODO -- can this be extracted from model?

        else:
            # create an FBA model from config
            self.stoichiometry = config['stoichiometry']
            self.reversible = config.get('reversible', [])
            self.external_molecules = config['external_molecules']
            self.objective = config['objective']
            self.flux_bounds = config.get('flux_bounds', {})
            self.molecular_weights = config.get('molecular_weights', {})
            self.exchange_bounds = config.get('exchange_bounds', {})
            self.default_upper_bound = config.get('default_upper_bound', DEFAULT_UPPER_BOUND)

            self.model = build_model(
                self.stoichiometry,
                self.reversible,
                self.objective,
                self.external_molecules,
                self.default_upper_bound)

            # apply constraints
            self.constrain_reaction_bounds(self.flux_bounds)
            self.set_exchange_bounds()

        self.exchange_bounds_keys = list(self.exchange_bounds.keys())

        # initialize solution
        self.solution = self.model.optimize()

        # set scaling factors on fluxes
        self.flux_scaling = 1
        if config.get('target_added_mass') is not None:
            target_added_mass = config['target_added_mass']
            self.flux_scaling_to_target_added(target_added_mass)

    def get_objective_composition(self):
        objective_composition = {}
        for reaction_id, coeff1 in self.objective.items():
            for mol_id, coeff2 in self.stoichiometry[reaction_id].items():
                if mol_id in objective_composition:
                    objective_composition[mol_id] += coeff1 * coeff2
                else:
                    objective_composition[mol_id] = coeff1 * coeff2
        return objective_composition

    def flux_scaling_to_target_added(self, target_added_mass):
        added_mass = self.get_added_mass()
        self.flux_scaling = target_added_mass / added_mass

    def get_added_mass(self, timestep=1):
        solution = self.model.optimize()
        objective_value = solution.objective_value * timestep
        added_mass = 0
        for reaction_id, coeff1 in self.objective.items():
            for mol_id, coeff2 in self.stoichiometry[reaction_id].items():
                if coeff2 < 0:
                    # molecule is used to make biomass (negative coefficient)
                    mw = self.molecular_weights[mol_id] * (units.g / units.mol)  # TODO -- clean up units
                    count = -coeff1 * coeff2 * objective_value
                    mol = count / AVOGADRO
                    mol_added_mass = mw * mol
                    added_mass += mol_added_mass.to('fg').magnitude

        return added_mass

    def minimal_external(self):
        '''get minimal external state'''
        max_growth = self.model.slim_optimize()
        max_exchange = minimal_medium(self.model, max_growth)
        return {ex[len(EXTERNAL_PREFIX):len(ex)]: value
            for ex, value in max_exchange.items()}

    def set_exchange_bounds(self, bounds={}):
        '''
        apply new_bounds for the defined molecules.
        reset unincluded molecules to their exchange_bounds.
        '''
        for external_mol, level in self.exchange_bounds.items():
            reaction = self.model.reactions.get_by_id(EXTERNAL_PREFIX + external_mol)

            if external_mol in bounds:
                level = bounds[external_mol] / self.flux_scaling

            if type(level) is list:
                reaction.upper_bound = level[1]
                reaction.lower_bound = level[0]
            elif isinstance(level, int) or isinstance(level, float):
                # reaction.upper_bound = level
                reaction.lower_bound = level

    def constrain_flux(self, bounds={}):
        '''add externally imposed constraints'''
        for reaction_id, bound in bounds.items():
            reaction = self.model.reactions.get_by_id(reaction_id)
            scaled_level = bound / self.flux_scaling

            if EXTERNAL_PREFIX in reaction_id and \
                any(substring in reaction_id for substring in self.exchange_bounds_keys):
                # exchanges use reverse flux
                scaled_level *= -1

            if reaction_id in self.tolerance:
                # use configured tolerance
                lower_tolerance, upper_tolerance = self.tolerance[reaction_id]
                reaction.upper_bound = upper_tolerance * scaled_level
                reaction.lower_bound = lower_tolerance * scaled_level
            else:
                # use default tolerance
                if bound >= 0:
                    reaction.upper_bound = self.default_tolerance[1] * scaled_level
                    reaction.lower_bound = self.default_tolerance[0] * scaled_level
                else:
                    reaction.upper_bound = self.default_tolerance[0] * scaled_level
                    reaction.lower_bound = self.default_tolerance[1] * scaled_level

    def constrain_reaction_bounds(self, bounds={}):
        reactions = self.get_reactions(list(bounds.keys()))
        for reaction_id, bound in bounds.items():
            reaction = reactions[reaction_id]
            scaled_bound = [b / self.flux_scaling for b in bound]
            reaction.lower_bound, reaction.upper_bound = scaled_bound

    def regulate_flux(self, reactions):
        '''regulate flux based on True/False activity values for each id in reactions dictionary'''
        for reaction_id, activity in reactions.items():
            reaction = self.model.reactions.get_by_id(reaction_id)

            if not activity:
                # no activity. reaction flux set to 0
                reaction.upper_bound = 0.0
                reaction.lower_bound = 0.0
            elif activity and reaction.bounds == (0.0, 0.0):
                # if new bounds need to be set
                if reaction_id in self.flux_bounds:
                    bounds = self.flux_bounds[reaction_id]
                    scaled_bounds = [b / self.flux_scaling for b in bounds]
                    reaction.lower_bound, reaction.upper_bound = scaled_bounds
                elif reaction_id in self.reversible:
                    reaction.upper_bound = self.default_upper_bound / self.flux_scaling
                    reaction.lower_bound = -self.default_upper_bound / self.flux_scaling
                else:
                    # set bounds based on default
                    reaction.upper_bound = self.default_upper_bound / self.flux_scaling
                    reaction.lower_bound = 0.0

    def objective_value(self):
        if self.solution and self.solution.objective_value:
            objective_value = self.solution.objective_value * self.flux_scaling
            return objective_value
        else:
            return float('nan')

    def optimize(self):

        if self.moma:
            self.solution = cobra.flux_analysis.moma(self.model, solution=self.solution)
        else:
            self.solution = self.model.optimize()

        return self.objective_value()

    def external_reactions(self):
        return [
            EXTERNAL_PREFIX + molecule
            for molecule in self.external_molecules]

    def internal_reactions(self):
        all_reactions = set(self.reaction_ids())
        return all_reactions - set(self.external_reactions())

    def read_fluxes(self, molecules):
        return {
            molecule: self.solution.fluxes[molecule] * self.flux_scaling
            for molecule in molecules}

    def read_internal_fluxes(self):
        return self.read_fluxes(self.internal_reactions())

    def read_exchange_reactions(self):
        '''leaves the prefix in exchange reactions keys'''
        return self.read_fluxes(self.external_reactions())

    def read_exchange_fluxes(self):
        '''removes the prefix from exchange reactions keys, leaving only the external molecule id'''
        external = self.external_reactions()
        levels = self.read_fluxes(external)
        return {
            molecule[len(EXTERNAL_PREFIX):len(molecule)]: level
            for molecule, level in levels.items()}

    def reaction_ids(self):
        return [reaction.id for reaction in self.model.reactions]

    def get_reactions(self, reactions=[]):
        if not reactions:
            reactions = self.reaction_ids()

        return {
            reaction: self.model.reactions.get_by_id(reaction)
            for reaction in reactions}

    def get_reaction_bounds(self, reactions=[]):
        return {
            reaction_key: (reaction.lower_bound, reaction.upper_bound)
            for reaction_key, reaction in self.get_reactions(reactions).items()}

    def __enter__(self):
        return self

    def __exit__(self ,type, value, traceback):
        return False



# def test_canonical():
    # fba = JsonFBA('../bigg_models/e_coli_core.json')
    # fba = JsonFBA('../bigg_models/iWFL_1372.json')
    # fba = FBA({'model_path': '../bigg_models/e_coli_core.json'})
    # fba = FBA({'model_path': '../bigg_models/iWFL_1372.json'})
    # return fba


# if __name__ == '__main__':

    # fba = test_canonical()

    # cobra.io.save_json_model(fba.model, 'demo_model.json')
    # print('MODEL: {}'.format(fba.model))
    # print('REACTIONS: {}'.format(fba.model.reactions))
    # print('METABOLITES: {}'.format(fba.model.metabolites))
    # print('GENES: {}'.format(fba.model.genes))
    # print('COMPARTMENTS: {}'.format(fba.model.compartments))
    # print('SOLVER: {}'.format(fba.model.solver))
    # print('EXPRESSION: {}'.format(fba.model.objective.expression))

    # print(fba.get_reactions())

    # df = pd.DataFrame(fba.get_reactions())

    # print(df)

    # df.to_csv('test.csv', index=False)



    # print(fba.optimize())
    # print(fba.model.summary())


