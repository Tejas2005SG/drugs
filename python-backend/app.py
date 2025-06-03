from flask import Flask, jsonify, request
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, MolStandardize, rdChemReactions
from rdkit import RDLogger
import numpy as np
import logging
import traceback
from pymongo import MongoClient
from itertools import combinations
from rxn4chemistry import RXN4ChemistryWrapper
import os
import time
from datetime import datetime
from pydantic import BaseModel, Field, ValidationError
from typing import List, Dict, Any, Union

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.warning')

app = Flask(__name__)
CORS(app, resources={r"/api/*": {"origins": ["http://localhost:5173"]}})

# MongoDB connection
mongo_uri = os.getenv('MONGO_URI', 'mongodb+srv://bhangaletejas003:G0yEjQa9yrTChtDU@h2skill.nnmre.mongodb.net/?retryWrites=true&w=majority&appName=h2skill')
mongo_client = MongoClient(mongo_uri)
db = mongo_client['test']
target_protein_collection = db['targetproteins']
reaction_responses_collection = db['reaction_responses']

# RXN4Chemistry configuration
RXN_API_KEY = os.getenv('RXN_API_KEY', 'apk-c01d9c918e326ad080245e040134343481fe4a7513d859ae56743d11740cf93e')
if not RXN_API_KEY or RXN_API_KEY == 'apk-c01d9c918e326ad080245e040134343481fe4a7513d859ae56743d11740cf93e':
    logger.error("RXN_API_KEY not set. Reaction prediction will rely on RDKit fallback.")
    rxn = None
else:
    try:
        rxn = RXN4ChemistryWrapper(api_key=RXN_API_KEY)
        project_name = "rxn4chemistry_tour"
        projects = rxn.get_projects()
        project_exists = any(project['name'] == project_name for project in projects)
        if not project_exists:
            rxn.create_project(project_name)
        rxn.set_project(project_name)
        logger.info(f"Using RXN4Chemistry project ID: {rxn.project_id}")
    except Exception as e:
        logger.error(f"Failed to initialize RXN4Chemistry: {str(e)}")
        rxn = None

# Comprehensive RDKit Reaction SMARTS - Greatly Expanded
reaction_smarts_map = {
    # Basic Organic Reactions
    'Esterification': {
        'smarts': '[C:1][OH:2].[C:3](=[O])[OH:4]>>[C:1][O:2][C:3](=[O]).[OH2:4]',
        'description': 'Alcohol + Carboxylic Acid -> Ester + Water',
        'priority': 8
    },
    'AmideFormation': {
        'smarts': '[N:1].[C:2](=[O])[OH:3]>>[N:1][C:2](=[O]).[OH2:3]',
        'description': 'Amine + Carboxylic Acid -> Amide + Water',
        'priority': 8
    },
    'Hydrolysis': {
        'smarts': '[C:1][O:2][C:3](=[O])>>[C:1][OH:2].[C:3](=[O])[OH]',
        'description': 'Ester -> Alcohol + Carboxylic Acid',
        'priority': 6
    },
    'AmideHydrolysis': {
        'smarts': '[C:1][C:2](=[O])[N:3]>>[C:1][C:2](=[O])[OH].[N:3]',
        'description': 'Amide -> Carboxylic Acid + Amine',
        'priority': 6
    },
    
    # Nucleophilic Substitution Reactions
    'SN2_Alcohol': {
        'smarts': '[C:1][Cl,Br,I:2].[OH2:3]>>[C:1][OH:3].[Cl,Br,I:2]',
        'description': 'Alkyl Halide + Water -> Alcohol + Halide',
        'priority': 7
    },
    'SN2_Amine': {
        'smarts': '[C:1][Cl,Br,I:2].[N:3]>>[C:1][N:3].[Cl,Br,I:2]',
        'description': 'Alkyl Halide + Amine -> Substituted Amine + Halide',
        'priority': 7
    },
    'SN2_Thiol': {
        'smarts': '[C:1][Cl,Br,I:2].[S:3][H]>>[C:1][S:3].[Cl,Br,I:2][H]',
        'description': 'Alkyl Halide + Thiol -> Thioether + Hydrogen Halide',
        'priority': 7
    },
    'SN2_Alkoxide': {
        'smarts': '[C:1][Cl,Br,I:2].[C:3][O-:4]>>[C:1][O:4][C:3].[Cl,Br,I:2]',
        'description': 'Alkyl Halide + Alkoxide -> Ether + Halide',
        'priority': 7
    },
    
    # Addition Reactions
    'HalogenAddition': {
        'smarts': '[C:1]=[C:2].[Cl:3][Cl:4]>>[C:1]([Cl:3])[C:2][Cl:4]',
        'description': 'Alkene + Halogen -> Dihalide',
        'priority': 8
    },
    'HydrogenHalideAddition': {
        'smarts': '[C:1]=[C:2].[H:3][Cl,Br,I:4]>>[C:1]([H:3])[C:2][Cl,Br,I:4]',
        'description': 'Alkene + Hydrogen Halide -> Alkyl Halide',
        'priority': 8
    },
    'WaterAddition': {
        'smarts': '[C:1]=[C:2].[OH2:3]>>[C:1]([OH:3])[C:2][H]',
        'description': 'Alkene + Water -> Alcohol (Hydration)',
        'priority': 7
    },
    'HydrogenAddition': {
        'smarts': '[C:1]=[C:2].[H:3][H:4]>>[C:1]([H:3])[C:2][H:4]',
        'description': 'Alkene + Hydrogen -> Alkane (Hydrogenation)',
        'priority': 8
    },
    
    # Elimination Reactions
    'Dehydration': {
        'smarts': '[C:1][C:2]([OH:3])[H:4]>>[C:1]=[C:2].[OH2:3]',
        'description': 'Alcohol -> Alkene + Water (Dehydration)',
        'priority': 6
    },
    'Dehydrohalogenation': {
        'smarts': '[C:1][C:2]([Cl,Br,I:3])[H:4]>>[C:1]=[C:2].[H:4][Cl,Br,I:3]',
        'description': 'Alkyl Halide -> Alkene + Hydrogen Halide',
        'priority': 6
    },
    
    # Oxidation Reactions
    'AlcoholToAldehyde': {
        'smarts': '[C:1][CH2:2][OH:3]>>[C:1][C:2](=[O:3])[H]',
        'description': 'Primary Alcohol -> Aldehyde',
        'priority': 7
    },
    'AlcoholToKetone': {
        'smarts': '[C:1][CH:2]([OH:3])[C:4]>>[C:1][C:2](=[O:3])[C:4]',
        'description': 'Secondary Alcohol -> Ketone',
        'priority': 7
    },
    'AldehydeToAcid': {
        'smarts': '[C:1][C:2](=[O:3])[H:4]>>[C:1][C:2](=[O:3])[OH:4]',
        'description': 'Aldehyde -> Carboxylic Acid',
        'priority': 7
    },
    
    # Reduction Reactions
    'CarbonylToAlcohol': {
        'smarts': '[C:1][C:2](=[O:3])[H:4].[H:5][H:6]>>[C:1][C:2]([OH:3])[H:4]',
        'description': 'Aldehyde -> Primary Alcohol',
        'priority': 7
    },
    'KetoneToAlcohol': {
        'smarts': '[C:1][C:2](=[O:3])[C:4].[H:5][H:6]>>[C:1][C:2]([OH:3])[C:4]',
        'description': 'Ketone -> Secondary Alcohol',
        'priority': 7
    },
    'CarboxylicAcidToAlcohol': {
        'smarts': '[C:1][C:2](=[O:3])[OH:4].[H:5][H:6]>>[C:1][C:2][OH:4]',
        'description': 'Carboxylic Acid -> Primary Alcohol',
        'priority': 6
    },
    
    # Aromatic Substitution
    'Nitration': {
        'smarts': '[c:1][H:2].[N:3](=[O:4])(=[O:5])[OH:6]>>[c:1][N:3](=[O:4])(=[O:5]).[OH2:6]',
        'description': 'Aromatic -> Nitro Compound',
        'priority': 7
    },
    'Halogenation_Aromatic': {
        'smarts': '[c:1][H:2].[Cl:3][Cl:4]>>[c:1][Cl:3].[H:2][Cl:4]',
        'description': 'Aromatic Halogenation',
        'priority': 7
    },
    'Friedel_Crafts_Alkylation': {
        'smarts': '[c:1][H:2].[C:3][Cl:4]>>[c:1][C:3].[H:2][Cl:4]',
        'description': 'Friedel-Crafts Alkylation',
        'priority': 6
    },
    'Friedel_Crafts_Acylation': {
        'smarts': '[c:1][H:2].[C:3](=[O:4])[Cl:5]>>[c:1][C:3](=[O:4]).[H:2][Cl:5]',
        'description': 'Friedel-Crafts Acylation',
        'priority': 6
    },
    
    # Cross-Coupling Reactions
    'Suzuki_Coupling': {
        'smarts': '[c:1][Br,I:2].[c:3][B:4]([OH:5])[OH:6]>>[c:1][c:3].[Br,I:2][B:4]([OH:5])[OH:6]',
        'description': 'Suzuki Cross-Coupling',
        'priority': 9
    },
    'Heck_Reaction': {
        'smarts': '[c:1][Br,I:2].[C:3]=[C:4]>>[c:1][C:3]=[C:4].[Br,I:2]',
        'description': 'Heck Reaction',
        'priority': 8
    },
    'Sonogashira_Coupling': {
        'smarts': '[c:1][Br,I:2].[C:3]#[C:4][H:5]>>[c:1][C:3]#[C:4].[Br,I:2][H:5]',
        'description': 'Sonogashira Coupling',
        'priority': 8
    },
    
    # Cycloaddition Reactions
    'Diels_Alder': {
        'smarts': '[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2][C:5][C:6][C:4][C:3]1',
        'description': 'Diels-Alder Cycloaddition',
        'priority': 9
    },
    'Cyclopropanation': {
        'smarts': '[C:1]=[C:2].[C:3][H:4][Cl:5]>>[C:1]1[C:2][C:3]1.[H:4][Cl:5]',
        'description': 'Cyclopropanation',
        'priority': 7
    },
    
    # Rearrangement Reactions
    'Claisen_Rearrangement': {
        'smarts': '[C:1]=[C:2][C:3][O:4][C:5]=[C:6]>>[C:1][C:2]=[C:3][C:4](=[O])[C:5]=[C:6]',
        'description': 'Claisen Rearrangement',
        'priority': 6
    },
    'Beckmann_Rearrangement': {
        'smarts': '[C:1][C:2](=[N:3][OH:4])[C:5]>>[C:1][N:3][C:2](=[O])[C:5]',
        'description': 'Beckmann Rearrangement',
        'priority': 6
    },
    
    # Specialized Reactions
    'ThiocarbonylHydrolysis': {
        'smarts': '[C:1](=[S:2])>>[C:1](=[O:2]).[H][S][H]',
        'description': 'Thiocarbonyl -> Carbonyl + H2S',
        'priority': 6
    },
    'ThioureaFormation': {
        'smarts': '[N:1][H].[C:2]=[S:3]>>[H][N:1][C:2](=[S:3])[N:1][H]',
        'description': 'Amine + Thiocarbonyl -> Thiourea Derivative',
        'priority': 7
    },
    'SecondaryAmineThiourea': {
        'smarts': '[N:1]([H])[C:4].[C:2]=[S:3]>>[C:4][N:1][C:2](=[S:3])[N:1][H]',
        'description': 'Secondary Amine + Thiocarbonyl -> Substituted Thiourea',
        'priority': 7
    },
    'Dehalogenation': {
        'smarts': '[c:1][I:2].[H][O:3][H]>>[c:1][H].[I:2][O:3][H]',
        'description': 'Aryl Iodide + Water -> Aryl Hydrogen + Hypoiodous Acid',
        'priority': 5
    },
    
    # Condensation Reactions
    'Aldol_Condensation': {
        'smarts': '[C:1][C:2](=[O:3])[H:4].[C:5][C:6](=[O:7])[H:8]>>[C:1][C:2](=[O:3])[C:4]([OH:8])[C:5][C:6](=[O:7])',
        'description': 'Aldol Condensation',
        'priority': 7
    },
    'Mannich_Reaction': {
        'smarts': '[C:1][C:2](=[O:3])[H:4].[N:5].[C:6]=[O:7]>>[C:1][C:2](=[O:3])[C:4]([N:5])[C:6][OH:7]',
        'description': 'Mannich Reaction',
        'priority': 7
    },
    'Knoevenagel_Condensation': {
        'smarts': '[C:1]=[O:2].[C:3][C:4](=[O:5])[OH:6]>>[C:1]=[C:3][C:4](=[O:5])[OH:6].[OH2:2]',
        'description': 'Knoevenagel Condensation',
        'priority': 7
    },
    
    # Grignard Reactions
    'Grignard_Carbonyl': {
        'smarts': '[C:1][Mg:2][Br:3].[C:4]=[O:5]>>[C:1][C:4]([OH:5])[Mg:2][Br:3]',
        'description': 'Grignard + Carbonyl -> Alcohol',
        'priority': 8
    },
    'Grignard_CO2': {
        'smarts': '[C:1][Mg:2][Br:3].[C:4]=[O:5]=[O:6]>>[C:1][C:4](=[O:5])[OH:6].[Mg:2][Br:3]',
        'description': 'Grignard + CO2 -> Carboxylic Acid',
        'priority': 8
    },
    
    # Protecting Group Chemistry
    'Silyl_Protection': {
        'smarts': '[C:1][OH:2].[Si:3]([C:4])([C:5])[Cl:6]>>[C:1][O:2][Si:3]([C:4])([C:5]).[H:6][Cl]',
        'description': 'Alcohol Silyl Protection',
        'priority': 6
    },
    'Acetal_Formation': {
        'smarts': '[C:1]=[O:2].[C:3][OH:4].[C:5][OH:6]>>[C:1]([O:4][C:3])([O:6][C:5])[H:2]',
        'description': 'Aldehyde/Ketone Acetal Protection',
        'priority': 6
    },
    
    # Heterocycle Formation
    'Pyrrole_Formation': {
        'smarts': '[C:1](=[O:2])[C:3](=[O:4]).[N:5][H:6][H:7]>>[c:1]1[c:3][nH:5][c:1][c:3]1',
        'description': 'Diketone + Ammonia -> Pyrrole',
        'priority': 7
    },
    'Imidazole_Formation': {
        'smarts': '[C:1](=[O:2])[C:3](=[O:4]).[N:5][H:6][H:7].[C:8]=[O:9]>>[c:1]1[n:5][c:8][n:5][c:3]1',
        'description': 'Diketone + NH3 + Aldehyde -> Imidazole',
        'priority': 7
    },
    
    # Multi-component Reactions
    'Ugi_Reaction': {
        'smarts': '[C:1]=[O:2].[N:3].[C:4](=[O:5])[OH:6].[C:7][N:8]#[C:9]>>[C:1]([N:3][C:4](=[O:5])[N:8]([C:7])[C:9])[OH:2]',
        'description': 'Ugi Four-Component Reaction',
        'priority': 8
    },
    'Biginelli_Reaction': {
        'smarts': '[C:1]=[O:2].[C:3][C:4](=[O:5])[C:6].[N:7][H:8][H:9]>>[C:1]1[N:7][C:4](=[O:5])[C:6][C:3][N:7]1',
        'description': 'Biginelli Multicomponent Reaction',
        'priority': 7
    }
}

# Enhanced functional group detection patterns
functional_group_patterns = {
    'alcohol': Chem.MolFromSmarts('[#6][OH1]'),
    'phenol': Chem.MolFromSmarts('[c][OH1]'),
    'primary_amine': Chem.MolFromSmarts('[#6][NH2]'),
    'secondary_amine': Chem.MolFromSmarts('[#6][NH1][#6]'),
    'tertiary_amine': Chem.MolFromSmarts('[#6][N]([#6])[#6]'),
    'amine': Chem.MolFromSmarts('[#7]'),
    'primary_amide': Chem.MolFromSmarts('[#6][C](=[O])[NH2]'),
    'secondary_amide': Chem.MolFromSmarts('[#6][C](=[O])[NH1][#6]'),
    'tertiary_amide': Chem.MolFromSmarts('[#6][C](=[O])[N]([#6])[#6]'),
    'amide': Chem.MolFromSmarts('[#6][C](=[O])[N]'),
    'aldehyde': Chem.MolFromSmarts('[#6][C](=[O])[H]'),
    'ketone': Chem.MolFromSmarts('[#6][C](=[O])[#6]'),
    'carbonyl': Chem.MolFromSmarts('[#6]=[O]'),
    'carboxylic_acid': Chem.MolFromSmarts('[#6](=[O])[OH]'),
    'ester': Chem.MolFromSmarts('[#6][O][C](=[O])'),
    'ether': Chem.MolFromSmarts('[#6][O][#6]'),
    'thiol': Chem.MolFromSmarts('[#6][SH1]'),
    'thioether': Chem.MolFromSmarts('[#6][S][#6]'),
    'disulfide': Chem.MolFromSmarts('[#6][S][S][#6]'),
    'aryl_halide': Chem.MolFromSmarts('[c][Cl,Br,I,F]'),
    'alkyl_halide': Chem.MolFromSmarts('[#6][Cl,Br,I,F]'),
    'boronic_acid': Chem.MolFromSmarts('[B](O)(O)'),
    'boronate_ester': Chem.MolFromSmarts('[B]([O][#6])([O][#6])'),
    'alkene': Chem.MolFromSmarts('[#6]=[#6]'),
    'alkyne': Chem.MolFromSmarts('[#6]#[#6]'),
    'aromatic': Chem.MolFromSmarts('c'),
    'benzene_ring': Chem.MolFromSmarts('c1ccccc1'),
    'nitro': Chem.MolFromSmarts('[N+](=O)[O-]'),
    'nitrile': Chem.MolFromSmarts('[#6]#[N]'),
    'isocyanate': Chem.MolFromSmarts('[N]=[C]=[O]'),
    'thiocarbonyl': Chem.MolFromSmarts('[#6]=[S]'),
    'thiourea': Chem.MolFromSmarts('[N][C](=[S])[N]'),
    'urea': Chem.MolFromSmarts('[N][C](=[O])[N]'),
    'sulfonamide': Chem.MolFromSmarts('[S](=O)(=O)[N]'),
    'phosphate': Chem.MolFromSmarts('[P](=O)([OH])([OH])[OH]'),
    'imidazole': Chem.MolFromSmarts('c1[nH]cnc1'),
    'pyridine': Chem.MolFromSmarts('c1ccncc1'),
    'furan': Chem.MolFromSmarts('c1ccoc1'),
    'pyrrole': Chem.MolFromSmarts('c1cc[nH]c1'),
    'thiophene': Chem.MolFromSmarts('c1ccsc1'),
    'indole': Chem.MolFromSmarts('c1ccc2c(c1)cc[nH]2'),
    'quinoline': Chem.MolFromSmarts('c1ccc2c(c1)cccn2'),
    'isoquinoline': Chem.MolFromSmarts('c1ccc2c(c1)cncc2')
}

# Pydantic models for schema validation
class Properties(BaseModel):
    smiles: str
    molecular_weight: float
    num_atoms: int
    logP: float
    tpsa: float
    num_h_donors: int
    num_h_acceptors: int
    num_rotatable_bonds: int
    has_stereochemistry: bool
    functional_groups: List[str]

class Reactant(BaseModel):
    smiles: str
    properties: Properties

class Product(BaseModel):
    smiles: str
    molecular_weight: float
    num_atoms: int
    logP: float
    tpsa: float
    num_h_donors: int
    num_h_acceptors: int
    num_rotatable_bonds: int
    has_stereochemistry: bool
    functional_groups: List[str]

class ReactionResult(BaseModel):
    reactionType: str
    description: str
    reactants: List[str]
    reactantGroups: List[List[str]]
    products: List[Product]
    confidence: float
    productSmiles: str

class FailedReaction(BaseModel):
    reactants: List[str]
    reason: str
    reactionType: str = "Unknown"

class Statistics(BaseModel):
    mean_mw: float
    std_mw: float
    min_mw: float
    max_mw: float

class ReactionResponse(BaseModel):
    reactants: List[Reactant]
    reactionResults: List[ReactionResult]
    failedReactions: List[FailedReaction]
    statistics: Statistics
    createdAt: datetime
    processedReactionTypes: Dict[str, List[str]]

def detect_functional_groups(mol):
    """Enhanced functional group detection."""
    if not mol:
        return []
    groups = []
    for group_name, pattern in functional_group_patterns.items():
        if pattern and mol.HasSubstructMatch(pattern):
            groups.append(group_name)
    logger.info(f"Detected functional groups for molecule: {groups}")
    return groups

def calculate_reaction_compatibility(groups1, groups2, reaction_type):
    """Calculate compatibility score between reactants for a specific reaction."""
    compatibility_rules = {
        'Esterification': {
            'required': [(['alcohol'], ['carboxylic_acid']), (['carboxylic_acid'], ['alcohol'])],
            'score': 0.9
        },
        'AmideFormation': {
            'required': [(['amine', 'primary_amine', 'secondary_amine'], ['carboxylic_acid']),
                         (['carboxylic_acid'], ['amine', 'primary_amine', 'secondary_amine'])],
            'score': 0.9
        },
        'Suzuki_Coupling': {
            'required': [(['aryl_halide'], ['boronic_acid', 'boronate_ester']),
                         (['boronic_acid', 'boronate_ester'], ['aryl_halide'])],
            'score': 0.95
        },
        'Heck_Reaction': {
            'required': [(['aryl_halide'], ['alkene']), (['alkene'], ['aryl_halide'])],
            'score': 0.9
        },
        'SN2_Alcohol': {
            'required': [(['alkyl_halide'], ['alcohol']), (['alcohol'], ['alkyl_halide'])],
            'score': 0.0
        },
        'Diels_Alder': {
            'required': [(['alkene'], ['alkene'])],
            'score': 0.85
        },
        'Grignard': {
            'required': [(['alkyl_halide'], ['carbonyl', 'aldehyde', 'ketone'])],
            'score': 0.9
        }
    }

    if reaction_type not in compatibility_rules:
        return 0.5
    rule = compatibility_rules[reaction_type]
    for req_pair in rule['required']:
        if (any(g in groups1 for g in req_pair[0]) and any(g in groups2 for g in req_pair[1])) or \
           (any(g in groups2 for g in req_pair[0]) and any(g in groups1 for g in req_pair[1])):
            return rule['score']
    return 0.1

def preprocess_molecule(mol, target_group=None):
    """Enhanced preprocessing with more transformations."""
    if not mol:
        return None
    try:
        groups = detect_functional_groups(mol)
        if target_group == 'carbonyl' and 'thiocarbonyl' in groups:
            logger.info("Preprocessing: Converting thiocarbonyl to carbonyl")
            rxn_smarts = rdChemReactions.ReactionFromSmarts(
                '[C:1](=[S:2])>>[C:1](=[O:2]).[H][S][H]'
            )
            products = rxn_smarts.RunReactants((mol,))
            if products:
                for prod_set in products:
                    for prod in prod_set:
                        if prod.HasSubstructMatch(Chem.MolFromSmarts('[#6]=[O]')):
                            return standardize_molecule(prod)
        if target_group == 'carboxylic acid' and 'amide' in groups:
            logger.info("Preprocessing: Hydrolyzing amide to carboxylic acid")
            rxn_smarts = rdChemReactions.ReactionFromSmarts(
                '[C:1][C:2](=[O])[N:3]>>[C:1][C:2](=[O])[OH].[N:3]'
            )
            products = rxn_smarts.RunReactants((mol,))
            if products:
                for prod_set in products:
                    for prod in prod_set:
                        if prod.HasSubstructMatch(Chem.MolFromSmarts('[#6](=[O])[OH]')):
                            return standardize_molecule(prod)
        if target_group == 'aryl_hydrogen' and 'aryl_halide' in groups:
            logger.info("Preprocessing: Dehalogenating aryl halide")
            rxn_smarts = rdChemReactions.ReactionFromSmarts(
                '[c:1][I:2].[H][O:3][H]>>[c:1][H].[I:2][O:3][H]'
            )
            products = rxn_smarts.RunReactants((mol,))
            if products:
                for prod_set in products:
                    for prod in prod_set:
                        if not prod.HasSubstructMatch(Chem.MolFromSmarts('[c][I,Br,Cl,F]')):
                            return standardize_molecule(prod)
        if target_group == 'alcohol' and 'ester' in groups:
            logger.info("Preprocessing: Hydrolyzing ester to alcohol")
            rxn_smarts = rdChemReactions.ReactionFromSmarts(
                '[C:1][O:2][C:3](=[O])>>[C:1][OH:2].[C:3](=[O])[OH]'
            )
            products = rxn_smarts.RunReactants((mol,))
            if products:
                for prod_set in products:
                    for prod in prod_set:
                        if prod.HasSubstructMatch(Chem.MolFromSmarts('[#6][OH]')):
                            return standardize_molecule(prod)
        return mol
    except Exception as e:
        logger.warning(f"Preprocessing failed: {str(e)}")
        return mol

def standardize_molecule(mol):
    """Enhanced standardization with better error handling."""
    if not mol:
        return None
    try:
        mol = Chem.RemoveHs(mol)
        Chem.SanitizeMol(mol)
        standardizer = MolStandardize.Standardizer()
        mol = standardizer.standardize(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        mol = Chem.AddHs(mol)
        return mol
    except Exception as e:
        logger.warning(f"Standardization failed: {str(e)}")
        try:
            Chem.SanitizeMol(mol)
            return mol
        except:
            return None

def compute_product_properties(mol):
    """Enhanced property computation with more descriptors."""
    try:
        Chem.SanitizeMol(mol)
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rotbonds = Descriptors.NumRotatableBonds(mol)
        has_stereo = any(atom.HasProp('_ChiralityPossible') for atom in mol.GetAtoms()) or \
                     any(bond.GetStereo() != Chem.BondStereo.STEREONONE for bond in mol.GetBonds())
        return {
            'smiles': Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True),
            'molecular_weight': round(mw, 2),
            'num_atoms': mol.GetNumAtoms(),
            'logP': round(logp, 2),
            'tpsa': round(tpsa, 2),
            'num_h_donors': hbd,
            'num_h_acceptors': hba,
            'num_rotatable_bonds': rotbonds,
            'has_stereochemistry': has_stereo,
            'functional_groups': detect_functional_groups(mol) or []
        }
    except Exception as e:
        logger.warning(f"Property computation failed: {str(e)}")
        return {'error': f'Property computation failed: {str(e)}', 'functional_groups': []}

def predict_reaction_with_rxn(reactant1_smiles, reactant2_smiles):
    """Enhanced RXN4Chemistry prediction with better error handling."""
    if not rxn:
        logger.error("RXN4ChemistryWrapper not initialized. Falling back to RDKit.")
        return None, 0.0
    try:
        reactants = f"{reactant1_smiles}.{reactant2_smiles}>>"
        logger.info(f"Predicting reaction for: {reactants}")
        response = rxn.predict_reaction(reactants)
        logger.info(f"RXN4Chemistry predict_reaction response: {response}")
        if 'prediction_id' not in response:
            logger.error("No prediction_id in RXN response")
            return None, 0.0
        time.sleep(3)
        prediction_id = response['prediction_id']
        results = rxn.get_predict_reaction_results(prediction_id)
        max_attempts = 8
        for attempt in range(max_attempts):
            if 'results' in results and results['results']:
                break
            if attempt < max_attempts - 1:
                time.sleep(2)
                results = rxn.get_predict_reaction_results(prediction_id)
        if 'results' not in results or not results['results']:
            logger.warning("No predictions returned from RXN4Chemistry")
            return None, 0.0
        predictions = results['results']
        if isinstance(predictions, list) and len(predictions) > 0:
            top_prediction = predictions[0]
            product_smiles = top_prediction.get('smiles', '')
            confidence = float(top_prediction.get('confidence', 0.0))
            if confidence < 0.1:
                logger.info(f"RXN confidence too low: {confidence}")
                return None, 0.0
            return product_smiles, confidence
        else:
            logger.warning("Invalid prediction format from RXN4Chemistry")
            return None, 0.0
    except Exception as e:
        logger.error(f"RXN4Chemistry prediction failed: {str(e)}")
        return None, 0.0

def predict_reaction_with_rdkit(mol1, mol2, smiles1, smiles2):
    """Enhanced RDKit prediction with prioritized reactions and better logic."""
    groups1 = detect_functional_groups(mol1)
    groups2 = detect_functional_groups(mol2)
    applicable_reactions = []
    for reaction_type, reaction_data in reaction_smarts_map.items():
        compatibility = calculate_reaction_compatibility(groups1, groups2, reaction_type)
        if compatibility > 0.3:
            priority = reaction_data.get('priority', 5)
            applicable_reactions.append((reaction_type, compatibility, priority))
    applicable_reactions.sort(key=lambda x: (x[1], x[2]), reverse=True)
    if not applicable_reactions:
        logger.info(f"No applicable RDKit reactions for groups: {groups1}, {groups2}")
        return None, 0.0, None
    for reaction_type, compatibility, priority in applicable_reactions:
        try:
            logger.info(f"Attempting RDKit reaction: {reaction_type} (compatibility: {compatibility:.2f})")
            rxn_smarts = rdChemReactions.ReactionFromSmarts(reaction_smarts_map[reaction_type]['smarts'])
            rxn_smarts.Initialize()
            num_reactants = rxn_smarts.GetNumReactantTemplates()
            if num_reactants == 1:
                reactants = (mol1,) if len(groups1) > len(groups2) else (mol2,)
            else:
                reactants = (mol1, mol2)
            products = rxn_smarts.RunReactants(reactants)
            if products and len(products) > 0:
                product_mols = []
                for product_set in products:
                    for prod in product_set:
                        if prod is not None:
                            try:
                                prod = standardize_molecule(prod)
                                if prod:
                                    product_mols.append(prod)
                            except:
                                continue
                if product_mols:
                    product_smiles = '.'.join(Chem.MolToSmiles(mol, isomericSmiles=True)
                                            for mol in product_mols if mol)
                    confidence = min(0.95, compatibility + (priority / 10.0))
                    logger.info(f"Successful RDKit reaction: {reaction_type}, products: {product_smiles}")
                    return product_mols, confidence, reaction_type
            else:
                logger.debug(f"RDKit reaction {reaction_type} produced no products")
        except Exception as e:
            logger.warning(f"RDKit reaction {reaction_type} failed: {str(e)}")
            continue
    logger.info("All RDKit reactions failed or produced no valid products")
    return None, 0.0, None

def process_product_smiles(product_smiles):
    """Enhanced product processing with better error handling."""
    if not product_smiles or product_smiles.strip() == "":
        return []
    products = []
    try:
        mol = Chem.MolFromSmiles(product_smiles, sanitize=False)
        if mol:
            Chem.SanitizeMol(mol)
            mol = standardize_molecule(mol)
            if mol:
                products.append(mol)
                return products
    except:
        pass
    separators = ['.', '>>>', '>>', '>', '|']
    parts = [product_smiles]
    for sep in separators:
        new_parts = []
        for part in parts:
            new_parts.extend(part.split(sep))
        parts = new_parts
    for part in parts:
        part = part.strip()
        if not part or len(part) < 2:
            continue
        try:
            part_mol = Chem.MolFromSmiles(part, sanitize=False)
            if part_mol:
                Chem.SanitizeMol(part_mol)
                part_mol = standardize_molecule(part_mol)
                if part_mol and part_mol.GetNumAtoms() > 0:
                    products.append(part_mol)
        except Exception as e:
            logger.debug(f"Failed to process product part '{part}': {str(e)}")
            continue
    return products

def has_reaction_been_processed(reactants, reaction_type):
    """Check if a reaction type has been processed for a reactant pair."""
    sorted_reactants = sorted(reactants)
    pair_key = ':'.join(sorted_reactants)
    query = {
        "processedReactionTypes": {"$elemMatch": {pair_key: {"$in": [reaction_type]}}}
    }
    return reaction_responses_collection.find_one(query) is not None

@app.route('/api/react', methods=['GET'])
def handle_reaction():
    try:
        # Get custom SMILES from query parameters
        smiles1 = request.args.get('smiles1')
        smiles2 = request.args.get('smiles2')
        ligand_smiles = []

        if smiles1 and smiles2:
            # Validate custom SMILES
            mol1 = Chem.MolFromSmiles(smiles1, sanitize=False)
            mol2 = Chem.MolFromSmiles(smiles2, sanitize=False)
            if not mol1 or not mol2:
                logger.error("Invalid custom SMILES provided")
                return jsonify({'error': 'Invalid SMILES strings provided'}), 400
            ligand_smiles = [smiles1, smiles2]
            logger.info(f"Using custom SMILES: {ligand_smiles}")
        else:
            # Fallback to database ligands
            latest_entry = target_protein_collection.find_one(sort=[('createdAt', -1)])
            if not latest_entry:
                logger.error("No TargetProtein entries found")
                return jsonify({'error': 'No TargetProtein entries found in database'}), 404
            ligands = latest_entry.get('TargetLigands', [])
            ligand_smiles = [
                ligand['LigandSmile'] for ligand in ligands
                if ligand.get('LigandSmile') and ligand['LigandSmile'] != 'Not applicable'
            ]
            if len(ligand_smiles) < 2:
                logger.error("Insufficient valid ligand SMILES")
                return jsonify({'error': 'At least two valid ligand SMILES are required'}), 400
            logger.info(f"Found {len(ligand_smiles)} valid ligand SMILES from database: {ligand_smiles}")

        all_reaction_results = []
        failed_reactions = []
        molecular_weights = []
        processed_pairs = set()
        processed_reaction_types = {}

        # Process only the provided SMILES pair or all combinations from database
        if smiles1 and smiles2:
            smiles_pairs = [(smiles1, smiles2)]
        else:
            smiles_pairs = list(combinations(ligand_smiles, 2))

        for smiles1, smiles2 in smiles_pairs:
            sorted_pair = sorted([smiles1, smiles2])
            pair_key = ':'.join(sorted_pair)
            if pair_key in processed_pairs:
                continue
            processed_pairs.add(pair_key)
            processed_reaction_types[pair_key] = []
            mol1 = Chem.MolFromSmiles(smiles1)
            mol2 = Chem.MolFromSmiles(smiles2)
            if not mol1 or not mol2:
                failed_reactions.append({
                    'reactants': [smiles1, smiles2],
                    'reason': 'Invalid SMILES string for one or both reactants',
                    'reactionType': 'Unknown'
                })
                logger.warning(f"Invalid SMILES: {smiles1}, {smiles2}")
                continue
            try:
                mol1 = standardize_molecule(mol1)
                mol2 = standardize_molecule(mol2)
                if not mol1 or not mol2:
                    failed_reactions.append({
                        'reactants': [smiles1, smiles2],
                        'reason': 'Molecule standardization failed',
                        'reactionType': 'Unknown'
                    })
                    continue
                Chem.SanitizeMol(mol1)
                Chem.SanitizeMol(mol2)
                groups1 = detect_functional_groups(mol1)
                groups2 = detect_functional_groups(mol2)
                logger.info(f"Processing {smiles1} (groups: {groups1}) + {smiles2} (groups: {groups2})")
                product_smiles, confidence = None, 0
                reaction_type = 'RXN4Chemistry_Prediction'
                if not has_reaction_been_processed([smiles1, smiles2], reaction_type):
                    product_smiles, confidence = predict_reaction_with_rxn(smiles1, smiles2)
                    processed_reaction_types[pair_key].append(reaction_type)
                else:
                    failed_reactions.append({
                        'reactants': [smiles1, smiles2],
                        'reason': 'RXN4Chemistry prediction skipped (already processed)',
                        'reactionType': reaction_type
                    })
                    logger.info(f"Skipped RXN4Chemistry for {smiles1} + {smiles2} (already processed)")
                if not product_smiles or confidence < 0.15:
                    logger.info(f"RXN4Chemistry failed or low confidence ({confidence:.2f}). Trying RDKit.")
                    product_mols, rdkit_confidence, rdkit_reaction_type = None, 0.0, None
                    if not has_reaction_been_processed([smiles1, smiles2], 'RDKit'):
                        product_mols, rdkit_confidence, rdkit_reaction_type = predict_reaction_with_rdkit(mol1, mol2, smiles1, smiles2)
                        if rdkit_reaction_type:
                            processed_reaction_types[pair_key].append(rdkit_reaction_type)
                    else:
                        failed_reactions.append({
                            'reactants': [smiles1, smiles2],
                            'reason': 'RDKit reaction skipped (already processed)',
                            'reactionType': 'RDKit'
                        })
                        logger.info(f"Skipped RDKit for {smiles1} + {smiles2} (already processed)")
                    if product_mols and rdkit_confidence > 0.3:
                        product_smiles = '.'.join(Chem.MolToSmiles(mol, isomericSmiles=True)
                                                for mol in product_mols if mol)
                        confidence = rdkit_confidence
                        reaction_type = rdkit_reaction_type
                        logger.info(f"RDKit prediction successful: {reaction_type} with confidence {confidence:.2f}")
                    else:
                        failed_reactions.append({
                            'reactants': [smiles1, smiles2],
                            'reason': f'No viable reaction predicted (RXN: {confidence:.2f}, RDKit failed)',
                            'reactionType': 'Unknown'
                        })
                        logger.info(f"No reaction for {smiles1} + {smiles2}")
                        continue
                else:
                    logger.info(f"RXN4Chemistry prediction successful with confidence {confidence:.2f}")
                    product_mols = process_product_smiles(product_smiles)
                if not product_mols:
                    failed_reactions.append({
                        'reactants': [smiles1, smiles2],
                        'reason': f'Could not parse product SMILES: {product_smiles}',
                        'reactionType': reaction_type
                    })
                    logger.warning(f"Failed to parse products: {product_smiles}")
                    continue
                product_info = []
                for mol in product_mols:
                    try:
                        if mol and mol.GetNumAtoms() > 0:
                            mol = standardize_molecule(mol)
                            if mol:
                                Chem.SanitizeMol(mol)
                                props = compute_product_properties(mol)
                                if 'error' not in props:
                                    product_info.append(props)
                                    molecular_weights.append(props['molecular_weight'])
                                else:
                                    logger.warning(f"Product property computation failed: {props['error']}")
                    except Exception as e:
                        logger.warning(f"Product processing failed: {str(e)}")
                        continue
                if product_info:
                    reaction_result = {
                        'reactionType': reaction_type,
                        'description': reaction_smarts_map.get(reaction_type, {}).get('description',
                                                                                     'Predicted by RXN4Chemistry'),
                        'reactants': [smiles1, smiles2],
                        'reactantGroups': [groups1, groups2],
                        'products': product_info,
                        'confidence': round(confidence, 3),
                        'productSmiles': product_smiles
                    }
                    all_reaction_results.append(reaction_result)
                    logger.info(f"Successfully processed reaction: {reaction_type}")
                else:
                    failed_reactions.append({
                        'reactants': [smiles1, smiles2],
                        'reason': 'All products failed property computation',
                        'reactionType': reaction_type
                    })
            except Exception as e:
                logger.error(f"Error processing pair {smiles1}, {smiles2}: {str(e)}")
                failed_reactions.append({
                    'reactants': [smiles1, smiles2],
                    'reason': f'Processing error: {str(e)}',
                    'reactionType': 'Unknown'
                })
        reactant_props = []
        for smiles in set(ligand_smiles):
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                try:
                    mol = standardize_molecule(mol)
                    if mol:
                        Chem.SanitizeMol(mol)
                        props = compute_product_properties(mol)
                        if 'error' not in props:
                            reactant_props.append({'smiles': smiles, 'properties': props})
                except Exception as e:
                    logger.warning(f"Failed to compute properties for {smiles}: {str(e)}")
        mw_stats = {
            'mean_mw': float(np.mean(molecular_weights)) if molecular_weights else 0.0,
            'std_mw': float(np.std(molecular_weights)) if molecular_weights else 0.0,
            'min_mw': float(np.min(molecular_weights)) if molecular_weights else 0.0,
            'max_mw': float(np.max(molecular_weights)) if molecular_weights else 0.0
        }
        response = {
            'reactants': reactant_props,
            'reactionResults': all_reaction_results,
            'failedReactions': failed_reactions,
            'statistics': mw_stats,
            'createdAt': datetime.utcnow(),
            'processedReactionTypes': processed_reaction_types
        }
        try:
            validated_response = ReactionResponse(**response)
        except ValidationError as e:
            logger.error(f"Response validation failed: {str(e)}")
            return jsonify({'error': f'Response validation failed: {str(e)}'}), 500
        try:
            reaction_responses_collection.insert_one(validated_response.dict())
            logger.info("Successfully saved reaction response to MongoDB")
        except Exception as e:
            logger.error(f"Failed to save response to MongoDB: {str(e)}")
            pass
        logger.info(f"Response generated with {len(all_reaction_results)} successful and {len(failed_reactions)} failed reactions")
        return jsonify(response)
    except Exception as e:
        logger.error(f"Server error: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': f'Server error: {str(e)}'}), 500

@app.route('/api/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return jsonify({
        'status': 'healthy',
        'timestamp': datetime.utcnow(),
        'reactions_available': len(reaction_smarts_map),
        'functional_groups_available': len(functional_group_patterns),
        'rxn4chemistry_available': rxn is not None
    })

@app.route('/api/reactions', methods=['GET'])
def get_available_reactions():
    """Get list of available reaction types."""
    reactions_info = []
    for reaction_type, data in reaction_smarts_map.items():
        reactions_info.append({
            'type': reaction_type,
            'description': data['description'],
            'priority': data.get('priority', 5),
            'smarts': data['smarts']
        })
    reactions_info.sort(key=lambda x: x['priority'], reverse=True)
    return jsonify({
        'total_reactions': len(reactions_info),
        'reactions': reactions_info
    })

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5001, debug=True)