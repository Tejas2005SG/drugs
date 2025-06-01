from flask import Flask, jsonify
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
mongo_client = MongoClient('mongodb+srv://bhangaletejas003:G0yEjQa9yrTChtDU@h2skill.nnmre.mongodb.net/?retryWrites=true&w=majority&appName=h2skill')
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

# RDKit Reaction SMARTS
reaction_smarts_map = {
    'Esterification': {
        'smarts': '[C:1][OH:2].[C:3](=[O])[OH:4]>>[C:1][O:2][C:3](=[O]).[OH2:4]',
        'description': 'Alcohol + Carboxylic Acid -> Ester + Water'
    },
    'AmideFormation': {
        'smarts': '[N:1].[C:2](=[O])[OH:3]>>[N:1][C:2](=[O]).[OH2:3]',
        'description': 'Amine + Carboxylic Acid -> Amide + Water'
    },
    'Hydrolysis': {
        'smarts': '[C:1][O:2][C:3](=[O])>>[C:1][OH:2].[C:3](=[O])[OH]',
        'description': 'Ester -> Alcohol + Carboxylic Acid'
    },
    'AmideHydrolysis': {
        'smarts': '[C:1][C:2](=[O])[N:3]>>[C:1][C:2](=[O])[OH].[N:3]',
        'description': 'Amide -> Carboxylic Acid + Amine'
    },
    'ThiocarbonylHydrolysis': {
        'smarts': '[C:1](=[S:2])>>[C:1](=[O:2]).[H][S][H]',
        'description': 'Thiocarbonyl -> Carbonyl + H2S'
    },
    'ThioureaFormation': {
        'smarts': '[N:1][H].[C:2]=[S:3]>>[H][N:1][C:2](=[S:3])[N:1][H]',
        'description': 'Amine + Thiocarbonyl -> Thiourea Derivative'
    },
    'SecondaryAmineThiourea': {
        'smarts': '[N:1]([H])[C:4].[C:2]=[S:3]>>[C:4][N:1][C:2](=[S:3])[N:1][H]',
        'description': 'Secondary Amine + Thiocarbonyl -> Substituted Thiourea'
    },
    'Dehalogenation': {
        'smarts': '[c:1][I:2].[H][O:3][H]>>[c:1][H].[I:2][O:3][H]',
        'description': 'Aryl Iodide + Water -> Aryl Hydrogen + Hypoiodous Acid'
    }
}

# Functional group detection patterns
functional_group_patterns = {
    'alcohol': Chem.MolFromSmarts('[#6][OH1]'),
    'amine': Chem.MolFromSmarts('[#7]'),
    'secondary_amine': Chem.MolFromSmarts('[#7H1]'),
    'carbonyl': Chem.MolFromSmarts('[#6]=[O]'),
    'carboxylic_acid': Chem.MolFromSmarts('[#6](=[O])[OH]'),
    'aryl_halide': Chem.MolFromSmarts('[c][Cl,Br,I]'),
    'boronic_acid': Chem.MolFromSmarts('[B](O)(O)'),
    'ester': Chem.MolFromSmarts('[#6][O][C](=[O])'),
    'alkyl_halide': Chem.MolFromSmarts('[#6][Cl,Br,I]'),
    'alkene': Chem.MolFromSmarts('[#6]=[#6]'),
    'alkyne': Chem.MolFromSmarts('[#6]#[#6]'),
    'nitro': Chem.MolFromSmarts('[N](=O)(=O)'),
    'amide': Chem.MolFromSmarts('[#6][C](=[O])[N]'),
    'thiocarbonyl': Chem.MolFromSmarts('[#6]=[S]'),
    'thiourea': Chem.MolFromSmarts('[N][C](=[S])[N]')
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

def detect_functional_groups(mol):
    """Detect functional groups in a molecule."""
    if not mol:
        return []
    groups = []
    for group_name, pattern in functional_group_patterns.items():
        if pattern and mol.HasSubstructMatch(pattern):
            groups.append(group_name)
    logger.info(f"Detected functional groups for molecule: {groups}")
    return groups

def preprocess_molecule(mol, target_group=None):
    """Preprocess a molecule to enable reactions."""
    if not mol:
        return None
    try:
        groups = detect_functional_groups(mol)
        if target_group == 'carbonyl' and 'thiocarbonyl' in groups:
            logger.info("Preprocessing: Converting thiocarbonyl to carbonyl")
            rxn_smarts = rdChemReactions.ReactionFromSmarts(reaction_smarts_map['ThiocarbonylHydrolysis']['smarts'])
            products = rxn_smarts.RunReactants((mol,))
            if products:
                for prod_set in products:
                    for prod in prod_set:
                        if prod.HasSubstructMatch(Chem.MolFromSmarts('[#6]=[O]')):
                            return standardize_molecule(prod)
        if target_group == 'carboxylic_acid' and 'amide' in groups:
            logger.info("Preprocessing: Hydrolyzing amide to carboxylic acid")
            rxn_smarts = rdChemReactions.ReactionFromSmarts(reaction_smarts_map['AmideHydrolysis']['smarts'])
            products = rxn_smarts.RunReactants((mol,))
            if products:
                for prod_set in products:
                    for prod in prod_set:
                        if prod.HasSubstructMatch(Chem.MolFromSmarts('[#6](=[O])[OH]')):
                            return standardize_molecule(prod)
        if target_group == 'aryl_hydrogen' and 'aryl_halide' in groups:
            logger.info("Preprocessing: Dehalogenating aryl iodide")
            rxn_smarts = rdChemReactions.ReactionFromSmarts(reaction_smarts_map['Dehalogenation']['smarts'])
            products = rxn_smarts.RunReactants((mol,))
            if products:
                for prod_set in products:
                    for prod in prod_set:
                        if not prod.HasSubstructMatch(Chem.MolFromSmarts('[c][I]')):
                            return standardize_molecule(prod)
        return mol
    except Exception as e:
        logger.warning(f"Preprocessing failed: {str(e)}")
        return mol

def standardize_molecule(mol):
    """Standardize molecule (tautomers, protonation, stereochemistry)."""
    if not mol:
        return None
    try:
        standardizer = MolStandardize.Standardizer()
        mol = standardizer.standardize(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        return mol
    except Exception as e:
        logger.warning(f"Standardization failed: {str(e)}")
        return mol

def compute_product_properties(mol):
    """Compute chemical properties for a molecule."""
    try:
        Chem.SanitizeMol(mol)
        return {
            'smiles': Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True),
            'molecular_weight': Descriptors.MolWt(mol),
            'num_atoms': mol.GetNumAtoms(),
            'logP': Descriptors.MolLogP(mol),
            'tpsa': Descriptors.TPSA(mol),
            'num_h_donors': Descriptors.NumHDonors(mol),
            'num_h_acceptors': Descriptors.NumHAcceptors(mol),
            'num_rotatable_bonds': Descriptors.NumRotatableBonds(mol),
            'has_stereochemistry': any(atom.HasProp('_ChiralityPossible') for atom in mol.GetAtoms()),
            'functional_groups': detect_functional_groups(mol) or []
        }
    except Exception as e:
        logger.warning(f"Property computation failed: {str(e)}")
        return {'error': f'Property computation failed: {str(e)}', 'functional_groups': []}

def predict_reaction_with_rxn(reactant1_smiles, reactant2_smiles):
    """Predict reaction products using RXN4ChemistryWrapper."""
    if not rxn:
        logger.error("RXN4ChemistryWrapper not initialized. Falling back to RDKit.")
        return None, 0.0

    try:
        reactants = f"{reactant1_smiles}>>{reactant2_smiles}"
        response = rxn.predict_reaction(reactants)
        logger.info(f"RXN4Chemistry predict_reaction response: {response}")
        
        time.sleep(5)
        prediction_id = response['prediction_id']
        results = rxn.get_predict_reaction_results(prediction_id)
        
        max_attempts = 10
        for _ in range(max_attempts):
            if 'results' in results:
                break
            time.sleep(5)
            results = rxn.get_predict_reaction_results(prediction_id)
        
        if 'results' not in results or not results['results']:
            logger.warning("No predictions returned from RXN4Chemistry")
            return None, 0.0
        
        top_prediction = results['results'][0]
        product_smiles = top_prediction['smiles']
        confidence = top_prediction.get('confidence', 0.0)
        return product_smiles, confidence
    except Exception as e:
        logger.error(f"RXN4Chemistry prediction failed: {str(e)}")
        return None, 0.0

def predict_reaction_with_rdkit(mol1, mol2, smiles1, smiles2):
    """Fallback to RDKit reaction prediction if RXN API fails."""
    groups1 = detect_functional_groups(mol1)
    groups2 = detect_functional_groups(mol2)
    applicable_reactions = []

    if 'alcohol' in groups1 and 'carboxylic_acid' in groups2 or 'alcohol' in groups2 and 'carboxylic_acid' in groups1:
        applicable_reactions.append('Esterification')
    if 'amine' in groups1 and 'carboxylic_acid' in groups2 or 'amine' in groups2 and 'carboxylic_acid' in groups1:
        applicable_reactions.append('AmideFormation')
    if 'ester' in groups1 or 'ester' in groups2:
        applicable_reactions.append('Hydrolysis')
    if 'amide' in groups1 or 'amide' in groups2:
        applicable_reactions.append('AmideHydrolysis')
    if 'amine' in groups1 and 'thiocarbonyl' in groups2 or 'amine' in groups2 and 'thiocarbonyl' in groups1:
        applicable_reactions.append('ThioureaFormation')
    if 'secondary_amine' in groups1 and 'thiocarbonyl' in groups2 or 'secondary_amine' in groups2 and 'thiocarbonyl' in groups1:
        applicable_reactions.append('SecondaryAmineThiourea')
    if 'aryl_halide' in groups1 or 'aryl_halide' in groups2:
        applicable_reactions.append('Dehalogenation')

    mol1_pre = mol1
    mol2_pre = mol2
    if 'thiocarbonyl' in groups1:
        mol1_pre = preprocess_molecule(mol1, 'carbonyl')
        groups1 = detect_functional_groups(mol1_pre)
    if 'thiocarbonyl' in groups2:
        mol2_pre = preprocess_molecule(mol2, 'carbonyl')
        groups2 = detect_functional_groups(mol2_pre)
    if 'amide' in groups1:
        mol1_pre = preprocess_molecule(mol1, 'carboxylic_acid')
        groups1 = detect_functional_groups(mol1_pre)
    if 'amide' in groups2:
        mol2_pre = preprocess_molecule(mol2, 'carboxylic_acid')
        groups2 = detect_functional_groups(mol2_pre)
    if 'aryl_halide' in groups1:
        mol1_pre = preprocess_molecule(mol1, 'aryl_hydrogen')
        groups1 = detect_functional_groups(mol1_pre)
    if 'aryl_halide' in groups2:
        mol2_pre = preprocess_molecule(mol2, 'aryl_hydrogen')
        groups2 = detect_functional_groups(mol2_pre)

    if 'alcohol' in groups1 and 'carboxylic_acid' in groups2 or 'alcohol' in groups2 and 'carboxylic_acid' in groups1:
        applicable_reactions.append('Esterification')
    if 'amine' in groups1 and 'carboxylic_acid' in groups2 or 'amine' in groups2 and 'carboxylic_acid' in groups1:
        applicable_reactions.append('AmideFormation')

    if not applicable_reactions:
        logger.info(f"No applicable RDKit reactions for groups: {groups1}, {groups2}")
        return None, 0.0, None

    for reaction_type in applicable_reactions:
        try:
            logger.info(f"Attempting RDKit reaction: {reaction_type}")
            rxn_smarts = rdChemReactions.ReactionFromSmarts(reaction_smarts_map[reaction_type]['smarts'])
            rxn_smarts.Initialize()
            reactants = (mol1_pre, mol2_pre) if len(rxn_smarts.GetReactants()) == 2 else (mol1_pre,)
            products = rxn_smarts.RunReactants(reactants)
            
            if products:
                product_mols = []
                for product_set in products:
                    product_mols.extend(product_set)
                product_smiles = '.'.join(Chem.MolToSmiles(mol, isomericSmiles=True) for mol in product_mols)
                confidence = 0.8
                logger.info(f"Successful RDKit reaction: {reaction_type}, products: {product_smiles}")
                return product_mols, confidence, reaction_type
            else:
                logger.warning(f"RDKit reaction {reaction_type} produced no products")
        except Exception as e:
            logger.warning(f"RDKit reaction {reaction_type} failed: {str(e)}")
            continue
    
    return None, 0.0, None

def process_product_smiles(product_smiles):
    """Convert product SMILES to molecule and split if needed."""
    if not product_smiles:
        return []
    
    products = []
    mol = Chem.MolFromSmiles(product_smiles, sanitize=False)
    if mol:
        try:
            Chem.SanitizeMol(mol)
            products.append(mol)
        except:
            pass
    
    if not products:
        for part in product_smiles.split('.'):
            part_mol = Chem.MolFromSmiles(part, sanitize=False)
            if part_mol:
                try:
                    Chem.SanitizeMol(part_mol)
                    products.append(part_mol)
                except:
                    continue
    
    return products

@app.route('/api/react', methods=['GET'])
def handle_reaction():
    try:
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

        logger.info(f"Found {len(ligand_smiles)} valid ligand SMILES: {ligand_smiles}")

        all_reaction_results = []
        failed_reactions = []
        molecular_weights = []
        processed_pairs = set()

        for smiles1, smiles2 in combinations(ligand_smiles, 2):
            pair_key = tuple(sorted([smiles1, smiles2]))
            if pair_key in processed_pairs:
                continue
            processed_pairs.add(pair_key)

            mol1 = Chem.MolFromSmiles(smiles1)
            mol2 = Chem.MolFromSmiles(smiles2)
            if not mol1 or not mol2:
                failed_reactions.append({
                    'reactants': [smiles1, smiles2],
                    'reason': 'Invalid SMILES string for one or both reactants'
                })
                logger.warning(f"Invalid SMILES: {smiles1}, {smiles2}")
                continue

            try:
                mol1 = standardize_molecule(mol1)
                mol2 = standardize_molecule(mol2)
                Chem.SanitizeMol(mol1)
                Chem.SanitizeMol(mol2)
                
                groups1 = detect_functional_groups(mol1)
                groups2 = detect_functional_groups(mol2)
                
                product_smiles, confidence = predict_reaction_with_rxn(smiles1, smiles2)
                time.sleep(0.2)
                reaction_type = 'RXN4Chemistry_Prediction'

                if not product_smiles or confidence < 0.01:
                    logger.info(f"RXN4Chemistry failed or low confidence ({confidence:.2f}). Falling back to RDKit.")
                    product_mols, confidence, rdkit_reaction_type = predict_reaction_with_rdkit(mol1, mol2, smiles1, smiles2)
                    if product_mols:
                        product_smiles = '.'.join(Chem.MolToSmiles(mol, isomericSmiles=True) for mol in product_mols)
                        reaction_type = rdkit_reaction_type
                    else:
                        failed_reactions.append({
                            'reactants': [smiles1, smiles2],
                            'reason': f'No reaction predicted (RXN confidence: {confidence:.2f}, RDKit failed)'
                        })
                        logger.info(f"No reaction for {smiles1} + {smiles2}")
                        continue
                else:
                    product_mols = process_product_smiles(product_smiles)

                if not product_mols:
                    failed_reactions.append({
                        'reactants': [smiles1, smiles2],
                        'reason': f'Could not parse product SMILES: {product_smiles}'
                    })
                    logger.warning(f"Failed to parse products: {product_smiles}")
                    continue
                
                product_info = []
                for mol in product_mols:
                    try:
                        mol = standardize_molecule(mol)
                        Chem.SanitizeMol(mol)
                        props = compute_product_properties(mol)
                        if 'error' not in props:
                            product_info.append(props)
                            molecular_weights.append(props['molecular_weight'])
                    except Exception as e:
                        logger.warning(f"Product processing failed: {str(e)}")
                        continue
                
                if product_info:
                    all_reaction_results.append({
                        'reactionType': reaction_type,
                        'description': reaction_smarts_map.get(reaction_type, {}).get('description', 'Predicted by RXN4Chemistry'),
                        'reactants': [smiles1, smiles2],
                        'reactantGroups': [groups1, groups2],
                        'products': product_info,
                        'confidence': confidence,
                        'productSmiles': product_smiles
                    })
                else:
                    failed_reactions.append({
                        'reactants': [smiles1, smiles2],
                        'reason': 'All products failed processing'
                    })

            except Exception as e:
                logger.error(f"Error processing pair {smiles1}, {smiles2}: {str(e)}")
                failed_reactions.append({
                    'reactants': [smiles1, smiles2],
                    'reason': f'Processing error: {str(e)}'
                })

        reactant_props = []
        for smiles in ligand_smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                try:
                    mol = standardize_molecule(mol)
                    Chem.SanitizeMol(mol)
                    props = compute_product_properties(mol)
                    if 'error' not in props:
                        reactant_props.append({'smiles': smiles, 'properties': props})
                except Exception as e:
                    logger.warning(f"Failed to compute properties for {smiles}: {str(e)}")
            else:
                logger.warning(f"Invalid SMILES: {smiles}")

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
            'createdAt': datetime.utcnow()
        }

        # Validate the response against the Pydantic schema
        try:
            validated_response = ReactionResponse(**response)
        except ValidationError as e:
            logger.error(f"Response validation failed: {str(e)}")
            return jsonify({'error': f'Response validation failed: {str(e)}'}), 500

        # Save the validated response to MongoDB
        try:
            reaction_responses_collection.insert_one(validated_response.dict())
            logger.info("Successfully saved validated reaction response to MongoDB")
        except Exception as e:
            logger.error(f"Failed to save response to MongoDB: {str(e)}")
            return jsonify({'error': f'Failed to save response to MongoDB: {str(e)}'}), 500

        logger.info(f"Response generated with {len(all_reaction_results)} successful and {len(failed_reactions)} failed reactions")
        return jsonify(response)

    except Exception as e:
        logger.error(f"Server error: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': f'Server error: {str(e)}'}), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5001, debug=True)