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
from rdchiral.main import rdchiralRun, rdchiralReaction, rdchiralReactants

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

# Drug-relevant reaction SMARTS
reaction_smarts_map = {
    'NucleophilicAddition': {
        'smarts': '[C:1]=[O:2].[N:3]>>[C:1]([O:2])[N:3]',
        'description': 'Carbonyl + Amine -> Secondary Alcohol/Amine Adduct',
        'drug_relevance': 'Common in biguanide and amine drug synthesis'
    },
    'AromaticSubstitution': {
        'smarts': '[c:1][Cl:2].[c:3][NH2:4]>>[c:1][c:3][NH2:4].[Cl:2]',
        'description': 'Aryl Chloride + Amine -> Aryl Amine',
        'drug_relevance': 'Used in kinase inhibitors and CNS drugs'
    },
    'Esterification': {
        'smarts': '[C:1][OH:2].[C:3](=O)[OH:4]>>[C:1][O:2][C:3](=O).[O:4]',
        'description': 'Alcohol + Carboxylic Acid -> Ester + Water',
        'drug_relevance': 'Esters improve lipophilicity'
    },
    'AmideFormation': {
        'smarts': '[C:1][NH2:2].[C:3](=O)[OH:4]>>[C:1][NH:2][C:3](=O).[O:4]',
        'description': 'Amine + Carboxylic Acid -> Amide + Water',
        'drug_relevance': 'Amides are prevalent in biologics'
    },
    'Reduction': {
        'smarts': '[C:1](=O)>>[C:1][OH]',
        'description': 'Ketone/Aldehyde -> Alcohol',
        'drug_relevance': 'Alcohols increase polarity'
    }
}

# Functional group SMARTS for validation
functional_group_smarts = {
    'NucleophilicAddition': ['[C:1]=[O:2]', '[N:3]'],
    'AromaticSubstitution': ['[c:1][Cl:2]', '[c:3][NH2:4]'],
    'Esterification': ['[C:1][OH:2]', '[C:3](=O)[OH:4]'],
    'AmideFormation': ['[C:1][NH2:2]', '[C:3](=O)[OH:4]'],
    'Reduction': ['[C:1](=O)']
}

# Functional group detection patterns
functional_group_patterns = {
    'alcohol': Chem.MolFromSmarts('[#6][OH1]'),
    'amine': Chem.MolFromSmarts('[#6][NH2]'),
    'carbonyl': Chem.MolFromSmarts('[#6]=[O]'),
    'carboxylic_acid': Chem.MolFromSmarts('[#6](=[O])[OH]'),
    'aryl_chloride': Chem.MolFromSmarts('[c][Cl]'),
    'biguanide': Chem.MolFromSmarts('[N:1]=C(N)N=C(N)N')
}

def detect_functional_groups(mol):
    """Detect functional groups in a molecule."""
    if not mol:
        return []
    groups = []
    for group_name, pattern in functional_group_patterns.items():
        if mol.HasSubstructMatch(pattern):
            groups.append(group_name)
    return groups

def validate_reactants(mol1, mol2, reaction_type):
    """Validate if reactants match required functional groups."""
    if reaction_type not in functional_group_smarts:
        return False, f"Unknown reaction type: {reaction_type}"

    required_groups = functional_group_smarts[reaction_type]
    mols = [mol1, mol2] if len(required_groups) == 2 else [mol1]

    for i, (mol, smarts) in enumerate(zip(mols, required_groups)):
        pattern = Chem.MolFromSmarts(smarts)
        if not mol or not pattern or not mol.HasSubstructMatch(pattern):
            group_name = smarts.split(':')[1].split(']')[0] if ':' in smarts else smarts
            return False, f"Reactant {i+1} lacks {group_name} for {reaction_type}"
    return True, ""

def compute_product_properties(mol):
    """Compute chemical properties for a molecule."""
    try:
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

def score_reaction(mol1, mol2, reaction_type):
    """Score reaction feasibility based on functional groups and complexity."""
    groups1 = detect_functional_groups(mol1)
    groups2 = detect_functional_groups(mol2)
    required_groups = {
        'NucleophilicAddition': ['carbonyl', 'amine'],
        'AromaticSubstitution': ['aryl_chloride', 'amine'],
        'Esterification': ['alcohol', 'carboxylic_acid'],
        'AmideFormation': ['amine', 'carboxylic_acid'],
        'Reduction': ['carbonyl']
    }

    score = 0.0
    if reaction_type in required_groups:
        needed = required_groups[reaction_type]
        if len(needed) == 2:
            if needed[0] in groups1 and needed[1] in groups2 or needed[0] in groups2 and needed[1] in groups1:
                score += 0.8
        elif needed[0] in groups1 or needed[0] in groups2:
            score += 0.8
    score -= 0.1 * (mol1.GetNumAtoms() + mol2.GetNumAtoms()) / 100
    return max(0.0, min(1.0, score))

def standardize_molecule(mol):
    """Standardize molecule (tautomers, protonation)."""
    if not mol:
        return None
    try:
        standardizer = MolStandardize.Standardizer()
        return standardizer.standardize(mol)
    except Exception as e:
        logger.warning(f"Standardization failed: {str(e)}")
        return mol

@app.route('/api/react', methods=['GET'])
def handle_reaction():
    try:
        # Fetch latest TargetProtein entry
        latest_entry = target_protein_collection.find_one(sort=[('createdAt', -1)])
        if not latest_entry:
            logger.error("No TargetProtein entries found")
            return jsonify({'error': 'No TargetProtein entries found in database'}), 404

        # Extract ligand SMILES
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

        # Process all pairwise combinations
        for smiles1, smiles2 in combinations(ligand_smiles, 2):
            pair_key = tuple(sorted([smiles1, smiles2]))
            if pair_key in processed_pairs:
                continue
            processed_pairs.add(pair_key)

            # Initialize and standardize molecules
            mol1 = Chem.MolFromSmiles(smiles1, sanitize=False)
            mol2 = Chem.MolFromSmiles(smiles2, sanitize=False)
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
                Chem.SanitizeMol(mol1, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)
                Chem.SanitizeMol(mol2, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)
            except Exception as e:
                failed_reactions.append({
                    'reactants': [smiles1, smiles2],
                    'reason': f'Molecule sanitization/standardization failed: {str(e)}'
                })
                logger.warning(f"Sanitization failed for {smiles1}, {smiles2}: {str(e)}")
                continue

            # Detect functional groups
            groups1 = detect_functional_groups(mol1)
            groups2 = detect_functional_groups(mol2)
            logger.info(f"Functional groups for {smiles1}: {groups1}")
            logger.info(f"Functional groups for {smiles2}: {groups2}")

            # Try all reactions
            for reaction_type in reaction_smarts_map:
                # Validate functional groups
                is_valid, reason = validate_reactants(mol1, mol2, reaction_type)
                if not is_valid:
                    failed_reactions.append({
                        'reactants': [smiles1, smiles2],
                        'reactionType': reaction_type,
                        'reason': reason
                    })
                    logger.info(f"Validation failed for {reaction_type}: {reason}")
                    continue

                # Score reaction feasibility
                confidence = score_reaction(mol1, mol2, reaction_type)
                if confidence < 0.3:
                    failed_reactions.append({
                        'reactants': [smiles1, smiles2],
                        'reactionType': reaction_type,
                        'reason': f'Low reaction feasibility (confidence: {confidence:.2f})'
                    })
                    logger.info(f"Low confidence for {reaction_type}: {confidence:.2f}")
                    continue

                # Run reaction with RDChiral
                try:
                    rxn = rdchiralReaction(reaction_smarts_map[reaction_type]['smarts'])
                    reactants = rdchiralReactants(f'{smiles1}.{smiles2}')
                    products = rdchiralRun(rxn, reactants, keep_isotopes=True)
                except Exception as e:
                    failed_reactions.append({
                        'reactants': [smiles1, smiles2],
                        'reactionType': reaction_type,
                        'reason': f'RDChiral reaction failed: {str(e)}'
                    })
                    logger.error(f"RDChiral failed for {reaction_type}: {str(e)}")
                    continue

                if not products:
                    failed_reactions.append({
                        'reactants': [smiles1, smiles2],
                        'reactionType': reaction_type,
                        'reason': 'No products formed, possibly due to steric hindrance or incompatible substituents'
                    })
                    logger.info(f"No products for {reaction_type} with {smiles1}, {smiles2}")
                    continue

                # Process products
                step_product_info = []
                for product_smiles in products:
                    mol = Chem.MolFromSmiles(product_smiles, sanitize=False)
                    if mol:
                        try:
                            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)
                            props = compute_product_properties(mol)
                            if 'error' not in props:
                                props['confidence'] = confidence
                                step_product_info.append([props])
                                molecular_weights.append(props['molecular_weight'])
                                logger.info(f"Product for {reaction_type}: {props['smiles']}")
                        except Exception as e:
                            logger.warning(f"Product processing failed: {str(e)}")
                            continue

                if step_product_info:
                    all_reaction_results.append({
                        'reactionType': reaction_type,
                        'description': reaction_smarts_map[reaction_type]['description'],
                        'drugRelevance': reaction_smarts_map[reaction_type]['drug_relevance'],
                        'reactants': [smiles1, smiles2],
                        'reactantGroups': [groups1, groups2],
                        'productSets': step_product_info,
                        'confidence': confidence
                    })

        if not all_reaction_results and not failed_reactions:
            logger.error("No reactions processed")
            return jsonify({'error': 'No reactions could be processed'}), 400

        # Statistics
        mw_stats = {
            'mean_mw': float(np.mean(molecular_weights)) if molecular_weights else 0.0,
            'std_mw': float(np.std(molecular_weights)) if molecular_weights else 0.0,
            'min_mw': float(np.min(molecular_weights)) if molecular_weights else 0.0,
            'max_mw': float(np.max(molecular_weights)) if molecular_weights else 0.0
        }

        # Reactant properties
        reactant_props = []
        for smiles in ligand_smiles:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol:
                try:
                    mol = standardize_molecule(mol)
                    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)
                    props = compute_product_properties(mol)
                    if 'error' not in props:
                        reactant_props.append({'smiles': smiles, 'properties': props})
                    else:
                        logger.warning(f"Failed to compute properties for {smiles}: {props['error']}")
                except Exception as e:
                    logger.warning(f"Processing failed for {smiles}: {str(e)}")
                    continue
            else:
                logger.warning(f"Invalid SMILES: {smiles}")

        response = {
            'reactants': reactant_props,
            'reactionResults': all_reaction_results,
            'failedReactions': failed_reactions,
            'statistics': mw_stats
        }
        logger.info(f"Response generated with {len(all_reaction_results)} successful and {len(failed_reactions)} failed reactions")
        return jsonify(response)

    except Exception as e:
        logger.error(f"Error: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': f'Error processing reaction: {str(e)}'}), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5001, debug=True)