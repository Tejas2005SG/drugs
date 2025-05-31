from flask import Flask, request, jsonify
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit import RDLogger
import numpy as np
import logging
import traceback

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.warning')

app = Flask(__name__)
CORS(app, resources={r"/api/*": {"origins": ["http://localhost:5173"]}})

# Expanded reaction SMARTS map for drug discovery
reaction_smarts_map = {
    'Condensation': {
        'smarts': '[C:1][OH:2].[C:3][OH:4]>>[C:1][O:2][C:3].[O:4]',
        'description': 'Alcohol + Alcohol -> Ether + Water (forms ethers, common in drug synthesis)',
        'drug_relevance': 'Ethers are stable linkers in drug molecules'
    },
    'Esterification': {
        'smarts': '[C:1][OH:2].[C:3](=O)[OH:4]>>[C:1][O:2][C:3](=O).[O:4]',
        'description': 'Alcohol + Carboxylic Acid -> Ester + Water (common in prodrugs)',
        'drug_relevance': 'Esters improve lipophilicity and bioavailability'
    },
    'AmideFormation': {
        'smarts': '[C:1][NH2:2].[C:3](=O)[OH:4]>>[C:1][NH:2][C:3](=O).[O:4]',
        'description': 'Amine + Carboxylic Acid -> Amide + Water (key for peptide drugs)',
        'drug_relevance': 'Amides are prevalent in biologics and small molecules'
    },
    'Reduction': {
        'smarts': '[C:1](=O)>>[C:1][OH]',
        'description': 'Ketone/Aldehyde -> Alcohol (reduces carbonyls)',
        'drug_relevance': 'Alcohols increase polarity, useful in metabolites'
    },
    'Oxidation': {
        'smarts': '[C:1][OH:2]>>[C:1](=O)',
        'description': 'Alcohol -> Ketone/Aldehyde (oxidizes alcohols)',
        'drug_relevance': 'Carbonyls are reactive intermediates in synthesis'
    },
    'SuzukiCoupling': {
        'smarts': '[c:1][Br:2].[c:3][B(OH)2:4]>>[c:1][c:3]',
        'description': 'Aryl Bromide + Aryl Boronic Acid -> Biaryl (cross-coupling)',
        'drug_relevance': 'Biaryls are common scaffolds in kinase inhibitors'
    },
    'NucleophilicSubstitution': {
        'smarts': '[C:1][Cl:2].[N:3]>>[C:1][N:3].[Cl:2]',
        'description': 'Alkyl Chloride + Amine -> Alkyl Amine (SN2 reaction)',
        'drug_relevance': 'Amines are key in CNS-active drugs'
    }
}

# Functional group SMARTS for validation
functional_group_smarts = {
    'Condensation': ['[C:1][OH:2]', '[C:3][OH:4]'],
    'Esterification': ['[C:1][OH:2]', '[C:3](=O)[OH:4]'],
    'AmideFormation': ['[C:1][NH2:2]', '[C:3](=O)[OH:4]'],
    'Reduction': ['[C:1](=O)'],
    'Oxidation': ['[C:1][OH:2]'],
    'SuzukiCoupling': ['[c:1][Br:2]', '[c:3][B(OH)2:4]'],
    'NucleophilicSubstitution': ['[C:1][Cl:2]', '[N:3]']
}

def validate_reactants(mol1, mol2, reaction_type):
    """Validate if reactants match the required functional groups."""
    if reaction_type not in functional_group_smarts:
        return False, f"Unknown reaction type: {reaction_type}"

    required_groups = functional_group_smarts[reaction_type]
    mols = [mol1, mol2] if len(required_groups) == 2 else [mol1]

    for i, (mol, smarts) in enumerate(zip(mols, required_groups)):
        pattern = Chem.MolFromSmarts(smarts)
        if not mol or not pattern or not mol.HasSubstructMatch(pattern):
            return False, f"Reactant {i+1} lacks required functional group for {reaction_type}"
    return True, ""

def get_applicable_reactions(mol1, mol2):
    """Identify all applicable reactions."""
    applicable_reactions = []
    for reaction_type in reaction_smarts_map:
        is_valid, _ = validate_reactants(mol1, mol2, reaction_type)
        if is_valid:
            applicable_reactions.append(reaction_type)
    return applicable_reactions

def compute_product_properties(mol):
    """Compute chemical properties for a molecule."""
    return {
        'smiles': Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True),
        'molecular_weight': Descriptors.MolWt(mol),
        'num_atoms': mol.GetNumAtoms(),
        'logP': Descriptors.MolLogP(mol),
        'tpsa': Descriptors.TPSA(mol),
        'num_h_donors': Descriptors.NumHDonors(mol),
        'num_h_acceptors': Descriptors.NumHAcceptors(mol)
    }

@app.route('/api/react', methods=['POST'])
def handle_reaction():
    try:
        # Parse input JSON
        data = request.get_json()
        if not data or not all(key in data for key in ['smiles1', 'smiles2']):
            logger.error("Missing required fields: smiles1, smiles2")
            return jsonify({'error': 'Missing required fields: smiles1, smiles2'}), 400

        smiles1, smiles2 = data['smiles1'], data['smiles2']
        logger.info(f"Received request: SMILES1={smiles1}, SMILES2={smiles2}")

        # Initialize molecules
        mol1 = Chem.MolFromSmiles(smiles1, sanitize=False)
        mol2 = Chem.MolFromSmiles(smiles2, sanitize=False)
        if not mol1 or not mol2:
            logger.error(f"Invalid SMILES: mol1={'valid' if mol1 else 'invalid'}, mol2={'valid' if mol2 else 'invalid'}")
            return jsonify({'error': 'Invalid SMILES strings'}), 400

        # Sanitize molecules
        try:
            Chem.SanitizeMol(mol1, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)
            Chem.SanitizeMol(mol2, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)
        except Exception as e:
            logger.error(f"Sanitization failed: {str(e)}")
            return jsonify({'error': f'Molecule sanitization failed: {str(e)}'}), 400

        logger.info(f"Sanitized mol1: {Chem.MolToSmiles(mol1, isomericSmiles=True)}")
        logger.info(f"Sanitized mol2: {Chem.MolToSmiles(mol2, isomericSmiles=True)}")

        # Identify applicable reactions
        applicable_reactions = get_applicable_reactions(mol1, mol2)
        if not applicable_reactions:
            logger.warning("No applicable reactions found")
            return jsonify({'error': 'No applicable reactions found for the input molecules'}), 400

        logger.info(f"Applicable reactions: {applicable_reactions}")

        all_reaction_results = []
        molecular_weights = []

        # Process each reaction
        for reaction_type in applicable_reactions:
            rxn = AllChem.ReactionFromSmarts(reaction_smarts_map[reaction_type]['smarts'])
            if not rxn:
                logger.error(f"Invalid reaction SMARTS for {reaction_type}")
                continue

            try:
                rxn.Initialize()
            except Exception as e:
                logger.error(f"Reaction initialization failed for {reaction_type}: {str(e)}")
                continue

            try:
                products = rxn.RunReactants([mol1, mol2])
            except Exception as e:
                logger.error(f"Reaction execution failed for {reaction_type}: {str(e)}")
                continue

            if not products:
                logger.warning(f"No products formed for {reaction_type}")
                continue

            logger.info(f"Products for {reaction_type}: {len(products)} product sets")

            # Process products
            step_product_info = []
            for product_set in products:
                product_props = []
                for mol in product_set:
                    if mol:
                        try:
                            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)
                            props = compute_product_properties(mol)
                            product_props.append(props)
                            molecular_weights.append(props['molecular_weight'])
                            logger.info(f"Generated product for {reaction_type}: {props['smiles']}")
                        except Exception as e:
                            logger.warning(f"Failed to process product for {reaction_type}: {str(e)}")
                            continue
                if product_props:
                    step_product_info.append(product_props)

            if step_product_info:
                all_reaction_results.append({
                    'reactionType': reaction_type,
                    'description': reaction_smarts_map[reaction_type]['description'],
                    'drugRelevance': reaction_smarts_map[reaction_type]['drug_relevance'],
                    'productSets': step_product_info
                })

        if not all_reaction_results:
            logger.error("No valid products formed")
            return jsonify({'error': 'No valid products formed for any reaction'}), 400

        # Statistical analysis of products
        mw_stats = {
            'mean_mw': float(np.mean(molecular_weights)) if molecular_weights else 0.0,
            'std_mw': float(np.std(molecular_weights)) if molecular_weights else 0.0,
            'min_mw': float(np.min(molecular_weights)) if molecular_weights else 0.0,
            'max_mw': float(np.max(molecular_weights)) if molecular_weights else 0.0
        }

        # Reactant properties
        reactant1_props = compute_product_properties(mol1)
        reactant2_props = compute_product_properties(mol2)

        response = {
            'reactants': [
                {'smiles': smiles1, 'properties': reactant1_props},
                {'smiles': smiles2, 'properties': reactant2_props}
            ],
            'reactionResults': all_reaction_results,
            'statistics': mw_stats
        }
        logger.info(f"Response: {response}")
        return jsonify(response)

    except Exception as e:
        logger.error(f"Reaction error: {str(e)}")
        traceback.print_exc()
        return jsonify({'error': f'Error processing reaction: {str(e)}'}), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5001, debug=True)