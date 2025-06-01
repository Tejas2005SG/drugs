// rdkit.js
import rdkitPkg from '@rdkit/rdkit';

let rdkit = null;

export async function initRDKit() {
  if (!rdkit) {
    try {
      // Handle different possible RDKit module structures
      if (typeof rdkitPkg.initRDKitModule === 'function') {
        rdkit = await rdkitPkg.initRDKitModule();
      } else if (typeof rdkitPkg === 'function') {
        rdkit = await rdkitPkg();
      } else if (rdkitPkg.RDKitModule) {
        rdkit = await rdkitPkg.RDKitModule();
      } else {
        throw new Error('No valid RDKit initialization method found');
      }
      console.log('RDKit initialized successfully');
      return rdkit;
    } catch (err) {
      console.error('RDKit initialization failed:', err);
      throw err;
    }
  }
  return rdkit;
}

export function getRDKit() {
  if (!rdkit) {
    throw new Error('RDKit not initialized. Call initRDKit() first.');
  }
  return rdkit;
}