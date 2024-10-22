import warnings
import tempfile
import os
from pathlib import Path
import glob
import re
import argparse
import parmed
from plip.structure.preparation import PDBComplex
from rdkit import Chem



def plip_mpro_merge_score(protein, ligand):
    """
    Generate scores from seperate protein & ligand files i.e. not in a complex, by first merging with parmed, then calling
    plip_mpro_score.

    Args:
        protein file
        ligand file

    Returns:
        List of scores/interaction strings.
    """
    with tempfile.TemporaryDirectory() as TD:

        # turn the rdkit.Mol into an SDF file
        if isinstance(ligand, Chem.Mol):
            ligand_path = os.path.join(TD, "ligand.sdf")
            with Chem.SDWriter(ligand_path) as SD:
                SD.write(ligand)
                ligand = ligand_path

        lig = parmed.load_file(ligand)
        if isinstance(lig, list):
            warnings.warn("The ligand was an array (SDF?). Using the first frame. ")
            lig = lig[0]
        protein = parmed.load_file(protein)
        system = protein + lig
        complex_path = os.path.join(TD, "complex.pdb")
        system.save(complex_path, renumber=False)
        return plip_mpro_score(str(complex_path))

def plip_mpro_score(complex_path, ref=None):
    """
    Get all the protein interactions from the interactions between the protein and the ligand. 
    If reference is provided, return interactions+tanimoto similarity, else just return interactions.

    Args:
        complex_path:

    Returns: A list of strings
    """
    complex = PDBComplex()
    complex.load_pdb(complex_path) # Load the PDB file into PLIP class
    complex.analyze()
    n_sites = len(complex.interaction_sets)
    interactions_by_site = []
    # assume there is only one ligand for now
    if n_sites != 1:
        print(f'Warning: detected {n_sites} interaction sites for {complex_path}.')
        if n_sites == 0:
            return []
        #raise ValueError("PLIP detected more (or less) than one ligand?!")
    print(complex.interaction_sets)
    print(complex.interaction_sets.values())
    print(list(complex.interaction_sets.values())[0])
    for site_key, site in complex.interaction_sets.items():
        print(f"Processing site: {site_key}")
        # pair key and values
        interactions = site #list(complex.interaction_sets.values())[0]

        # take all the interactions
        hydrophobic_contacts = ["hydrophobic_" + c.restype + "_" + str(c.resnr) for c in interactions.hydrophobic_contacts]
        # extract protein donors
        hdonors = ["hdonor_" + d.restype + "_" + str(d.resnr) + "_" + str(d.d.atomicnum) for d in interactions.hbonds_pdon]
        # extract protein acceptors
        hacceptors = ["hacceptor_" + a.restype + "_" + str(a.resnr) + "_" + str(a.a.atomicnum) for a in interactions.hbonds_ldon]
        pistacking = ["pistacking_" + r.restype + "_" + str(r.resnr) for r in interactions.pistacking]
        saltbridge = ["saltbridge_" + r.restype + "_" + str(r.resnr) for r in  interactions.saltbridge_pneg]
        waterbridge = ["waterbridge_" + r.restype + "_" + str(r.resnr) for r in interactions.water_bridges]
        pication = ["pication_" + r.restype + "_" + str(r.resnr) for r in interactions.pication_paro]
        halogen_bond = ["halogenbond_" + r.restype + "_" + str(r.resnr) for r in interactions.halogen_bonds]
        metal_complex = ["metalcomplex_" + r.restype + "_" + str(r.resnr) for r in interactions.metal_complexes]

        protein_interaction_fingerprints = (hydrophobic_contacts +
                                            hdonors +
                                            hacceptors +
                                            pistacking +
                                            saltbridge +
                                            waterbridge +
                                            pication +
                                            halogen_bond +
                                            metal_complex)

        if ref:
            intersection = len(mpro_crystal_structures_interactions.intersection(protein_interaction_fingerprints))
            count = len(mpro_crystal_structures_interactions) + len(protein_interaction_fingerprints)
            tanimoto_distance = intersection / (count - intersection)
            interactions_by_site.append(protein_interaction_fingerprints, tanimoto_distance)
        else:
            interactions_by_site.append(protein_interaction_fingerprints)
    return interactions_by_site
parser = argparse.ArgumentParser()
parser.add_argument('data_dir', help='directory to process')
args = parser.parse_args()


print(f'Processing ./{args.data_dir}')
fs = glob.glob(f'{args.data_dir}/*')
# match all files ending in _h
pattern = r'.*pdb'


def glob_re(pattern, strings):
    return filter(re.compile(pattern).match, strings)


complexes = list(glob_re(pattern, fs))

interactions_list = []
# get all interactions as list of lists
for complex in complexes:
    interactions_list.append(plip_mpro_score(complex))

# flatten list of lists so set() can be used
interactions = [interaction for sublist in interactions_list for interaction in sublist]
print(interactions)
# all unique interactions
interactions_set = set(interactions)

print(interactions_set)
mpro_crystal_structures_interactions = {'hacceptor_THR_25_8',
                                         'hdonor_THR_25_8',
                                         'hydrophobic_THR_25',
                                         'hacceptor_HIS_41_8',
                                         'hydrophobic_HIS_41',
                                         'pistacking_HIS_41',
                                         'hacceptor_CYS_44_8',
                                         'hydrophobic_PHE_140',
                                         'hacceptor_ASN_142_8',
                                         'hdonor_ASN_142_7',
                                         'hdonor_GLY_143_7',
                                        'hacceptor_SER_144_8',
                                        'hdonor_SER_144_7',
                                        'hdonor_SER_144_8',
                                        'hdonor_CYS_145_7',
                                        'hacceptor_HIS_163_7',
                                        'hydrophobic_MET_165',
                                        'hacceptor_GLU_166_8',
                                        'hdonor_GLU_166_7',
                                        'hdonor_GLU_166_8',
                                        'hydrophobic_GLU_166',
                                        'saltbridge_GLU_166',
                                        'hydrophobic_PRO_168',
                                        'hydrophobic_ASP_187',
                                        'hacceptor_GLN_189_8',
                                        'hdonor_GLN_189_7',
                                        'hydrophobic_GLN_189'}
mpro_xstal = set(mpro_crystal_structures_interactions)
# find differences
unique_to_set1 = mpro_xstal - interactions_set
unique_to_set2 = interactions_set - mpro_xstal

#print(unique_to_set1, 'og_mpro_xstal')
print(interactions_set)
#print(unique_to_set2, 'interactions')

# results
unique_to_set1, unique_to_set2
with open('interactions.dat', 'w') as f:
    f.write(str(interactions_set))
# pymol select all residues present in interaction set
# pymol select resid 212 resid 224 resid 218 resid 59 resid 162 resid 223 resid 159 resid 223 resid 155 resid 164 resid 223 resid 220 resid 218 resid 159 resid 163 resid 155 resid 211
