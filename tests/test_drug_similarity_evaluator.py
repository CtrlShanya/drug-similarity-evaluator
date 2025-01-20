from drug_similarity_evaluator.__main__ import Molecule
from rdkit import Chem
from rdkit.Chem import AllChem


def test_basic_molecule_storage():
    name = "Aspirin"
    smiles = "CC (= O ) OC1 = CC = CC = C1C (= O ) O"
    molecule = Molecule(name=name, smiles=smiles)

    expected_mol = Chem.MolFromSmiles(smiles)
    expected_fp = AllChem.GetMorganFingerprintAsBitVect(expected_mol, radius=2)
    assert molecule.fingerprint == expected_fp


def test_basic_tanimoto_calculation():
    name = "Aspirin"
    smiles = "CC (= O ) OC1 = CC = CC = C1C (= O ) O"

    mol_1 = Molecule(name=name, smiles=smiles)
    mol_2 = Molecule(name=name, smiles=smiles)

    assert mol_1.calculate_similarity(mol_2) == 1.0
