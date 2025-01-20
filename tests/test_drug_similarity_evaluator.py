import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from drug_similarity_evaluator.__main__ import Molecule

ASPIRIN = "Aspirin"
ASPIRIN_SMILES = "CC (= O ) OC1 = CC = CC = C1C (= O ) O"


def test_basic_molecule_storage() -> None:
    molecule = Molecule(name=ASPIRIN, smiles=ASPIRIN_SMILES)

    expected_mol = Chem.MolFromSmiles(ASPIRIN_SMILES)
    expected_fp = AllChem.GetMorganFingerprintAsBitVect(expected_mol, radius=2)
    assert molecule.fingerprint == expected_fp


@pytest.mark.parametrize(
    "names, smiles, exp_tani",
    [
        ([ASPIRIN, ASPIRIN], [ASPIRIN_SMILES, ASPIRIN_SMILES], 1.0),
        ([ASPIRIN, "test_mol"], [ASPIRIN_SMILES, ""], 0.0),
    ],
)
def test_basic_tanimoto_calculation(
    names: list[str], smiles: list[str], exp_tani: float
) -> None:

    mol_1 = Molecule(name=names[0], smiles=smiles[0])
    mol_2 = Molecule(name=names[1], smiles=smiles[1])

    assert mol_1.calculate_similarity(mol_2) == exp_tani
