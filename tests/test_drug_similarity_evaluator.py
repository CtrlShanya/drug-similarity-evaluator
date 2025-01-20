from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from drug_similarity_evaluator.__main__ import DrugLibrary, Molecule

ASPIRIN = "Aspirin"
ASPIRIN_SMILES = "CC(=O)OC1=CC=CC=C1C(=O)O"
RESOURCE_DIR = Path(__file__).parent / "resources"
EXAMPLE_FILE = RESOURCE_DIR / "example_file_1.txt"


@pytest.fixture()
def drug_library() -> DrugLibrary:
    drug_lib = DrugLibrary()
    drug_lib.from_file(EXAMPLE_FILE.__str__())
    return drug_lib


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


def test_basic_drug_library_formation(drug_library: DrugLibrary) -> None:
    expected_names = ["Aspirin", "Ibuprofen", "Paracetamol", "Caffeine", "Nicotine"]
    drug_library_names = [x.name for x in drug_library.molecules_list]

    assert len(drug_library.molecules_list) == 5
    assert drug_library_names == expected_names


def test_get_top_n_molecules(drug_library: DrugLibrary) -> None:
    query_mol = Molecule(name=ASPIRIN, smiles=ASPIRIN_SMILES)

    top_n_mols = drug_library.find_top_similar(query_molecule=query_mol, top_n=3)
    result_names = [x[0] for x in top_n_mols]
    expected_names = ["Paracetamol", "Ibuprofen", "Nicotine"]

    assert result_names == expected_names
