import typing as t

import typer
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity

app = typer.Typer()


class Molecule:
    name: str
    smiles: str
    fingerprint: t.Any

    def __init__(self, name: str, smiles: str) -> None:
        self.name = name
        self.smiles = smiles
        molecule = Chem.MolFromSmiles(smiles)
        self.fingerprint = (
            AllChem.GetMorganFingerprintAsBitVect(molecule, radius=2)
            if molecule
            else None
        )

    def calculate_similarity(self, other_molecule: "Molecule") -> float:
        if not self.fingerprint or not other_molecule.fingerprint:
            return 0.0
        return TanimotoSimilarity(self.fingerprint, other_molecule.fingerprint)


class DrugLibrary:
    molecules_list: list[Molecule]

    def find_top_similar(self, top_n: int) -> None: ...


@app.command()
def drug_similarity_evaluator() -> None: ...


if __name__ == "__main__":
    typer.run(app)
