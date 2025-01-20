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
    molecules_list: list[Molecule] = []

    @classmethod
    def from_file(cls, file_path: str) -> None:
        try:
            with open(file_path, "r") as f:
                for line in f:
                    name, smiles = line.strip().split(",")
                    mol = Molecule(name.strip(), smiles.strip())
                    cls.add_molecule(molecule=mol)
        except FileNotFoundError:
            print(f"Error: File '{file_path}' not found.")

    @classmethod
    def add_molecule(cls, molecule: Molecule) -> None:
        if not cls._get_molecule_by_name(molecule.name):
            cls.molecules_list.append(molecule)

    @classmethod
    def _get_molecule_by_name(cls, name: str) -> t.Optional[Molecule]:
        for mol in cls.molecules_list:
            if mol.name == name:
                return mol
        return None

    def _get_tanimoto_similarities(
        self, query_molecule: Molecule
    ) -> list[list[str | float]]:
        similarities: list[list[str | float]] = []

        for mol in self.molecules_list:
            if mol.name != query_molecule.name:
                tanimoto = query_molecule.calculate_similarity(mol)
                similarities.append([mol.name, tanimoto])

        return similarities

    def find_top_similar(
        self, query_molecule: Molecule, top_n: int
    ) -> list[list[str | float]]:
        similarities = self._get_tanimoto_similarities(query_molecule=query_molecule)
        similarities.sort(key=lambda x: x[1], reverse=True)
        return similarities[:top_n]


@app.command()
def drug_similarity_evaluator(
    file_path: str = typer.Option(..., envvar="FILE_PATH"),
    query_name: str = typer.Option(..., envvar="QUERY_MOL_NAME"),
    query_smiles: str = typer.Option(..., envvar="QUERY_MOL_SMILES"),
    top_n: int = typer.Option(..., envvar="TOP_N"),
) -> None:
    drug_library = DrugLibrary()
    drug_library.from_file(file_path=file_path)

    query_molecule = Molecule(name=query_name, smiles=query_smiles)

    similarity_list = drug_library.find_top_similar(
        query_molecule=query_molecule, top_n=top_n
    )
    print("Similar Molecules are: ", similarity_list)


if __name__ == "__main__":
    typer.run(app)
