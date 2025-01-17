import typing as t

import typer

app = typer.Typer()


class Molecule:
    name: str
    smiles: str
    fingerprint: t.Any

    def calculate_similarity(self, other_molecule: "Molecule") -> None: ...


class DrugLibrary:
    molecules_list: list[Molecule]

    def find_top_similar(self, top_n: int) -> None: ...


@app.command()
def drug_similarity_evaluator() -> None: ...


if __name__ == "__main__":
    typer.run(app)
