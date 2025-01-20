# Drug Similarity Evaluator
This project focuses on comparing molecule similarities using Tanimotot Similarity, provided by the package [rdkit](https://www.rdkit.org/), given the SMILES of the various molecules.
It stores molecules in a DrugLibrary and can give the top "N" similar molecules to a given query.

## Set-Up
Please ensure that the `poetry` package (https://python-poetry.org/) is installed for your given python3 version. 
(For now, this project only works for python3.10-3.11)

To install all requirements and set up a poetry environment to run code:
`poetry install`

## Running the App
To run the app, a file with various molecule names and subsequent SMILES is required.
An example can be found in the `tests/resources` folder.
The name and SMILES of a query molecule must also be provided.

To list all the important attributes needed to run the app:
`poetry run drug-similarity --help`

Example run:
`poetry run drug-similarity --file-path \tests\resources\example_file_1.txt --query-name Aspirin --query-smiles "CC(=O)OC1=CC=CC=C1C(=O)O" --top-n 2`

NOTE: Please make sure the query SMILES is encapsulated by double inverted commas.

## Development Mode
To run the formatters and linters:
1. `poetry run black check .`
2. `poetry run mypy .`
3. `poetry run ruff check`

To run the tests:
`poetry run pytest .`
