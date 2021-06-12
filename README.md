# Programming Practices Demo

This GitHub repository accompanies my July 2021 tutorial on 'Good Programming
Practices' for new members of the Patel Group.

Using an example Python program to calculate the RMSD of a small polymer
and a small peptide, this repository (and the accompanying slides) cover the
following topics:

## Modular programming (functions, classes, and objects)
- No modularity:
    - `examples/rmsd_no_modular/`

- Modularity using functions:
    - `examples/rmsd_functions/`

- Object-oriented programming:
    - `examples/rmsd_OOP/`
    - `rmsd_demo/rmsd.py`

## Code style
`rmsd_demo/rmsd.py` follows Google's Python code style conventions and PEP8.

## Documentation
Scripts to generate documentation from `rmsd_demo`'s docstrings are in `docsource/`.
GitHub pages documentation is hosted from the `docs/` folder.

## Testing
`tests/unit/test_rmsd_poly.py` and `tests/unit/test_rmsd_complete.py` demonstrate the importance of
unit tests and regression tests. `tests/integration/test_poly.py` and `tests/integration/test_peptide.py` demonstrate the use of 'gold standard' integration tests.

## Version control
See accompanying slides.

## Project structure and packaging
See accompanying slides.

## Continuous integration
See GitHub Actions for this repository.

## Prototyping with Jupyter notebook
`examples/rmsd_jupyter/` demonstrates the use of Jupyter notebooks for quick
prototyping.
