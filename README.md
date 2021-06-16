# Programming Practices Demo

[![Actions Status](https://img.shields.io/github/workflow/status/apallath/prog_practices_demo/build_test)](https://github.com/apallath/prog_practices_demo/actions)

This GitHub repository accompanies my July 2021 tutorial on 'Good Programming
Practices' for new members of the Patel Group.

By using an example Python program to calculate the radius of gyration of a small polymer
and a relatively large protein, this repository (and the accompanying slides) cover the
following topics:

## Modular programming (functions, classes, and objects)
- Modularity using functions:
    - `examples/rg_functions/`

- Object-oriented programming:
    - `analysis_demo/gmx_analysis.py`
    - `examples/rg/`

## Testing
`tests/unit/test_rgyr.py` and `tests/unit/test_selection_parsers.py` demonstrate unit testing.
`tests/integration/test_poly.py` and `tests/integration/test_protein.py` demonstrate the use of 'gold master' integration tests.

## Documentation
Scripts to generate documentation from `analysis_demo`'s docstrings are in `docsource/`.
GitHub pages documentation is hosted from the `docs/` folder.

## Continuous integration
See GitHub Actions for this repository.
