# Programming Practices Demo

This GitHub repository accompanies my July 2021 tutorial on 'Good Programming
Practices' for new members of the Patel Group.

By using an example Python program to calculate the radius of gyration of a small polymer
and a relatively large protein, this repository (and the accompanying slides) cover the
following topics:

## Modular programming (functions, classes, and objects)
- No modularity:
    - `examples/rg_no_modular/`

- Modularity using functions:
    - `examples/rg_functions/`

- Object-oriented programming:
    - `examples/rg/`
    - `analysis_demo/gmx_analysis.py`

## Testing
`tests/unit/test_rgyr.py` demonstrates unit testing.
`tests/integration/test_poly.py` and `tests/integration/test_protein.py` demonstrate the use of 'gold master' integration tests.

## Documentation
Scripts to generate documentation from `analysis_demo`'s docstrings are in `docsource/`.
GitHub pages documentation is hosted from the `docs/` folder.

## Continuous integration
See GitHub Actions for this repository.
