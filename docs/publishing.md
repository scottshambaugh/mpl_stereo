# Installation, Testing, & Publishing

## Basic Installation

```
pip install mpl_stereo
```

## Installing from Source

```
git clone https://github.com/scottshambaugh/mpl_stereo.git
cd mpl_stereo
pip install poetry
poetry install
```

## Running Tests and Type Checking

```
poetry run coverage run --source=mpl_stereo -m pytest && poetry run coverage report -m 
poetry run mypy src tests/test*
```

## Releasing a New Version and Publishing to PyPi

1) Update `CHANGELOG.md`
2) Update the version in `pyproject.toml`
3) Update and install the package
    ```
    poetry update
    poetry install
    ```
4) Run tests, type checking, and linting locally
    ```
    poetry run coverage run --source=mpl_stereo -m pytest && poetry run coverage report -m 
    poetry run mypy src tests/test*
    poetry run flake8
    ```
5) Run plotting tests manually
6) Commit any changes and push up the main branch
7) Wait for [CI tests](https://github.com/scottshambaugh/mpl_stereo/actions) to pass
8) [Create a new release](https://github.com/scottshambaugh/mpl_stereo/releases), creating a new tag and including a changelog:    
    ```
    **Changelog**: https://github.com/scottshambaugh/mpl_stereo/blob/main/CHANGELOG.md    
    **Full Diff**: https://github.com/scottshambaugh/mpl_stereo/compare/v0.x.x...v0.x.x
    ```
9) Build wheels: `poetry build`
10) Publish to PyPi: `poetry publish`
11) Check that [the package](https://pypi.org/project/mpl_stereo/) has updated
