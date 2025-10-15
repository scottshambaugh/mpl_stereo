# Installation, Testing, & Publishing

## Basic Installation

```
pip install mpl_stereo
```

## Installing from Source

```
git clone https://github.com/scottshambaugh/mpl_stereo.git
cd mpl_stereo
uv sync --group dev
```

## Running Tests and Type Checking

```
uv run coverage run --source=mpl_stereo -m pytest && uv run coverage report -m 
uv run mypy src tests/test*
```

## Releasing a New Version and Publishing to PyPi

1) Update `CHANGELOG.md`
2) Update the version in `pyproject.toml`
3) Update and install the package
    ```
    uv sync --upgrade --group dev
    ```
4) Run tests, type checking, and linting locally
    ```
    uv run coverage run --source=mpl_stereo -m pytest && uv run coverage report -m 
    uv run mypy src tests/test*
    uv run flake8 src tests docs
    ```
5) Run plotting tests manually
6) Commit any changes and push up the main branch
7) Wait for [CI tests](https://github.com/scottshambaugh/mpl_stereo/actions) to pass
8) [Create a new release](https://github.com/scottshambaugh/mpl_stereo/releases), creating a new tag and including a changelog:    
    ```
    **Changelog**: https://github.com/scottshambaugh/mpl_stereo/blob/main/CHANGELOG.md    
    **Full Diff**: https://github.com/scottshambaugh/mpl_stereo/compare/v0.x.x...v0.x.x
    ```
    This will automatically publish the release to [PyPI](https://pypi.org/project/mpl-stereo/).
