# Repository Guidelines

## Project Structure & Module Organization
- Python package lives in `gopher/`; `gopher/gopher.py` is the CLI entry point wired to `gopher` console script.
- Shared utilities and domain logic sit in modules like `annotations.py`, `enrichment.py`, `graph_search.py`, and `normalize.py`.
- Tests are split by scope: unit tests in `tests/unit_tests/`, integration/system checks in `tests/system_tests/`, and fixtures/data in `tests/data/`.
- Documentation sources are in `docs/`; MkDocs configuration is at `mkdocs.yml`.

## Build, Test, and Development Commands
- Install for development (creates editable install with extras):
  `uv sync --dev` or `pip install -e .[dev]`
- Run the test suite:
  `uv run pytest` or `pytest`
- Lint and format with Ruff (configured in `pyproject.toml`):
  `uv run ruff check .` and `uv run ruff format .`
- Pre-commit hooks (lint/format/metadata checks):
  `pre-commit run --all-files`
- Build docs locally:
  `uv run mkdocs serve`

## Coding Style & Naming Conventions
- Ruff enforces PEP8-style rules; line length is 79. Target Python 3.10+.
- Prefer type hints and docstrings for public functions; docstring checks are relaxed for tests and `__init__.py`.
- Use snake_case for functions/variables, CapWords for classes, and UPPER_SNAKE_CASE for constants.
- Keep modules focused; place CLI-only code in `gopher/gopher.py` and reusable logic elsewhere.

## Testing Guidelines
- Write fast, deterministic tests; mock network/filesystem when practical.
- Place new unit tests beside related code in `tests/unit_tests/`; heavier end-to-end flows go in `tests/system_tests/`.
- Name tests descriptively (`test_handles_missing_annotations`), and prefer given/when/then structure in assertions.
- Aim to cover new branches and edge cases before opening a PR; run `pytest -q` to keep output compact.

## Commit & Pull Request Guidelines
- Commit messages in this repo are short, imperative phrases (e.g., `Update pre-commit hooks`, `Migrate to uv and ruff`). Match that style.
- Keep commits focused; include only related changes plus updated tests.
- PRs should state the problem, the solution, and validation (tests, screenshots for visual changes). Link related issues and note breaking changes explicitly.
- Ensure lint, tests, and (if touched) docs build succeed before requesting review.

## Security & Configuration Notes
- Avoid committing secrets or private keys; pre-commit includes a detector, but double-check before pushing.
- Network calls (e.g., ontology downloads) should remain in well-scoped helper functions in `gopher/` so they can be mocked in tests.
