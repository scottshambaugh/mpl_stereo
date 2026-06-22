import matplotlib

# Use the non-interactive Agg backend so the tests run headless (e.g. on CI,
# where no display / a broken Tk install would otherwise crash figure creation).
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest


@pytest.fixture(autouse=True)
def _close_figures():
    """Close any figures a test opened, to avoid leaking them across tests."""
    yield
    plt.close("all")
