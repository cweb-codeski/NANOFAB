"""Functional check of the reduced linearized-BTE law against vaartstra2019_thesis Fig. 7 (a = 0.052).

This tests our implementation against the authors' own plotted curve at one accommodation value.
It is not experimental validation and does not verify the Meland-Ytrehus form at other values of a.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from nanofab.evaporation import SONE_C4, schrage_heat_flux, sone_bte_heat_flux
from nanofab.fluids import UNIVERSAL_GAS_CONSTANT, Fluid

CSV = Path(__file__).resolve().parents[1] / "data" / "literature" / "vaartstra2019_thesis_fig7_digitized.csv"
T_SAT = 373.15
A = 0.052


@pytest.fixture(scope="module")
def fig7():
    data = np.loadtxt(CSV, delimiter=",", comments="#", skiprows=2)
    return data[(data[:, 0] > 1.5) & (data[:, 0] < 9.8)]


@pytest.fixture(scope="module")
def water():
    return Fluid("Water")


def test_reduced_bte_matches_thesis_fig7_within_digitization_error(fig7, water):
    """Reading error is ~1% of the 500 W/cm^2 axis span (line thickness), i.e. ~5 W/cm^2 absolute."""
    model = np.array([sone_bte_heat_flux(T_SAT + dT, T_SAT, A, water) / 1e4 for dT in fig7[:, 0]])
    assert np.max(np.abs(model - fig7[:, 1])) < 5.0
    upper = fig7[:, 0] > 3.0
    assert np.max(np.abs(model[upper] / fig7[upper, 1] - 1)) < 0.01


def test_alternative_accommodation_forms_do_not_match(fig7, water):
    """C4/a, or evaluating the Knudsen-layer density at the far-field pressure, miss the curve badly."""
    R_s = UNIVERSAL_GAS_CONSTANT / water.molar_mass
    dT = fig7[:, 0]
    T = T_SAT + dT
    dps = np.array([water.p_sat(t) - water.p_sat(T_SAT) for t in T])
    hlv = np.array([water.h_lv(t) for t in T])
    base = dps * np.sqrt(2 / (R_s * T)) * hlv / 1e4
    c4_over_a = base / (SONE_C4 / A)
    far_field = base / (SONE_C4 + 2 * math.sqrt(math.pi) * (1 - A) / A) * (water.p_sat(T_SAT) / (water.p_sat(T_SAT) + dps))
    assert np.median(c4_over_a / fig7[:, 1]) > 1.4
    assert np.median(far_field / fig7[:, 1]) < 0.9


def test_schrage_at_same_a_is_a_few_percent_lower(fig7, water):
    i = np.argmin(abs(fig7[:, 0] - 8.0))
    sch = schrage_heat_flux(T_SAT + fig7[i, 0], T_SAT, A, water) / 1e4
    bte = sone_bte_heat_flux(T_SAT + fig7[i, 0], T_SAT, A, water) / 1e4
    assert 0.95 < sch / bte < 0.99
