"""GFEX-EEG: Cohort-specific weight-zoo loader (Python)

Retrieves cohort-specific FNO weight presets from ``weight_zoo.json``.
Each preset pins the geometric hyperparameters (rho, beta, D_standard)
and the path to the cohort-specific Tier 1.5 MLP weights file.

Three production presets ship in V1.1.5: ``LEMON_Polhemus_Adult``,
``WH_Neuromag70_Adult``, ``CapTrak_Adult``. The ``default`` and ``lemon``
aliases map to ``LEMON_Polhemus_Adult``.

Usage
-----
    from geodesic_rescue_py import load_cohort_preset
    preset = load_cohort_preset('LEMON_Polhemus_Adult')
    # -> {'tag': ..., 'rho': 0.248383, 'beta': 0.235926, 'D_standard': 0.1388, ...}

Citation
--------
If you use this code in research, please cite both the software archive
and the accompanying manuscript:

    [Software]
    McMahon, M., Schukat, M., & Barrett, E. (2026).
    GFEX-EEG Toolbox [Software].
    Zenodo. https://doi.org/10.5281/zenodo.20580899

    [Paper]
    McMahon, M., Schukat, M., & Barrett, E. (Submitted).
    GFEX-EEG: Geodesic recovery of anatomical fiducials for MRI-free
    EEG source imaging.

See also CITATION.cff in the repository root (machine-readable).
Current version: see ``geodesic_rescue_py.__version__`` (single source
of truth in ``__init__.py``).

Repository
----------
https://github.com/michaelmcmahon/GFEX-EEG

License
-------
MIT — see LICENSE in the repository root.
"""
# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2026 Michael McMahon, University of Galway

import json
import os

_PACKAGE_DATA_ZOO = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'data', 'weight_zoo.json'
)


def load_cohort_preset(tag, zoo_path=None):
    """Retrieve a cohort-specific FNO weight preset from the zoo.

    Parameters
    ----------
    tag : str
        Cohort tag (e.g. ``'LEMON_Polhemus_Adult'``, ``'lemon'``,
        ``'default'``). Aliases are resolved transparently, case-insensitive.
    zoo_path : str, optional
        Path to ``weight_zoo.json``. Defaults to the JSON bundled inside the
        ``geodesic_rescue_py`` package under ``data/weight_zoo.json``.

    Returns
    -------
    dict
        Dict with keys: ``tag``, ``rho``, ``beta``, ``D_standard``,
        ``status``, and ``description`` when available.

    Raises
    ------
    FileNotFoundError
        If ``weight_zoo.json`` is not located.
    KeyError
        If the tag is not in the zoo (listing available production tags).
    RuntimeError
        If the tag resolves to an entry with ``status == 'pending_tune'``.
    """
    if zoo_path is None:
        zoo_path = _PACKAGE_DATA_ZOO
    if not os.path.isfile(zoo_path):
        raise FileNotFoundError(f"weight_zoo.json not found at {zoo_path}")

    with open(zoo_path, 'r', encoding='utf-8') as fh:
        zoo = json.load(fh)

    # Resolve aliases (case-insensitive)
    resolved = tag
    aliases = zoo.get('aliases', {})
    for alias_name, alias_target in aliases.items():
        if alias_name.lower() == tag.lower():
            resolved = alias_target
            break

    weight_sets = zoo.get('weight_sets', {})
    if resolved not in weight_sets:
        production_tags = [
            k for k, v in weight_sets.items()
            if v.get('status') == 'production'
        ]
        raise KeyError(
            f'Cohort tag "{tag}" not found in weight zoo.\n'
            f'Available production tags: {", ".join(production_tags)}\n'
            f'Available aliases: see weight_zoo.json'
        )

    entry = weight_sets[resolved]

    if entry.get('status') == 'pending_tune':
        note = entry.get('note', '')
        raise RuntimeError(
            f'Cohort "{resolved}" is marked pending_tune '
            f'(no production weights yet).\n'
            f'Note: {note}\n'
            f'Use Tier 2 FNO tuner on pilot data with MRI HTJ ground truth, '
            f'or fall back to the default LEMON preset.'
        )

    return {
        'tag': resolved,
        'rho': entry['rho'],
        'beta': entry['beta'],
        'D_standard': entry['D_standard'],
        'status': entry['status'],
        'description': entry.get('description'),
        'mlp_weights_file': entry.get('mlp_weights_file'),
        'mlp_holdout_mean_err_mm': entry.get('mlp_holdout_mean_err_mm'),
    }
