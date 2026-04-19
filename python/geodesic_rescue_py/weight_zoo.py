# /*******************************************************************************
# * GFEX-EEG - Geodesic fiducial extrapolation for MRI-free EEG source imaging
# * Version: 1.0.0
# * Repository: https://github.com/michaelmcmahon/GFEX-EEG
# * License:  MIT License
# * Authors: Michael McMahon / University of Galway
# *******************************************************************************/

"""Weight-zoo loader for GFEX-EEG cohort-specific FNO presets."""

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
