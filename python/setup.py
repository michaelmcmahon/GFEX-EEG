"""GFEX-EEG: Python package manifest (setuptools)

Installs ``geodesic_rescue_py`` — the Python implementation of the
GFEX-EEG toolbox for MRI-free EEG source-imaging fiducial recovery.

See also CITATION.cff in the repository root (machine-readable).
Repository: https://github.com/michaelmcmahon/GFEX-EEG
License: MIT — see LICENSE in the repository root.
"""
# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2026 Michael McMahon, University of Galway

from setuptools import setup, find_packages

setup(
    name="geodesic-rescue-py",
    version="1.1.5",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "mne",
    ],
    package_data={
        "geodesic_rescue_py": [
            "data/ICBM152_scalp.mat",
            "data/weight_zoo.json",
            "data/mlp/*.mat",
        ],
    },
    description="Python implementation of the Geodesic Rescue Toolbox for EEG fiducial extrapolation.",
    long_description=open("README.md").read() if os.path.exists("README.md") else "",
    long_description_content_type="text/markdown",
    python_requires=">=3.7",
)
