Source code for my PhD thesis
=============================

Setting up
----------

* Create a python virtual environment: ``python3.9 -m venv .venv && source .venv/bin/activate``
* Install NMR-EsPy: ``pip install --upgrade pip && pip install nmrespy``
* Install the ``utils`` package: ``pip install -e ./utils``
* Configure the ``matplotlibrc`` file: ``./script/ud_mplrc.sh``

Compilation of the thesis can then be done by running ``./scripts/compile.sh``

Figures can be produced by running the relevant ``make_figure.py``
script; these are found within the ``./figures/`` directory. *Always run these
from the root directory*.
