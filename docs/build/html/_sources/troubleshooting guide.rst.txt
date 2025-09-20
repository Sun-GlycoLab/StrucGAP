Troubleshooting guide
========================
StrucGAP successfully on Python 3.7, 3.8, 3.9, 3.10, but for stability we recommend **Python 3.9 or 3.10**.

It's worth noting that most errors are caused by the versions of **NumPy** and **pandas**.

| For Python 3.7 and 3.8, we strongly recommend using **NumPy 1.18.1**.
| For Python 3.9 and 3.10, we strongly recommend using **NumPy 1.26.4**.
| For all Python versions, we recommend using **pandas 1.3.5**.

Q1: ModuleNotFoundError: No module named 'statsmodels'

A1:

.. code-block:: bash

    pip install statsmodels


Q2: ImportError: Missing optional dependency 'openpyxl'

A2:

.. code-block:: bash

    pip install openpyxl


Q3: TypeError: loop of ufunc does not support argument 0 of type float which has no callable log2 method

A3:

.. code-block:: bash

    pip install numpy==1.26.4
    pip install pandas==1.3.5


Q4: ModuleNotFoundError: No module named 'importlib.metadata'

Cause: This package is included in the Python standard library only since version 3.8.

A4:
Upgrade Python to ≥3.8 (recommended), or

.. code-block:: bash

    pip install importlib-metadata==6.7.0


Q5: During installation, I see errors mentioning Rust or Cargo. Do StrucGAP or gseapy depend on Rust?

Cause: Neither StrucGAP nor gseapy requires Rust directly. However, in some environments pip may try to compile dependencies from source (instead of using pre-built wheels). When that happens, Rust may be required if the dependency’s source build involves it.

A5:
Alternatively, ensure that your environment has access to the official Python package index (PyPI) so wheels can be downloaded. If compilation still occurs, you may need to install Rust temporarily, but this is not normally required.

If your problem is not listed here, please open a GitHub Issue with your error message and environment details (Python version, OS, and package versions).
------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



