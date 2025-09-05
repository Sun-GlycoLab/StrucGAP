# StrucGAP

**Structural and site-specific Glycoproteomics Analysis Platform**

**StrucGAP** is a modular toolkit for downstream analysis of glycoproteomics data. It includes support for preprocessing, glycan structure analysis, quantification, network visualization, functional annotation, and more.

View the [**documentation**](https://strucgap.readthedocs.io/en/latest/index.html) for more information.

## Installation

```bash
pip install strucgap

## Troubleshooting

StrucGAP works successfully on Python **3.7, 3.8, 3.9, 3.10**,  
but for stability we recommend **Python 3.9 or 3.10**.

>  Most errors are caused by incompatible versions of **NumPy** and **pandas**.

- For **Python 3.7 and 3.8** → use **NumPy 1.18.1**  
- For **Python 3.9 and 3.10** → use **NumPy 1.26.4**  
- For all Python versions → use **pandas 1.3.5**

---

### Q1: `ModuleNotFoundError: No module named 'statsmodels'`

**A1:**
```bash
pip install statsmodels
```

### Q2: `ImportError: Missing optional dependency 'openpyxl'`

**A2:**
```bash
pip install openpyxl
```

### Q3: `TypeError: loop of ufunc does not support argument 0 of type float which has no callable log2 method`

**A3:**
```bash
pip install numpy==1.26.4
pip install pandas==1.3.5
```

### Q4: `ModuleNotFoundError: No module named 'importlib.metadata'`
### Cause: This package is included in the Python standard library only since version 3.8.

**A4:**
```bash
pip install importlib-metadata==6.7.0
```

### If your problem is not listed here, please open a GitHub Issue with your error message and environment details (Python version, OS, and package versions).




