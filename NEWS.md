### Release 0.24 (06Dec2022)

Allow N and other IUPAC codes that might be in a primer to Y and R, specifically 2 base ones W, S, M, K


### Release 0.23 (19Oct2022)

Fix problem with multiple sample named files


### Release 0.22 (31Mar2020)

Add 0 method to identify sample from even a forward barcode (any port in a storm)


### Release 0.21 (03Dec2018)

Sensitive to use of same barcode/index in forward and reverse lists:
* Method 1 allows use of same index in both lists
* Method 2 reports an error for this
* Method 3 reports a warning, but continues

Add **-info all** and **-info both** options:
* **-info both** combines the forward and reverse indexes and reports edit distance for all
* **-info all** shows primer edit distances and both index set edit distances

Improved error handling for sample file

### Release 0.19 (12Jun2018)

   Added minibar.py as release 0.19 to github 12Jun2018
