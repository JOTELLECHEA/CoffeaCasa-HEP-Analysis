# CoffeaCasa-HEP-Analysis
<img src="https://github.com/JOTELLECHEA/CoffeaCasa-HEP-Analysis/blob/main/title.png" width="100%">

This project will calulate the ratio of $W^+$ and $W^-$. The focus will be on leptonic decays from either $W^+$ and $W^-$.
The W Bosons decay as follows (l = electron or muon): $W^- \rightarrow l^-v$ and $W^+ \rightarrow l^+\bar{v}$

Steps that will take:

1. Calculating the Z Boson Invariant mass. (This is to enusre we are accesing the data correctly).
2. Find a data set on CMS Opendata for $W^{+/-}$ canidates.
3. Calculate the invariant mass of $W^{+/-}$.
4. Calculate the ratios of $W^{+/-}$.

## Feynman Diagrams (electron) 

### $W^+ \rightarrow e^+\bar{v}_e$
<img src="https://github.com/JOTELLECHEA/CoffeaCasa-HEP-Analysis/blob/main/Pictures/wplus.jpg" alt="wplus" width="300"/>

### $W^- \rightarrow e^-v_e$
<img src="https://github.com/JOTELLECHEA/CoffeaCasa-HEP-Analysis/blob/main/Pictures/wminus.jpg" alt="wminus" width="300"/>

## How to use Analysis notebook

1. Need to have access to https://coffea-opendata.casa/
2. ``` git clone https://github.com/JOTELLECHEA/CoffeaCasa-HEP-Analysis.git ```
3.  ```cd /CoffeaCasa-HEP-Analysis/Analysis ```
4. Open the Jupyter notebook to run the analysis.
5. The ```functions.py``` file has all the functions needed to run the analysis. 
