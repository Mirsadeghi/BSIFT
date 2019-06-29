# BSIFT
BSIFT: Boosting SIFT Using Principal Component Analysis. M. Fotouhi, S. E. Mirsadeghi, S. Kasaei, K. Faez

This repository contains implementation of the above paper.

# Quick Start
Run BSIFT_demo.m to see result of BSIFT on the leuven dataset.

# Dependency
VLFeat -- Vision Lab Features Library 
https://github.com/vlfeat/vlfeat

# Results
TABLE III. 	TOTAL AND CORRECT MATCHES FOR GRAFFITI SAMPLE </br>/
| Method | PCA Energy | Equivalent Dimension | I1 -> I2 | I1 -> I3 | I1 -> I4 | I1 -> I5 | I1 -> I6 |
| ------ | ---------- | -------------------- | -------- | -------- | -------- | -------- | -------- |
| BSIFT	 | 50 %	|               8	           |  957/1384| 305/725  |	50/380  | 9/298    |	0/197   |
| BSIFT	 | 60 %	|               13           |  986/1473| 367/882  |	74/475  | 10/366   |	3/293   |
| BSIFT	 | 70 %	|               19           | 1005/1533| 421/995  |	94/572  | 12/441   |	3/350   |
| BSIFT	 | 80 %	|               30           | 1010/1554| 441/1039 |	105/609 | 15/468   |	3/390   |
| BSIFT	 | 90 %	|               51           | 1014/1568| 448/1050 |	111/663 | 18/491   |	2/414   |
| SIFT	 | 100 %|               128          |  997/1505| 417/1001 |	98/612  | 14/501   |	2/438   |
