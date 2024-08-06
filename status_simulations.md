the data comes from 3 runs: gm2 window, tau, run3

# correlator with same mu

ens       |  mul    |  mus          |    muc                    |  mul1     | sea-loop-l 
----------|---------|---------------|---------------------------|-----------|--
cA.53.24  | 0.00530 | 0.010  0.020  | 0.265  0.290  0.300       |           |
cA.40.24  | 0.00400 | 0.010  0.020  | 0.265  0.290  0.300       |           |
cA.30.32  | 0.00300 | 0.010  0.020  | 0.265  0.290  0.300       |           |
cB.72.64  | 0.00072 | 0.018  0.020  | 0.210  0.230  0.250       | 0.0006675 | yes   
cB.72.96  | 0.00072 | 0.018  0.019  | 0.210  0.230  0.250       | 0.0006675 | yes
cC.06.80  | 0.00060 | 0.016  0.018  | 0.175  0.195  0.215       | 0.0005850 | yes
cC.06.112 | 0.00060 | 0.016  0.018  |                           | 
cD.54.96  | 0.00054 | 0.013  0.014  | 0.165  0.175  0.175       | 0.0004964 | yes
cE.44.112 | 0.00044 | 0.011  0.012  | 0.130  0.140  0.150       | 4.3100e-04| ??

correlator of mus comes from tau run
muc from gm2_window
mul from gm2_window

# Mesons



ens       |  Pi   |  K                 |    D                |  Ds
----------|-------|--------------------|---------------------|--- 
cA.53.24  | l1    |     no             |                     |
cA.40.24  | l1    |     no             |                     |
cA.30.32  | l1    |     no             |                     |
cB.72.64  | l1 l2 | l1-(s1-s2)         |                     | 
cB.72.96  | l1 l2 | l1-(s1-s2)         |                     | 
cC.06.80  | l1 l2 | l1-(s1-s2)         |                     | 
cC.06.112 | l1    | ???                | ?????               | 
cD.54.96  | l1 l2 | (l1)-(s1-s2)  (l1)-sa     | (l1)-(c1-c2-c3)    | sa-(c1-c2-c3)
cE.44.112 | l1 l2 | (l1)-(s1-s2)  (l1-l2)-sa  | (l1-l2)-(c1-c2-c3) | sa-(c1-c2-c3)

l1-(s1-s2)    from tau run
(l1-l2)-sa    comes from run3 cE.44.112
l2 for cB.72.64, cB.72.96, cC.06.80, cD.54.96  from gm2_window
l2 for cE.44.112  from run3 
l2 for cD.54.96  from run3 
sa-(c1-c2-c3)  sa is a stange value, opefully the physical one, we need a new strange sb-(c1-c2-c3)


ens        |  sa
-|-
cD.54.96 |   1.3557e-02
cE.44.112 |  1.1759e-02
