const MAX_CONTRACTION = 12
export MAX_CONTRACTION

#                        | s | p |    d     |    f    |
#                        ------------------------------ 
const axial_norm_fact = [ 1.0 1.0    1.0        1.0   ; 
                          0.0 1.0 sqrt(3.0)  sqrt(5.0);
                          0.0 1.0 sqrt(3.0)  sqrt(5.0);
                          0.0 0.0    1.0     sqrt(5.0);
                          0.0 0.0 sqrt(3.0) sqrt(15.0);
                          0.0 0.0    1.0     sqrt(5.0);
                          0.0 0.0    0.0        1.0   ;
                          0.0 0.0    0.0     sqrt(5.0);
                          0.0 0.0    0.0     sqrt(5.0);
                          0.0 0.0    0.0        1.0
                        ]
export axial_norm_fact

