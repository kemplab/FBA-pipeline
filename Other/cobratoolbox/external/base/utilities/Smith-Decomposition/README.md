# Smith-Decomposition

The Smith normal form (also called Smith Canonical form or Invariant Factor theorem) is a diagonal matrix D that contains the invariant factors of any A matrix of size n × m over a field F (in the attached implementation it is provided for the ring of integers Z and rings of polynomials F[x]).

```
 D = |d1 0 ... 0 ... 0|= TAS
     |0 d2 ... 0 ... 0| 
     |: : ...  : ... :| 
     |0 0 ... dr ... 0| 
     |: : ...  : ... :| 
     |0 0 ...  0 ... 0|
```

where d1 , ..., dr ∈ F are monic, dj |dj+1 for 1 ≤ k ≤ r − 1. T is a product of elementary row unimodular matrices, and S is a product of elementary column unimodular matrices.

Provided are two functions: 
```
[T,D,S]=smithFormInt(A) 
```
for integer matrices and: 
```
[T,D,S]=smithFormPoly(A) 
```
for polynomial matrices.

Run the 
```
smithDemo.m 
```
for usage and example matrices.

*Author:* Nadia Figueroa
