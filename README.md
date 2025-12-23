# DSAEA-PS

**Dual Surrogate-Assisted Evolutionary Algorithm Based on Parallel Search  
for Expensive Multi/Many-Objective Optimization**

---

## 📌 Introduction

This repository provides the **official implementation** of **DSAEA-PS**, a dual surrogate-assisted evolutionary algorithm based on parallel search, proposed for solving **computationally expensive multi-objective and many-objective optimization problems**.

DSAEA-PS simultaneously leverages **approximation-based** and **classification-based** surrogate models, enabling cooperative utilization of **solution quality** and **uncertainty information**. Moreover, a **parallel search framework** built upon **heterogeneous multi-objective evolutionary algorithms (MOEAs)** is introduced to enhance exploration capability in complex decision spaces.

The algorithm has been validated on widely-used benchmark problems and a real-world **five-objective blended-wing-body underwater glider design problem** involving time-consuming CFD and structural simulations.

---

## 📄 Reference

If you use this code in your research, please cite the following paper:

> Jiangtao Shen, Peng Wang, Ye Tian, Huachao Dong,  
> **A dual surrogate assisted evolutionary algorithm based on parallel search for expensive multi/many-objective optimization**,  
> *Applied Soft Computing*, 2023.  
> DOI: https://doi.org/10.1016/j.asoc.2023.110879

### BibTeX
```bibtex
@article{Shen2023DSAEAPS,
  title   = {A dual surrogate assisted evolutionary algorithm based on parallel search for expensive multi/many-objective optimization},
  author  = {Shen, Jiangtao and Wang, Peng and Tian, Ye and Dong, Huachao},
  journal = {Applied Soft Computing},
  year    = {2023},
  pages   = {110879},
  doi     = {10.1016/j.asoc.2023.110879}
}
