[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Black-Box CoVaR and Its Gradient Estimation

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
[Black-Box CoVaR and Its Gradient Estimation](https://doi.org/10.1287/ijoc.2024.0833) by Hao Cao, Jian-Qiang Hu, and Jiaqiao Hu. 

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2024.0833

https://doi.org/10.1287/ijoc.2024.0833.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{cao2025,
  author =        {Hao Cao and Jian-Qiang Hu and Jiaqiao Hu},
  publisher =     {INFORMS Journal on Computing},
  title =         {{Black-Box CoVaR and Its Gradient Estimation}},
  year =          {2025},
  doi =           {10.1287/ijoc.2024.0833.cd},
  url =           {https://github.com/INFORMSJoC/2024.0833},
  note =          {Available for download at https://github.com/INFORMSJoC/2024.0833},
}  
```

## Description

The goal of this software is to demonstrate a stochastic approximation algorithm 
for CoVaR and its gradient estimation, which is applicable in various complicated cases, 
including black-box and streaming data scenarios.

## Structure

The structure of this repository is as follows:
- `KBSA.py`: The code of the proposed kernel-based stochastic approximation (KBSA) algorithm.
- `implementation/`: The directory containing implementation examples.
- 'results/': The directory storing results for independent replications of the algorithm.

## Requirements
The code is tested in the environment of python 3.12.4 with Windows 11.  
Numerous package requirements exist, including but not limited to: 
`numpy`, `pandas`, `matplotlib`, `scipy`, `operator`, `IPython`, 
`random`, `itertools`, `math`, `pickle`, `joblib`, `seaborn`, and `statsmodels`.

To meet GitHub's file size limit, in the `results/` directory, dataset files, 
e.g., `toy3.pkl` and `toy3_20d.pkl`, are compressed; please download and extract  
the files into the `results/` directory, before previewing.
