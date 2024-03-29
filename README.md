# Transcriptomic analysis of the SmLE-PZQ-ER and SmLE-PZQ-ES populations

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5297217.svg)](https://doi.org/10.5281/zenodo.5297217)

This repository contains the notebook and scripts used to analyze the gene and isoform expression of juvenile and adult worms from the SmLE-PZQ-ER and SmLE-PZQ-ES populations in order to better understand the origin of the praziquantel (PZQ) resistance. SmLE-PZQ-ER and ES are specific populations enriched in resistant (ER) or sensitive (ES) alleles to PZQ. This analysis is part of the work described in [Genetic analysis of praziquantel response in schistosome parasites implicates a Transient Receptor Potential channel](https://doi.org/10.1101/2021.06.09.447779). 

After the initial publication in 2021, an improved reference genome (version 10), along with a revised annotation, was published in 2023. We reanalyzed our transcriptomic data using this improved version and updated the present repository. The revised analysis is part of the short note [A single locus determines praziquantel response in *Schistosoma mansoni*](https://doi.org/10.1101/2023.11.01.565202).

## Abstract

Mass treatment with praziquantel (PZQ) monotherapy is the mainstay for schistosomiasis treatment. This drug shows imperfect cure rates in the field and parasites showing reduced PZQ response can be selected in the laboratory, but the extent of resistance in *Schistosoma mansoni* populations is unknown. We examined the genetic basis of variation in PZQ response in a *S. mansoni* population (SmLE-PZQ-R) selected with PZQ in the laboratory: 35% of these worms survive high dose (73 µg/mL) PZQ treatment. We used genome wide association to map loci underlying PZQ response. The major chr. 3 peak contains a transient receptor potential (*Sm.TRPM<sub>PZQ</sub>*) channel (Smp\_246790), activated by nanomolar concentrations of PZQ. PZQ response shows recessive inheritance and marker-assisted selection of parasites at a single *Sm.TRPM<sub>PZQ</sub>* SNP enriched populations of PZQ-resistant (PZQ-ER) and sensitive (PZQ-ES) parasites showing >377 fold difference in PZQ response. The PZQ-ER parasites survived treatment in rodents better than PZQ-ES.  Resistant parasites show 2.25-fold lower expression of *Sm.TRPM<sub>PZQ</sub>* than sensitive parasites. Specific chemical blockers of *Sm.TRPM<sub>PZQ</sub>* enhanced PZQ resistance, while *Sm.TRPM<sub>PZQ</sub>* activators increased sensitivity. A single SNP in *Sm.TRPM<sub>PZQ</sub>* differentiated PZQ-ER and PZQ-ES lines, but mutagenesis showed this was not involved in PZQ-response, suggesting linked regulatory changes. We surveyed *Sm.TRPM<sub>PZQ</sub>* sequence variation in 259 parasites from the New and Old World revealing one nonsense mutation that results in a truncated protein with no PZQ-binding site. Our results demonstrate that *Sm.TRPM<sub>PZQ</sub>* underlies variation in PZQ response in *S. mansoni* and provides an approach for monitoring emerging PZQ-resistance alleles in schistosome elimination programs.

## Prerequisites

Two dependencies must be installed before running the Jupyter notebook:
* A [conda](https://docs.conda.io/en/latest/) distribution, like [miniconda](https://docs.conda.io/en/latest/miniconda.html),
* [Jupyter Notebook](https://jupyter.readthedocs.io/en/latest/install.html) with a [bash kernel](https://github.com/takluyver/bash_kernel).

Once this is done, run the first cells of the Jupyter notebook from this repository to install a dedicated conda environment and related R packages.
