# RILLIE: RNA In Silico Evolution via LLM and Inverse folding

✨ Welcome to the official repository for "RILLIE: RNA In Silico Evolution via LLM and Inverse folding". 

-- 🔗 Read our paper: [ArXiv](https://arxiv.org/abs/****.*****) --

## Overview

🚀 We introduce RILLIE, a general RNA foundation model that integrate sequence and structure information to evolve RNA in a zero-shot fashion. RILLIE integrating a large language model with an inverse folding model can generate functional RNA sequences aligning with natural evolutionary patterns at the sequence level while preserving the structural integrity of key functional regions. Using RILLIE, we successfully evolved two engineered RNA aptamers, Broccoli and Pepper, with a high success rate, low sequence similarity, improved binding affinity and fluorescence in live cell.

![RILLIE Overview](./image/RILLIE.png)

## 🗂 Contents

- [Installation](#installation)
- [Benchmark](#Benchmark)
  - [Dataset Description](#dataset-description)
  - [Model Description](#model-description)
  - [Generate RNA Structure](#generate-rna-structure)
  - [Zero-shot ncRNA fitness prediction](#data-processing)
  - [High-fitness sensitivity analysis](#data-processing)
- [Using RILLIE](#using-RILLIE)
  - [CheckPoints](#checkpoints)
  - [General RNA Evolution](#general-rna-evolution)
  - [Multi-Round Evolution (Optional)](#multi-round-evolution-(Optional))
<!-- - [Citation](#citation) -->
- [Acknowledgements](#acknowledgements)

## Installation
### If you prefer a faster setup, you can use the provided RILLIE.yaml file:
```bash

conda env create -f RILLIE.yaml -y
conda activate RILLIE

```

### You can also install the environment either by following the step-by-step instructions below.
```bash
# Create a conda environment
conda create -y -n RILLIE python=3.9
conda activate RILLIE

# Install PyTorch and CUDA dependencies
#conda install pytorch==2.3.0 torchvision==0.18.0 torchaudio==2.3.0 pytorch-cuda=12.1 -c pytorch -c nvidia
pip install torch==2.3.0+cu121 torchvision==0.18.0+cu121 torchaudio==2.3.0+cu121

# Install other dependencies
pip install rna-fm==0.2.2
pip install ml-collections==0.1.1
pip install simtk==0.1.0
pip install openmm==8.1.1
pip install torchdrug==0.2.1
pip install torch-geometric==2.5.3
pip install equiformer-pytorch
pip install rdkit-pypi==2021.3.4
pip install biopython==1.84
pip install pandas==2.2.2
pip install scikit-learn==1.5.1
pip install prody==2.4.1

```

<!--```bash
if you have problem like this error information: ImportError: cannot import name 'packaging' from 'pkg_resources' (/xcfhome/ypxia/anaconda3/envs/gernabind/lib/python3.8/site-packages/pkg_resources/__init__.py) with torchdrug, you can try this command below:
pip install --upgrade packaging
pip install --upgrade pip setuptools
```
-->

## Benchmark
### Dataset Description
We collected 6 ncRNA DMS datasets including tRNA, RNA aptamer and ribozyme from previous papers or private data:
- tRNA datasets:
  
   [Pair-tRNA](https://www.nature.com/articles/s41586-018-0170-7)

   [tRNA](http://genesdev.cshlp.org/lookup/doi/10.1101/gad.245936.114)
- RNA aptamer datasets:

   [Pepper](https://www.nature.com/articles/s41587-019-0249-1)

   [Clivia](https://www.nature.com/articles/s41592-023-01997-7)

    [Okra](https://www.nature.com/articles/s41589-024-01629-x)
- ribozyme datasets:

   [Glms](https://www.nature.com/articles/s41467-020-15540-1)
### Model Description
Our benchmark includes following models:
- RNA language models:

  [AIDO.RNA(1.6B)](https://github.com/genbio-ai/AIDO)

  [RiNALMo](https://github.com/lbcb-sci/RiNALMo)

  [RNAFM](https://github.com/ml4bio/RNA-FM)

  [RNAMSM](https://pmc.ncbi.nlm.nih.gov/articles/PMC9477273/)

- DNA language models:

  [Evo 1](https://hariboss.pasteur.cloud/)

  [Nucleotide Transformer](https://hariboss.pasteur.cloud/)

  [Grover](https://hariboss.pasteur.cloud/)

  [GENA](https://hariboss.pasteur.cloud/)

- RNA inverse-folding models:
  
  [RhoDesign](https://github.com/ml4bio/RhoDesign)

### Generate RNA Structure
We use [RhoFold](https://github.com/ml4bio/RhoFold) and [AlphaFold3](https://github.com/google-deepmind/alphafold3) to generate RNA 3D Structure 

We use [Chai](https://github.com/chaidiscovery/chai-lab) to generate RNA 2D structure (as the input of RILLIE).

### Data Processing
You can process data through the following steps:
```bash
python data_utils/process_data.py --fasta example/a.fasta --smile example/mol.txt --RhoFold_path your_RhoFold_project_path --RhoFold_weight RhoFold_model_weight_path
```
And the processed data will be saved in ./data folder as "new_data.pkl" file.

### Dataset Download
We process the processed data with Robin & Biosensor dataset. You can download the processed data from [Zenodo](https://zenodo.org/records/14808549).

## Using RILLIE

### CheckPoints
Download the model weights and put into the "Model" folder, which contains the model checkpoint. You can direct run the scripts in ./Model folder to ger the model weights.

```bash
bash Model/get_weights.sh
```

### General RNA Evolution
Otherwise, you can also use our model to get RNA target binding sites prediction.
You can run the file below, so that you can get the RNA_binding.csv about RNA.

```bash
python inference_binding_site.py
```
### Multi-Round Evolution (Optional)
Based on the wet-lab testing result, we can discarding harmful mutations to increase the success rate, while introducing new mutations to help direct evolution escape local optima and discover global optima.  This approach enables efficient directed evolution without retraining the model and is specially useful when the tested variants are very few.
```bash
python train_model.py --dataset Robin --split_method random --model_output_path Model/
```

## License
No Commercial use of either the model nor generated data, details to be found in license.md.

<!-- ## Citation
🤗 If you find our code useful or our work relevant, please consider citing:

```
@article{gerna2025,
  title={GerNA-Bind: Geometric-informed RNA-ligand Binding Specificity Prediction with Deep Learning},
  author={Yunpeng Xia, Jiayi Li, Jiahua Rao, Dong-Jun Yu, Xiucai Chen and Shuangjia Zheng},
  journal={arXiv preprint arXiv:****.*****},
  year={2025}
}
``` -->

## Acknowledgements

Our work builds upon [AIDO.RNA(1.6B)](https://github.com/genbio-ai/AIDO),[RiNALMo](https://github.com/lbcb-sci/RiNALMo),[RNAFM](https://github.com/ml4bio/RNA-FM),[RNAMSM](https://pmc.ncbi.nlm.nih.gov/articles/PMC9477273/),[Evo 1](https://hariboss.pasteur.cloud/),[Nucleotide Transformer](https://hariboss.pasteur.cloud/),[Grover](https://hariboss.pasteur.cloud/),[GENA](https://hariboss.pasteur.cloud/) ,[RhoDesign](https://github.com/ml4bio/RhoDesign),[RhoFold](https://github.com/ml4bio/RhoFold)Thanks for their excellent work and open-source contributions.
