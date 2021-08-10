# EvoCov
<img align="right" src="evocado.png" width=150px>
This is EvoCov, a pipeline designed for analysis of SARS-CoV-2 sequences from GISAID. The pipeline can be run interactively or by default, with a view to using SARS-CoV-2 sequence information to make an evolutionarily aware estimate of efficient epitopes on the spike protein for antibody design. While the code is open source, the data intended for use with the pipeline may only be obtained with express permission from [GISAID](https://www.gisaid.org/).

## Installation
Clone the github repository to your machine to use the EvoCov package. Before using, you should check that any python dependencies are installed.

```bash
git clone https://github.com/ciarajudge/EvoCov.git
```

```bash
pip install -r requirements.txt
```

## Preparation of site-wise mutation rates using Treecov
The prediction aspect of the pipeline makes use of estimated site-wise mutation rates from analysis of a SARS-CoV-2 phylogenetic tree with baseml, a [Phylogenetic Analysis by Maximum Likelihood](http://abacus.gene.ucl.ac.uk/software/paml.html) program. To generate these rates, you must download and compile paml and place it in the ./treecov/ directory, where the path to the baseml executable is ./treecov/paml/bin/baseml. It is important that the folders are named correctly. You must also download a phylogenetic tree on GISAID by clicking Audacity on the platform, and place the file global.tree in the ./treecov/ directory. To run the treecov pipeline to generate the rates, navigate to the treecov directory and use the command:
```bash
python treepipe.py /absolute/path/to/GISAID/fasta/file
```
This initiates the process of iterative sampling and analysis of the phylogenetic tree 100 times, in 10 batches of 10. These batch sizes, or the number of batches, can be adjusted by changing the number of loops in the code in subsampletree.R (for batch size) and treepipe.py (for no. of batches).

## Default Usage of Evocov
Navigate to the cloned repository and call the package along with the file paths of your latest GISAID unmasked sequence file and metadata file. This will initiate a default run of the pipeline, including handling of any exceptions or options. This includes the final step of the pipeline where the results are piped to a PDF using R. 

```bash
python -m evocov /path/to/sequencefile_masked.fa /path/to/metadata.tsv
```

If you'd like to be notified by text when the pipeline is complete, pass a third argument with a valid mobile number (no plus signs or brackets) for example: 353877910680 where the country code is +353 and the phone number is 0877910680.

```bash
python -m evocov /path/to/sequencefile_masked.fa /path/to/metadata.tsv 353877910680
```

## Interactive Usage
Navigate to the cloned repository and call the package using the below command.

```bash
python -m evocov
```

Running the pipeline in this manner will create an interactive session where you will be able to select file names for the output, and give the names of the variants you want included in the analysis. Following epitope scoring you will also be given the option to use R to generate an output PDF with the key findings of the pipeline.

## Pipeline Structure
Below is a flowchart outlining the rough pipeline structure.
![Image](pipelineflowchart.jpg)

## Things to note

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)

