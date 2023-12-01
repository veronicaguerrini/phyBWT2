# phyBWT2

phyBWT2 is a new alignment-, assembly-, and reference-free method that builds a partition tree without relying on the pairwise comparison of sequences, thus avoiding to use a distance matrix to infer phylogeny.

It applies the properties of the Extended Burrows-Wheeler Transform (EBWT) to the idea of decomposition for phylogenetic inference. 
In particular, it hinges the combinatorial properties of the *EBWT positional clustering* framework recently introduced , overcoming the limitations of employing *k*-mers with a priori fixed size. 
Finally, phyBWT2 infers the tree structure by comparing all the sequences simultaneously, instead of performing their pairwise comparisons.

Let *S={S<sub>1</sub>,...,S<sub>n</sub>}* be the input collection of sequences, where each *S<sub>i</sub>* is a multiset of strings representing an organism (e.g. sequencing reads, contigs, genome). The tool phyBWT takes as input the following data structures:
- the extended Burrows–Wheeler transform (eBWT), or multi-string BWT, of collection *S*;
- the longest common prefix array (LCP) of collection *S*;
- the color document array (CDA) of collection *S*, obtained from the document array (DA) of collection *S*.

### Install

```sh
git clone https://github.com/veronicaguerrini/phyBWT2
cd phyBWT2
```

### Compile
phyBWT2 can handle datasets of different types (short reads, long reads, contigs, or entire genomes). 

For short reads, phyBWT2 must be compiled by using

```sh
make
```
while for sequences longer than 250 bp, set SHORT=0

```sh
make SHORT=0
```

### Preprocessing steps

The required data structures eBWT, LCP and CDA can be built independently from phyBWT2. 
This is a good feature that allows the user to choose the most appropriate tool according to the resources available and the dataset composition (short reads or longer sequences).

For instance, to build .ebwt, .lcp, and .da files from scratch from a single fasta file, one could use BCR [https://github.com/giovannarosone/BCR_LCP_GSA] for short reads, and for longer sequences, gsufsort [https://github.com/felipelouza/gsufsort]. Note that gsufsort tool returns the output files with slightly different filename extensions.

To install and compile both BCR and gsufsort for the preprocessing step, in addition to phyBWT2, one could run

```sh
./Install.sh
```
Note that by default the above script compiles phyBWT2 for short reads. To correctly compile phyBWT2 for longer sequences, the parameter SHORT inside the script must be set to 0.

To obtain the color document array CDA from the document array DA (file fastaName.da), one could use

```sh
./create_cda fastaName fileInfo
```
where fileInfo is a tab-separated file that stores per line the number of sequences of each *S<sub>i</sub>* of the collection *S*, according to the format:

*S<sub>i</sub>&emsp;num*

Alternatively, one could run the following script that builds up the three data structures eBWT, LCP and CDA starting from a collection of FASTA files (one for each multiset *S<sub>i</sub>* of the collection *S*) stored in the directory input_directory.

```sh
./Preprocessing.sh input_directory fastaName
```

The script outputs the three data structures eBWT, LCP and CDA (files fastaName.ebwt, fastaName.lcp and fastaName.cda), and a tab-separated file (fastaName.txt) that stores per line the number of sequences of each *S<sub>i</sub>* of the collection *S* (*i.e.* fileInfo). 
By default the above script computes the data structures eBWT, LCP and DA by using BCR (preferred choice for constructing the eBWT of short reads). To build them for longer sequences, please change the parameter SHORT inside the script and set it to 0.

### Run

Our alignment-, assembly-, and reference-free method can be run by using:

```sh
./phyBWT2 fastaName fileInfo output k_min tau
```
where fastaName is the base name of files fastaName.ebwt, fastaName.lcp and fastaName.cda, while fileInfo is a tab-separated file that stores per line the number of sequences of each *S<sub>i</sub>* of the collection *S*.

### Quick test

By setting SHORT=0 in Preprocessing.sh and compiling phyBWT2 with parameter SHORT=0, the input datastructures can be built by

```sh
./Preprocessing.sh ./example/hiv_sequences/ HIV-1_data
```
and phyBWT2 can be run by using

```sh
./phyBWT2 HIV-1_data.fasta HIV-1_data.txt Tree_HIV-1_data.new 16 0.6
```

## References

Guerrini V., Conte A., Grossi R., Liti G., Rosone G., Tattini L., phyBWT: Alignment-Free Phylogeny via eBWT Positional Clustering. In *Proceedings of the 22nd International Workshop on Algorithms in Bioinformatics* WABI 2022. LIPIcs, vol 242, pp 23:1--23:19. doi: [10.4230/LIPIcs.WABI.2022.23](https://doi.org/10.4230/LIPIcs.WABI.2022.23)

Guerrini V., Conte A., Grossi R., Liti G., Rosone G., Tattini L., phyBWT2: phylogeny reconstruction via eBWT positional clustering. (2023) Algorithms Mol Biol 18, 11. doi: [10.1186/s13015-023-00232-4](https://doi.org/10.1186/s13015-023-00232-4)

--

<small> Supported by PNRR project “THE—Tuscany Health Ecosystem” — Spoke 6 “Precision medicine & personalized healthcare”, funded by the European Commission under the NextGeneration EU programme.
