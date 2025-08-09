# Input Files
A list of file names, with this syntax, without extension:
    filename1 filename2 [...]

# Output files

Gsufsort output:
- Document array (.da) associated to the BWT of the collection of EDSs;
- file containing the BWT of the collection of EDSs (.bwt);

EDS-BWT output:
- a bitvector (.bitvector) for each EDS;
- a FASTA (.fasta) for each EDS;
- an information file about empty symbols (.info) for each EDS file;

phygentree output: 
- General Document Array (.bin) associated to the BWT of the collection of EDSs;
- General Bitvector associated to the concatenation of the collection of EDSs;
- matrix containing distances between EDSs (.phy);
- phylogenetic tree of the EDSs (.txt)

# Compile

- download code or clone this repository and move in the created directory:
```
git clone https://github.com/alessiosantini02/phygentree.git
```

- install EDS-BWT in phygentree directory:
```
cd phygentree/
git clone --recursive https://github.com/giovannarosone/EDS-BWT.git 
```

- install gsufsort in EDS-BWT directory and compile:
```
cd EDS-BWT
git clone https://github.com/felipelouza/gsufsort.git
cd gsufsort
make TERMINATOR=0 DNA=1
```

- compile EDS-BWT:
```
cd ..
make
```

- compile phygentree:
```
cd ..
make
```

- install rapidNJ (not necessarily in phygentree directory)
```
git clone https://github.com/somme89/rapidNJ.git
cd rapidNJ-master
make
```
and move the executable file from rapidNJ-master/bin to phygentree/

# Run

move in phygentree folder and type the following command followed by the paths of the EDS files (without extension)

```
./phygentree.sh edsfile1 edsfile2 [...]
```
If EDS files are not well formed or you are not sure about that, use the EDS-BWT function stringCheck like this:
```
cd EDS-BWT/
./stringCheck input.eds output
```
It will correct the EDS file and will save it as output.eds

# Improvements

Here's a list of things that could be improved by anyone who will use this tool:

- Files are read one symbol at a time, it could be optimized by reading blocks;
- names of species in distance matrixes are specified in code, defining an array named variants_name, it could be passed as argument;
- in functions compute_distance and compute_distance_sap should be added a case for N characters, that is considered as a normal symbol of the alphabet, so if two sequences present an N in the same position they are considered equal in that position, and this is wrong in general.
