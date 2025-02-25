# Input Files
A list of file names, with this syntax, without extension:
    filename1 filename2 ...

# Output files

Gsufsort output:
- Document array (.da)
- .bwt

EDS-BWT output:
- Bitvector (.bitvector)
- (FASTA (.fasta))
- .info

phygentree output: 
- 


# Compile
[mettere istruzioni di compilazione di gsufsort e di eds-bwt, e poi il comando make in questa cartella]
- install EDS-BWT
```git clone --recursive https://github.com/giovannarosone/EDS-BWT.git 
cd EDS-BWT```

- install and compile gsufsort:
```git clone https://github.com/felipelouza/gsufsort.git
cd gsufsort
make TERMINATOR=0 DNA=1
cd ..```

- compile EDS-BWT
```make
# Run
move in the folder```

```./phygentree.sh edsfile1 edsfile2 [...]```
