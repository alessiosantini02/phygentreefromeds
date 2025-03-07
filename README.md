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
