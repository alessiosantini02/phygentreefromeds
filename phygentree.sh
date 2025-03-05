#!/bin/bash

EDSBWTPATH="EDS-BWT"
GSUFPATH="EDS-BWT/gsufsort"
OUTPUT1="conc_file"
OUTPUT2="distance_matrix"
OUTPUT3="phygentree"

#salvataggio di tutti i filename passati come parametro in un array
file_array=("$@")

#funzione che fa fasta e i bv di ogni eds
for file in "${file_array[@]}"; do
    echo "Processo il file: $file"
    $EDSBWTPATH/eds_to_fasta $file.eds $file
done

#funzione che concatena i fasta e i bitvector
./fasta_bv_concat "${file_array[@]}" $OUTPUT1

#bwt del fasta concatenato
$GSUFPATH/gsufsort $OUTPUT1.fasta --da --bwt --output $OUTPUT1
rm $OUTPUT1.fasta 

#funzione che calcola il gda leggendo il da
./compute_gda $OUTPUT1
rm $OUTPUT1.bitvector.bin

#funzione che calcola la distanza fra le eds sulla bwt
./compute_distance_bwt $OUTPUT1 $OUTPUT2
rm eds_number.aux

#funzione che costruisce l'albero filogenetico in base alla matrice distanze (nome funzione non definitivo)
#./neighbour_join $OUTPUT2 $OUTPUT3
./rapidnj $OUTPUT2.phy > $OUTPUT3.txt