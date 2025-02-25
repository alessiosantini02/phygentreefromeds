/*

// Apertura file in lettura
FILE *gda_file = fopen("gda_output.bin", "rb");
vector<gda_element> gda_read(gda.size());  // Allocazione spazio

// Lettura del file
fread(gda_read.data(), sizeof(gda_element), gda.size(), gda_file);
fclose(gda_file);

// Stampa per verifica
for (const auto &elem : gda_read) {
    printf("EDS: %d, Dollaro: %d\n", elem.eds_index, elem.dollar_index);
}
*/

/*
- aprire la bwt e caricarla in un array di char
- aprire il gda e caricarlo in un array di gda_element
- fare un ciclo for annidato che:
    - carica in un altro array di char gli elementi appartenenti a due stringhe
    - conta le alternanze di colori
    - salva il valore in una matrice
- salva la matrice in un file (vedere formato, forse phylip che è lo standard in bionformatica)
*/