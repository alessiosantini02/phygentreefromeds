#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

void fasta_concat(const string& concatfasta, string edsfilenames[], int filenumber);
void bv_concat(const string& concatbv, string edsfilenames[], int filenumber);


int main(int argc, char *argv[]){
    //salvataggio dei nomi passati come parametro in un array
    int filenumber = argc-2;
    string edsfilenames[filenumber];

    for (int i = 0; i < filenumber; i++)
    {
        edsfilenames[i]=argv[i+1];
        //edsfilenames[i]+= ".fasta";
    }

    //creazione fasta
    string concatfasta=argv[filenumber+1];
    concatfasta+=".fasta";
    fasta_concat(concatfasta, edsfilenames, filenumber);

    //creazione bitvector
    string concatbv=argv[filenumber+1];
    concatbv+=".bitvector";
    bv_concat(concatbv, edsfilenames, filenumber);
    
    return 0;
}



void fasta_concat(const string& concatfasta, string edsfilenames[], int filenumber){
    FILE *destinationfile = fopen(concatfasta.c_str(), "a");
    

    //ciclo che appende il contenuto di ogni elemento di edsfilenames nel file creato
    for (int i = 0; i < filenumber; i++)
    {
        char buffer[1024];
        size_t bytes;
        FILE *src=fopen((edsfilenames[i]+".fasta").c_str(), "r");
        while ((bytes = fread(buffer, 1, sizeof(buffer), src)) > 0) {
            fwrite(buffer, 1, bytes, destinationfile);
        }
        fclose(src);
    }
}

void bv_concat(const string& concatbv, string edsfilenames[], int filenumber){
    FILE *destinationfile = fopen(concatbv.c_str(), "ab");

    //ciclo che appende il contenuto di ogni bitvector nel file creato e memorizza le dimensioni delle eds in un altro file
    vector<int> eds_sizes(filenumber);
    for (int i = 0; i < filenumber; i++)
    {
        char buffer[1024];
        size_t bytes;
        FILE *src=fopen((edsfilenames[i]+".bitvector").c_str(), "rb");
        while ((bytes = fread(buffer, 1, sizeof(buffer), src)) > 0) {
            fwrite(buffer, 1, bytes, destinationfile);
            if (eds_sizes[i] == 0) //la dimensione è nel primo byte del bitvector, quindi ce lo scrivo e poi dopo non verrà modificato
            {
                eds_sizes[i]=buffer[0];
            }
            
        }
        fclose(src);
    }

    //salvataggio dei valori di eds_size in un file binario
    FILE *sizes = fopen((concatbv+".bin").c_str(), "wb");
    fwrite(eds_sizes.data(), sizeof(int), eds_sizes.size(), sizes);
    fclose(sizes);
}