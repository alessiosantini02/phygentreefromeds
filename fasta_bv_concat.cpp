#include <iostream>
#include <fstream>

using namespace std;

void fasta_concat(const string& concatfasta, string edsfilenames[], int filenumber);
void bv_concat(const string& concatbv, string edsfilenames[], int filenumber);


int main(int argc, char *argv[]){
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

    //ciclo che appende il contenuto di ogni elemento di edsfilenames nel file creato
    for (int i = 0; i < filenumber; i++)
    {
        char buffer[1024];
        size_t bytes;
        FILE *src=fopen((edsfilenames[i]+".bitvector").c_str(), "rb");
        while ((bytes = fread(buffer, 1, sizeof(buffer), src)) > 0) {
            fwrite(buffer, 1, bytes, destinationfile);
        }
        fclose(src);
    }
}