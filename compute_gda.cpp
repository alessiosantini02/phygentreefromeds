#include <iostream>
#include <fstream>
#include <vector>

#include "gda_element.hpp"

using namespace std;

int binary_search_eds_index(int string_index_da, const vector<int>& eds_sizes);

int main(int argc, char *argv[]){
    //apertura file con le dimensioni scritte
    string sizes_filename=argv[1];
    FILE *sizes = fopen((sizes_filename+".bitvector.bin").c_str(), "rb");
    
    //scorrimento del file e salvataggio degli interi in un array
    vector<int> eds_sizes;
    int temp;
    while (fread(&temp, sizeof(int), 1, sizes) > 0)
    {
        if (eds_sizes.size()==0)
        {
            eds_sizes.push_back(temp);
        }else{
            int inserimento = temp+eds_sizes[eds_sizes.size()-1];
            eds_sizes.push_back(inserimento);
        }
    }
    fclose(sizes);
    
    //creazione array GDA
    vector<gda_element> gda;
    
    //apertura del file del DA
    string da_filename = argv[1];
    FILE *da = fopen((da_filename+".4.da").c_str(), "rb");
    
    //scorrimento del file del da leggendo gli elementi e calcolando direttamente il gda
    int buffer;
    while (fread(&buffer, sizeof(int), 1, da)>0)
    {
        gda_element temp;
        
        //calcolo della coppia GDA per ogni elemento di DA
        temp.eds_index = binary_search_eds_index(buffer, eds_sizes);
        if (temp.eds_index == 0)
        {
            temp.dollar_index = buffer;
        }else
        {
            temp.dollar_index = buffer-(eds_sizes[temp.eds_index-1]);
        }
        
        //inserimento nell'array gda
        gda.push_back(temp);
    }
    fclose(da);
    
    //salvataggio del gda in un nuovo file
    string gda_filename = argv[1];
    FILE *gda_file = fopen((gda_filename+"_gda.bin").c_str(), "wb");
    
    size_t written = fwrite(gda.data(), sizeof(gda_element), gda.size(), gda_file);
    if (written != gda.size()) {
        perror("Errore nella scrittura del file GDA");
    }
    fclose(gda_file);

    //salvataggio del numero di eds in un file ausiliario
    FILE *eds_number = fopen("eds_number.aux", "wb");    
    int gda_size = eds_sizes.size();
    fwrite(&gda_size, sizeof(int), 1, eds_number);
    fclose(eds_number);
    return 0;
}

int binary_search_eds_index(int string_index_da, const vector<int>& eds_sizes) {
    int left = 0, right = eds_sizes.size() - 1;
    
    while (left <= right)
    {
        int middle = (left+right)/2;

        if (string_index_da < eds_sizes[middle])
        {
            if (middle == 0 || string_index_da >= eds_sizes[middle-1])
            {
                //se middle è la dimensione della prima stringa
                //o se l'indice è compreso fra dimensione della middle-esima e della middle-1-esima
                return middle;
            }else
            {
                right=middle;
            }            
        }else
        {
            if (middle == eds_sizes.size()-1 || string_index_da <= eds_sizes[middle+1])
            {
                //se l'indice è compreso fra la middle-esima e la middle+1-esima dimensione
                return middle+1; 
            }else
            {
                left = middle;
            }
        }
    }
    
    return -1;  //non trovato
}