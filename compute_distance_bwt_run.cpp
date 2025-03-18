#include <iostream>
#include <fstream>
#include <vector>

#include "gda_element.hpp"

using namespace std;

int main(int argc, char* argv[]){
    //lettura del numero di EDS da un file ausiliario
    FILE *eds_number_file = fopen("eds_number.aux", "rb");
    int eds_number; //lettura da un file ausiliario in cui ho scritto la dimensione dell'array in compute_gda da eliminare dopo questa funzione
    fread(&eds_number, sizeof(int), 1, eds_number_file);
    fclose(eds_number_file);

    //inzializzazione matrice distanze
    vector<vector<double>> distance_matrix(eds_number, vector<double>(eds_number));
    
    //apertura file GDA e BWT
    string filename = argv[1];
    FILE *gda_file = fopen((filename+"_gda.bin").c_str(), "rb");
    FILE *bwt_file = fopen((filename+".bwt").c_str(), "rb");

    //ciclo for che scorre il triangolo inferiore della matrice
    for (size_t i = 0; i < eds_number; i++)
    {
        for (size_t j = 0; j < eds_number; j++)
        {
            if (i==j)
            {
                distance_matrix[i][j] = 0; //per convenzione, per com'è definita la distanza non è esattamente 0
            }else if (i>j)
            {
                char previous_symbol='b'; //per vedere se il carattere corrente appartiene allo stesso run di prima, inizializzo con un carattere non nell'alfabeto
                int gda_difference=0; //contatore cumulativo su cui si fa +1 o -1 in base al valore del gda
                int symbols_counter=0; //per la normalizzazione
                int distance_counter=0; //contatore della distanza fra la EDS i e la j

                char bwt_symbol_buffer; 
                int gda_buffer;
                while (fread(&bwt_symbol_buffer, 1, sizeof(bwt_symbol_buffer), bwt_file) > 0)
                {
                    fread(&gda_buffer, sizeof(int), 1, gda_file);

                    if (bwt_symbol_buffer == previous_symbol)
                    {
                        //guardo il valore del gda, se è i o j faccio +1 su un contatore cumulativo
                        if (gda_buffer == i)
                        {
                            gda_difference++;
                            symbols_counter++;
                        }else if (gda_buffer == j)
                        {
                            gda_difference--;
                            symbols_counter++;
                        }
                    }else
                    {
                        distance_matrix[i][j]+=abs(gda_difference);
                        gda_difference=0;
                        previous_symbol=bwt_symbol_buffer;
                    }
                }

                //conteggio ultimo run
                //guardo il valore del gda, se è i o j faccio +1 su un contatore cumulativo
                if (gda_buffer == i)
                {
                    gda_difference++;
                    symbols_counter++;
                }else if (gda_buffer == j)
                {
                    gda_difference--;
                    symbols_counter++;
                }
                distance_matrix[i][j]+=abs(gda_difference);
                
                //inserimento distanza fra EDS i e EDS j nella matrice
                distance_matrix[i][j]/=symbols_counter;
                distance_matrix[j][i]=distance_matrix[i][j];
                
                rewind(bwt_file);
                rewind(gda_file);                
            }
            
            
        }
        
    }

    //salvataggio della matrice in phylip per il momento, poi aggiungere opzione per csv
    string name_distance_matrix_file = argv[2];
    ofstream phylip_file((name_distance_matrix_file+".phy").c_str());
    phylip_file << distance_matrix.size() << "\n"; //numero di elementi

    for (size_t i = 0; i < distance_matrix.size(); i++) {
        phylip_file << "E" << i << "  ";
        for (double d : distance_matrix[i]) {
            phylip_file << d << "  ";
        }
        phylip_file << "\n";
    }

    phylip_file.close();
    
    return 0;
}