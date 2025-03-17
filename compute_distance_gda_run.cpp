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
    
    //apertura file GDA
    string gdafilename = argv[1];
    FILE *gda_file = fopen((gdafilename+"_gda.bin").c_str(), "rb");

    //scorrimento del GDA per calcolare le distanze a due a due fra le EDS
    for (size_t i = 0; i < eds_number; i++)
    {
        int j;
        for (j = 0; j < eds_number; j++)
        {
            int symbols_number=0; //il numero di simboli delle due eds di cui si calcola la distanza, si usa come denominatore per la normalizzazione

            if (i==j)
            {
                distance_matrix[i][j] = 0; //per convenzione, per com'è definita la distanza non è esattamente 0
            }else if(i>j) //il calcolo della distanza si fa solo sul triangolo inferiore della matrice, per poi copiare il valore nella posizione corrispondente nel triangolo superiore
            {
                int current_run_counter=1; //contatore per la lunghezza del run attuale, inizializzato a 1 così alla prima iterazione viene aggiunto al contatore 1-1=0
                int current_string=0; //0 solo all'inizio, 1 se sta contando un run della stringa i, 2 se j

                int eds_index_temp;
                while (fread(&eds_index_temp, sizeof(int), 1, gda_file) > 0)
                {
                    //lettura di un altro elemento per buttarlo via perchè il GDA è memorizzato come sequenza di interi in binario, non è strutturato a coppie
                    int bin;
                    fread(&bin, sizeof(int), 1, gda_file);

                    //controllo il valore del GDA per controllare se fa parte del run corrente o inizia un nuovo run
                    if (eds_index_temp == i)
                    {
                        if (current_string != 1)
                        {
                            distance_matrix[i][j] += current_run_counter-1;
                            current_string = 1;
                            current_run_counter = 1;
                        }else
                        {
                            current_run_counter++;
                        }
                        symbols_number++;
                    }else if (eds_index_temp == j)
                    {
                        if (current_string != 2)
                        {
                            distance_matrix[i][j] += current_run_counter-1;
                            current_string = 2;
                            current_run_counter = 1;
                        }else
                        {
                            current_run_counter++;
                        }
                        symbols_number++;
                    }
                }

                rewind(gda_file);

                if (eds_index_temp == i || eds_index_temp == j)
                {
                    distance_matrix[i][j] += current_run_counter-1;// conteggio dell'ultimo run del gda
                }

                distance_matrix[i][j] = distance_matrix[i][j] / symbols_number; //normalizzazione per la somma
                distance_matrix[j][i] = distance_matrix[i][j]; //per completare la matrice che è simmetrica
            }
        }
    }

    fclose(gda_file);

    for (size_t i = 0; i < distance_matrix.size(); i++)
    {
        for (size_t j = 0; j < distance_matrix.size(); j++)
        {
            cout << distance_matrix[i][j]<< " ";
        }
        cout << endl;
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