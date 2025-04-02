#include <iostream>
#include <fstream>
#include <vector>
#include <limits>

#include "gda_element.hpp"

using namespace std;
void print_distance_matrix(vector<vector<double>> distance_matrix);
void print_as_list(vector<vector<double>>distance_matrix);
void store_distance_matrix_on_phylip_file(vector<vector<double>> distance_matrix, string name_distance_matrix_file);
void normalizza_matrice_distanze(vector<vector<double>> &distance_matrix);

int main(int argc, char* argv[]){
    //lettura del numero di EDS da un file ausiliario
    FILE *eds_number_file = fopen("eds_number.aux", "rb");
    int eds_number; //lettura da un file ausiliario in cui ho scritto la dimensione dell'array in compute_gda da eliminare dopo questa funzione
    fread(&eds_number, sizeof(int), 1, eds_number_file);
    fclose(eds_number_file);

    //inzializzazione matrice distanze
    vector<vector<double>> distance_matrix_gda(eds_number, vector<double>(eds_number, 0.0));
    fill(distance_matrix_gda.begin(), distance_matrix_gda.end(), vector<double>(eds_number, 0.0));

    //apertura file PDA e caricamento in un array
    string filename = argv[1];
    FILE *gda_file = fopen((filename+"_gda.bin").c_str(), "rb");
    vector<gda_element> pda;
    int eds_index; int string_index;
    while (fread(&eds_index, sizeof(int), 1, gda_file)>0)
    {
        fread(&string_index, sizeof(int), 1, gda_file);
        gda_element temp;
        temp.eds_index=eds_index;
        temp.dollar_index=string_index;
        pda.push_back(temp);
    }

    /*string variants_name[31]={"19A", "19B", "20A", "20B", "20C", "20D", "20E", "20F", "20G", "20H (Beta)", "20I (Alpha)", "20J (Gamma)", "21A (Delta)", "21B (Kappa)", "21C (Epsilon)", "21D (Eta)",
     "21F (Iota)", "21G (Lambda)", "21H (Mu)", "21I (Delta)", "21J (Delta)", "21K (BA.1)",
      "21L (BA.2)", "22A (BA.4)", "22B (BA.5)", "22C (BA.2.12.1)", "22D (BA.2.75)", "22E (BQ.1)", "22F (XBB)", "23A (XBB.1.5)", "23B (XBB.1.16)"};
*/
    string variants_name[31]={"19A", "19B", "21A (Delta)", "21I (Delta)", "21J (Delta)"};

    //scorrimento del triangolo inferiore delle matrici per calcolare le distanze a due a due
    for (size_t i = 0; i < eds_number; i++)
    {
        int j;
        for (j = 0; j<i; j++)
        {
            //variabili per la distanza sul pda
            int current_run_counter=1; //contatore per la lunghezza del run attuale, inizializzato a 1 così alla prima iterazione viene aggiunto al contatore 1-1=0
            int current_string=0; //0 solo all'inizio, 1 se sta contando un run della stringa i, 2 se j
            int symbols_number=0;

            //apertura file SAP delle eds i e j
            string sapfilename="/home/alessio/Scrivania/fastaconcatenatiadueadue/";
            sapfilename=sapfilename+variants_name[j]+variants_name[i]+".fasta.sap";
            FILE *sapfile=fopen(sapfilename.c_str(), "r");
            vector<int> saparray;
            char buffer;
            while (fread(&buffer, sizeof(char), 1, sapfile)>0){
                int current_sap_value = buffer - '0';
                saparray.push_back(current_sap_value);
            }

            vector<int> pda_projected;
            for (size_t k = 0; k < pda.size(); k++)
            {
                if (pda[k].eds_index == i || pda[k].eds_index == j)
                {
                    pda_projected.push_back(pda[k].eds_index);
                }
            }

            for(int l=0; l<saparray.size(); l++)
            {
                if (saparray[l]==0 || saparray[l-1]==0)
                {//posso contare
                    //cout<<"valore sap =0, entro nell'if"<<endl;
                    if (pda_projected[l] == i)
                    {
                        //cout<<"pda[k]=i, entro nell'if"<<endl;
                        if (current_string != 1)
                        {
                            distance_matrix_gda[i][j] += current_run_counter-1;
                            //cout<<"salvo nella matrice"<<current_run_counter-1<<endl;
                            current_string = 1;
                            current_run_counter = 1;
                        }else
                        {
                            //cout<<"aggiungo uno al contatore del run"<<endl;
                            current_run_counter++;
                        }
                        symbols_number++;
                    }else if (pda_projected[l] == j)
                    {
                        //cout<<"pda[k]=j, entro nell'if"<<endl;
                        if (current_string != 2)
                        {
                            distance_matrix_gda[i][j] += current_run_counter-1;
                            //cout<<"salvo nella matrice"<<current_run_counter-1<<endl;
                            current_string = 2;
                            current_run_counter = 1;
                        }else
                        {
                            current_run_counter++;
                            //cout<<"aggiungo uno al contatore del run"<<endl;
                        }
                        symbols_number++;
                    }
                }else
                {
                    //cout<<"azzero i contatore"<<endl;
                    current_run_counter=1; //contatore per la lunghezza del run attuale, inizializzato a 1 così alla prima iterazione viene aggiunto al contatore 1-1=0
                    current_string=0; //0 solo all'inizio, 1 se sta contando un run della stringa i, 2 se j
                }
            }

            rewind(gda_file);

            //conteggio ultimo run
            distance_matrix_gda[i][j] += current_run_counter-1;

            distance_matrix_gda[j][i] = distance_matrix_gda[i][j]; //per completare la matrice che è simmetrica
        }
    }

    fclose(gda_file);
    
    normalizza_matrice_distanze(distance_matrix_gda);

    cout<<"Matrice PDA-SAP"<<endl;
    print_distance_matrix(distance_matrix_gda);
    print_as_list(distance_matrix_gda);

    store_distance_matrix_on_phylip_file(distance_matrix_gda, "risultati/distance_matrix_pda_sap");

    return 0;
}

void print_distance_matrix(vector<vector<double>> distance_matrix){
    for (size_t i = 0; i < distance_matrix.size(); i++)
    {
        for (size_t j = 0; j < distance_matrix.size(); j++)
        {
            cout << distance_matrix[i][j]<< " ";
        }
        cout << endl;
    }
}

void print_as_list(vector<vector<double>>distance_matrix){
    /*string variants_name[31]={"19A", "19B", "20A", "20B", "20C", "20D", "20E", "20F", "20G", "20H (Beta)", "20I (Alpha)", "20J (Gamma)", "21A (Delta)", "21B (Kappa)", "21C (Epsilon)", "21D (Eta)",
     "21F (Iota)", "21G (Lambda)", "21H (Mu)", "21I (Delta)", "21J (Delta)", "21K (BA.1)",
      "21L (BA.2)", "22A (BA.4)", "22B (BA.5)", "22C (BA.2.12.1)", "22D (BA.2.75)", "22E (BQ.1)", "22F (XBB)", "23A (XBB.1.5)", "23B (XBB.1.16)"};
*/
    string variants_name[31]={"19A", "19B", "21A (Delta)", "21I (Delta)", "21J (Delta)"};
    
    for (size_t i = 0; i < distance_matrix.size(); i++)
    {
        for (size_t j = 0; j < i; j++)
        {
            string nodo1=variants_name[i];
            string nodo2=variants_name[j];
            double distance=distance_matrix[i][j];
            cout<<"('"<<nodo1<<"','"<<nodo2<<"',"<<distance<<"),"<<endl;
        }
    }
    
}

void store_distance_matrix_on_phylip_file(vector<vector<double>> distance_matrix, string name_distance_matrix_file){
    ofstream phylip_file((name_distance_matrix_file+".phy").c_str());
    phylip_file << distance_matrix.size() << "\n"; //numero di elementi

    string variants_name[31]={"19A", "19B", "20A", "20B", "20C", "20D", "20E", "20F", "20G", "20H", "20I", "20J", "21A", "21B", "21C", "21D", "21F", "21G", "21H", "21I", "21J", "21K", "21L", "22A", "22B", "22C", "22D", "22E", "22F", "23A", "23B"};

    for (size_t i = 0; i < distance_matrix.size(); i++) {
        phylip_file <<  variants_name[i]/*"E" << i */<< "  ";
        for (double d : distance_matrix[i]) {
            phylip_file << d << "  ";
        }
        phylip_file << "\n";
    }

    phylip_file.close();
}

void normalizza_matrice_distanze(vector<vector<double>> &distance_matrix){
    double max=0;
    
    for (size_t i = 0; i < distance_matrix.size(); i++)
    {
        for (size_t j = 0; j<i; j++)
        {
            if (distance_matrix[i][j]>max)
            {
                max=distance_matrix[i][j];
            }
        }
    }

    for (size_t i = 0; i < distance_matrix.size(); i++)
    {
        for (size_t j = 0; j<i; j++)
        {
            distance_matrix[i][j]=(distance_matrix[i][j]) / (max);
            distance_matrix[j][i]=distance_matrix[i][j];
        }
    }
}