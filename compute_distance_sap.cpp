#include <iostream>
#include <fstream>
#include <vector>
#include <limits>

#include "gda_element.hpp"

using namespace std;
void print_distance_matrix(vector<vector<double>> distance_matrix);
void print_as_list(vector<vector<double>>distance_matrix_gda);

int main(int argc, char* argv[]){
    //lettura del numero di EDS da un file ausiliario
    FILE *eds_number_file = fopen("eds_number.aux", "rb");
    int eds_number; //lettura da un file ausiliario in cui ho scritto la dimensione dell'array in compute_gda da eliminare dopo questa funzione
    fread(&eds_number, sizeof(int), 1, eds_number_file);
    fclose(eds_number_file);

    //inzializzazione matrici distanze
    vector<vector<double>> distance_matrix_gda(eds_number, vector<double>(eds_number));

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

    string variants_name[31]={"19A", "19B", "20A", "20B", "20C", "20D", "20E", "20F", "20G", "20H (Beta)", "20I (Alpha)", "20J (Gamma)", "21A (Delta)", "21B (Kappa)", "21C (Epsilon)", "21D (Eta)",
     "21F (Iota)", "21G (Lambda)", "21H (Mu)", "21I (Delta)", "21J (Delta)", "21K (BA.1)",
      "21L (BA.2)", "22A (BA.4)", "22B (BA.5)", "22C (BA.2.12.1)", "22D (BA.2.75)", "22E (BQ.1)", "22F (XBB)", "23A (XBB.1.5)", "23B (XBB.1.16)"};
    
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

    //normalizzazione
    double max=0;
    double min=std::numeric_limits<double>::max();
    
    for (size_t i = 0; i < distance_matrix_gda.size(); i++)
    {
        for (size_t j = 0; j<i; j++)
        {
            if (distance_matrix_gda[i][j]>max)
            {
                max=distance_matrix_gda[i][j];
            }else if (distance_matrix_gda[i][j]<min)
            {
                min=distance_matrix_gda[i][j];
            }   
        }
    }
    cout<<"massimo: "<<max<<" minimo: "<<min<<endl;

    for (size_t i = 0; i < distance_matrix_gda.size(); i++)
    {
        for (size_t j = 0; j<i; j++)
        {
            distance_matrix_gda[i][j]= ((distance_matrix_gda[i][j]-min) / (max -min));
            distance_matrix_gda[j][i]=distance_matrix_gda[i][j];
        }
    }

    print_distance_matrix(distance_matrix_gda);
    print_as_list(distance_matrix_gda);

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

void print_as_list(vector<vector<double>>distance_matrix_gda){
    string variants_name[31]={"19A", "19B", "20A", "20B", "20C", "20D", "20E", "20F", "20G", "20H (Beta)", "20I (Alpha)", "20J (Gamma)", "21A (Delta)", "21B (Kappa)", "21C (Epsilon)", "21D (Eta)",
     "21F (Iota)", "21G (Lambda)", "21H (Mu)", "21I (Delta)", "21J (Delta)", "21K (BA.1)",
      "21L (BA.2)", "22A (BA.4)", "22B (BA.5)", "22C (BA.2.12.1)", "22D (BA.2.75)", "22E (BQ.1)", "22F (XBB)", "23A (XBB.1.5)", "23B (XBB.1.16)"};

    for (size_t i = 0; i < distance_matrix_gda.size(); i++)
    {
        for (size_t j = 0; j < i; j++)
        {
            string nodo1=variants_name[i];
            string nodo2=variants_name[j];
            double distance=distance_matrix_gda[i][j];
            cout<<"('"<<nodo1<<"','"<<nodo2<<"',"<<distance<<"),"<<endl;
        }
    }
    
}