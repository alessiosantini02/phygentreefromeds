#include <iostream>
#include <fstream>
#include <vector>
#include <limits>

#include "gda_element.hpp"

using namespace std;

void calcolo_distanza_gda(int i, int j, vector<vector<double>> &distance_matrix_gda, vector<int> gda, int eds_sizes);
void calcolo_distanza_bwt(int i, int j, vector<vector<double>> &distance_matrix_bwt, vector<int> gda, vector<char> bwt, int eds_sizes);
void normalizza_matrice_distanze(vector<vector<double>> &distance_matrix);
void print_distance_matrix(vector<vector<double>> distance_matrix);
void store_distance_matrix_on_phylip_file(vector<vector<double>> distance_matrix, string name_distance_matrix_file);
void print_as_list (vector<vector<double>> distance_matrix);

int main(int argc, char* argv[]){

    //lettura del numero di EDS da un file ausiliario
    FILE *eds_number_file = fopen("eds_number.aux", "rb");
    int eds_number; //lettura da un file ausiliario in cui ho scritto la dimensione dell'array in compute_gda da eliminare dopo questa funzione
    fread(&eds_number, sizeof(int), 1, eds_number_file);
    fclose(eds_number_file);

#if (DOLLAR_COUNT==0)
    //apertura file con le dimensioni scritte
    string eds_dim_filename=argv[1];
    FILE *sizes = fopen((eds_dim_filename+".bitvector.bin").c_str(), "rb");
    
    //scorrimento del file e salvataggio degli interi in un array
    int eds_sizes;
    int temp;
    while (fread(&temp, sizeof(int), 1, sizes) > 0)
    {
        eds_sizes+=temp;
    }
    fclose(sizes);
#endif

    //inzializzazione matrici distanze
    vector<vector<double>> distance_matrix_gda(eds_number, vector<double>(eds_number, 0.0));
    vector<vector<double>> distance_matrix_bwt(eds_number, vector<double>(eds_number, 0.0));
    fill(distance_matrix_gda.begin(), distance_matrix_gda.end(), vector<double>(eds_number, 0.0));
    fill(distance_matrix_bwt.begin(), distance_matrix_bwt.end(), vector<double>(eds_number, 0.0));

    //apertura file GDA e BWT e salvataggio in array
    string filename = argv[1];
    FILE *gda_file = fopen((filename+"_gda.bin").c_str(), "rb");
    FILE *bwt_file = fopen((filename+".bwt").c_str(), "rb");
    vector<int> gda;
    vector<char> bwt;
    char bwt_symbol_buffer; 
    int gda_buffer;
    while (fread(&bwt_symbol_buffer, 1, sizeof(bwt_symbol_buffer), bwt_file) > 0){
        bwt.push_back(bwt_symbol_buffer);
        
        fread(&gda_buffer, sizeof(int), 1, gda_file);
        gda.push_back(gda_buffer);

        //valore da buttare
        int bin;
        fread(&bin, sizeof(int), 1, gda_file);
    }
    fclose(gda_file);
    fclose(bwt_file);


    //scorrimento del triangolo inferiore delle matrici per calcolare le distanze a due a due
    for (size_t i = 0; i < eds_number; i++)
    {
        int j;
        for (j = 0; j<=i; j++)
        {
            if (i==j)
            {
                distance_matrix_gda[i][j] = 0; 
                distance_matrix_bwt[i][j] = 0;
            }else if (i>j)
            {
                calcolo_distanza_gda(i,j, distance_matrix_gda, gda, eds_sizes);
                calcolo_distanza_bwt(i,j,distance_matrix_bwt, gda, bwt, eds_sizes);
            }
        }
    }

    normalizza_matrice_distanze(distance_matrix_bwt);
    normalizza_matrice_distanze(distance_matrix_gda);
    
    cout<<"inizio matrice bwt"<<endl;
    print_distance_matrix(distance_matrix_bwt);
    print_as_list(distance_matrix_bwt);
    cout<<endl;

    cout<<"inizio matrice gda"<<endl;
    print_distance_matrix(distance_matrix_gda);
    print_as_list(distance_matrix_gda);
    cout<<endl;

    //salvataggio della matrice in phylip per il momento, poi aggiungere opzione per csv
    string name_distance_matrix_file = argv[2];
    store_distance_matrix_on_phylip_file(distance_matrix_bwt, name_distance_matrix_file+"_bwt");
    store_distance_matrix_on_phylip_file(distance_matrix_gda, name_distance_matrix_file+"_gda");

    return 0;
}

void calcolo_distanza_gda(int i, int j, vector<vector<double>> &distance_matrix_gda, vector<int> gda, int eds_sizes){
    //variabili per la distanza sul gda
    int current_run_counter=1; //contatore per la lunghezza del run attuale, inizializzato a 1 così alla prima iterazione viene aggiunto al contatore 1-1=0
    int current_string=0; //0 solo all'inizio, 1 se sta contando un run della stringa i, 2 se j
    int symbols_number=0;

    
    #if (DOLLAR_COUNT==0)
        int k=0;
    #endif
    
    int counter=0;
    for(counter=0; counter<gda.size(); counter++)
    {
    #if (DOLLAR_COUNT==0)
        if (k>eds_sizes)
        {
    #endif

        //DISTANZA SUI RUN DEL GDA
        //controllo il valore del GDA per controllare se fa parte del run corrente o inizia un nuovo run
        if (gda[counter] == i)
        {
            if (current_string != 1)
            {
                distance_matrix_gda[i][j] += current_run_counter-1;
                current_string = 1;
                current_run_counter = 1;
            }else
            {
                current_run_counter++;
            }
            symbols_number++;
        }else if (gda[counter] == j)
        {
            if (current_string != 2)
            {
                distance_matrix_gda[i][j] += current_run_counter-1;
                current_string = 2;
                current_run_counter = 1;
            }else
            {
                current_run_counter++;
            }
            symbols_number++;
        }

#if (DOLLAR_COUNT==0)
    }
    k++;
#endif
    }

    //conteggio dell'ultimo run del gda
    if (gda[counter] == i || gda[counter] == j)
    {
        distance_matrix_gda[i][j] += current_run_counter-1;
    }

    //distance_matrix_gda[i][j] = distance_matrix_gda[i][j] / symbols_number; //normalizzazione per la somma
    distance_matrix_gda[j][i] = distance_matrix_gda[i][j]; //per completare la matrice che è simmetrica
}

void calcolo_distanza_bwt(int i, int j, vector<vector<double>> &distance_matrix_bwt, vector<int> gda, vector<char> bwt, int eds_sizes){
    //variabili per la distanza sulla bwt
    char previous_symbol='b'; //per vedere se il carattere corrente appartiene allo stesso run di prima, inizializzo con un carattere non nell'alfabeto
    int gda_difference=0; //contatore cumulativo su cui si fa +1 o -1 in base al valore del gda
    int symbols_counter=0; //per la normalizzazione
    int distance_counter=0; //contatore della distanza fra la EDS i e la j

    #if (DOLLAR_COUNT==0)
        int k=0;
    #endif
    
    int counter=0;
    for(counter=0; counter<bwt.size(); counter++)
    {
    #if (DOLLAR_COUNT==0)
        if (k>eds_sizes)
        {
    #endif
        //DISTANZA SUI RUN DELLA BWT
        if((gda[counter]==i || gda[counter]==j)){
            #if DOLLAR_INTO_PREVIOUS_RUN==1
            if (bwt[counter] == previous_symbol || (previous_symbol != 'b' && bwt[counter]=='#'))
            #else
            if (bwt[counter] == previous_symbol)
            #endif
            {                 
                //guardo il valore del gda, se è i o j faccio +1 su un contatore cumulativo
                if (gda[counter] == i)
                {
                    gda_difference++;
                    symbols_counter++;
                }else if (gda[counter] == j)
                {
                    gda_difference--;
                    symbols_counter++;
                }
            }else
            {
                distance_matrix_bwt[i][j]+=abs(gda_difference);
                gda_difference=0;
                previous_symbol=bwt[counter];
            }
        }
    #if (DOLLAR_COUNT==0)
    }
    k++;
    #endif
    }

    //conteggio ultimo run
    //guardo il valore del gda, se è i o j faccio +1 su un contatore cumulativo
    if (gda[counter] == i)
    {
        gda_difference++;
        symbols_counter++;
    }else if (gda[counter] == j)
    {
        gda_difference--;
        symbols_counter++;
    }
    distance_matrix_bwt[i][j]+=abs(gda_difference);
    
    //inserimento distanza fra EDS i e EDS j nella matrice
    //distance_matrix_bwt[i][j]/=symbols_counter;
    distance_matrix_bwt[j][i]=distance_matrix_bwt[i][j];
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

void print_as_list(vector<vector<double>> distance_matrix){
    string variants_name[31]={"19A", "19B", "20A", "20B", "20C", "20D", "20E", "20F", "20G", "20H (Beta)", "20I (Alpha)", "20J (Gamma)", "21A (Delta)", "21B (Kappa)", "21C (Epsilon)", "21D (Eta)",
     "21F (Iota)", "21G (Lambda)", "21H (Mu)", "21I (Delta)", "21J (Delta)", "21K (BA.1)",
      "21L (BA.2)", "22A (BA.4)", "22B (BA.5)", "22C (BA.2.12.1)", "22D (BA.2.75)", "22E (BQ.1)", "22F (XBB)", "23A (XBB.1.5)", "23B (XBB.1.16)"};
    
    //string variants_name[31]={"19A", "19B", "21A (Delta)", "21I (Delta)", "21J (Delta)"};
    
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

void print_distance_matrix(vector<vector<double>> distance_matrix){
    for (size_t i = 0; i < distance_matrix.size(); i++)
    {
        cout<<"[";
        for (size_t j = 0; j < distance_matrix.size(); j++)
        {
            cout << distance_matrix[i][j]<< ", ";
        }
        cout <<"]," <<endl;
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