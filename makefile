all: fasta_bv_concat compute_gda compute_distance_bwt

fasta_bv_concat: fasta_bv_concat.cpp
		g++ -std=c++11 -fdiagnostics-color=always -g fasta_bv_concat.cpp -o fasta_bv_concat

compute_gda: compute_gda.cpp
		g++ -std=c++11 -fdiagnostics-color=always -g compute_gda.cpp -o compute_gda

compute_distance_bwt: compute_distance_bwt.cpp
		g++ -std=c++11 -fdiagnostics-color=always -g compute_distance_bwt.cpp -o compute_distance_bwt

clean:
		rm -f fasta_bv_concat compute_gda compute_distance_bwt
