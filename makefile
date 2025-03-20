all: fasta_bv_concat compute_gda compute_distance_bwt_run compute_distance_gda_run compute_distance

fasta_bv_concat: fasta_bv_concat.cpp
		g++ -std=c++11 -fdiagnostics-color=always -g fasta_bv_concat.cpp -o fasta_bv_concat

compute_gda: compute_gda.cpp
		g++ -std=c++11 -fdiagnostics-color=always -g compute_gda.cpp -o compute_gda


#DOLLAR_COUNT ?= 0 #se =0 (valore di default) la prima parte della bwt con tutti i caratteri che precedono dollari non viene contata per la distanza
#DOLLAR_INTO_PREVIOUS_RUN ?= 1 # se =1 (valore di default) i dollari verranno accorpati al run che li precede nel calcolo della distanza sui blocchi della bwt, fungeranno da caratteri jolly quindi
compute_distance_bwt_run: compute_distance_bwt_run.cpp
		g++ -std=c++11 -fdiagnostics-color=always -g compute_distance_bwt_run.cpp -o compute_distance_bwt_run -DDOLLAR_COUNT=$(DOLLAR_COUNT)

compute_distance_gda_run: compute_distance_gda_run.cpp
		g++ -std=c++11 -fdiagnostics-color=always -g compute_distance_gda_run.cpp -o compute_distance_gda_run -DDOLLAR_COUNT=$(DOLLAR_COUNT)

compute_distance: compute_distance.cpp
		g++ -std=c++11 -fdiagnostics-color=always -g compute_distance.cpp -o compute_distance -DDOLLAR_COUNT=$(DOLLAR_COUNT) -DDOLLAR_INTO_PREVIOUS_RUN=$(DOLLAR_INTO_PREVIOUS_RUN)


clean:
		rm -f fasta_bv_concat compute_gda compute_distance_bwt_run compute_distance_gda_run compute_distance
