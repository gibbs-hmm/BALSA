CC = c++

all:	file db

file:	BALSAfile

db:	BALSAdatabase

CFLAGS = -g -O3

BALSAfile:	BALSA_file.cpp BALSA_file_get_seq.cpp BALSA_file_score_sums.cpp BALSA_file_align.cpp *.h
	$(CC) $(CFLAGS)  BALSA_file.cpp BALSA_file_get_seq.cpp BALSA_file_score_sums.cpp BALSA_file_align.cpp -o ../bin/BALSAfile
	@echo "../bin/BALSAfile created"

BALSAdatabase:	BALSA_database.cpp BALSA_database_score_sums.cpp BALSA_database_align.cpp *.h
	$(CC) BALSA_database.cpp BALSA_database_score_sums.cpp BALSA_database_align.cpp -o ../bin/BALSAdatabase
	@echo "../bin/BALSAdatabase created"
