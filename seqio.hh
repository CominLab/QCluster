#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;


// SeqRecord class --- a fasta record, accno + sequence
class SeqRecord
{
private:
	string header, description, sequence, quality;
public:
	SeqRecord(const SeqRecord&);
	SeqRecord(string, string, string, string);
	~SeqRecord();
	string seq() const {return sequence;}
	string id() const {return header;}
	string desc() const {return description;}
	string qual() const {return quality;}
	int length() const {return sequence.length();}
	SeqRecord toupper() const; // capitalize the sequence
	SeqRecord rc() const; // reverse complement
};


// Sequentially read records from a FASTA file, avoiding the need to read 
// huge files into memory at once; before reading a sequence using the 
// next() method one needs to make sure the stream is in a good state using
// the good() method 
class RecordGenerator
{
private:
	ifstream ifs;
	// we make the following variables instance attributes in order
	// to allow for persistence between the call of the constructor 
	// and the subsequent call to the next() method
	string id, desc;
	bool is_good;
	int indice_riga;
public:
	RecordGenerator(string);
	~RecordGenerator();
	bool good(); // indicates whether the fasta file
	// is good for reading another sequence
	SeqRecord next(); // read the next sequence
};


// Function for reading the vector of sequences into memory.
// Vector to read the sequence records into is passed by reference,
// and it gets appended with the newly read records --- it allows to 
// cumulatively read sequences from several files.
void FastqRead(string fname, vector<SeqRecord>& recs);


// Function for counting reads in a fasta file
int CountReads(string fname);
