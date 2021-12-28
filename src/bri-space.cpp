#include<iostream>

#include "br_index.hpp"

using namespace std;
using namespace bri;

void help(){
	cout << "bri-space: breakdown of index space usage" << endl;
	cout << "Usage:       bri-space <index>" << endl;
	cout << "   <index>   index file (with extension .bri)" << endl;
	exit(0);
}

int main(int argc, char** argv){

	if(argc != 2)
		help();

	br_index<> idx;
	idx.load_from_file(argv[1]);

	auto space = idx.print_space();

	cout << "TOTAL space: " << space << " bytes" << endl;


}