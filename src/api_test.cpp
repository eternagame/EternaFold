#include "contrafold_api.hpp"
#include <string>
#include <iostream>
#include <string.h>

void test_pfunc();
// void test_predict_struct();
// void test_c_fold();
// void test_c_energy_of_structure();

void test_pfunc(){
	char seq[] = "GGGGGAAAAAACCCCC";
	char structure[] = "(((((......)))))";
	float Z = pfunc(seq,structure);
	std::cout << "Log partition coefficient: " << Z << std::endl;

}

void test_predict_struct(){
	char seq[] = "GGGGGAAAAAACCCCC";
	char* pred_struct = new char[strlen(seq)];
	pred_struct = predict_struct(seq);
	std::cout << "MEA predicted struct: " << pred_struct << std::endl;
}

void test_c_fold(){
	char seq[] = "GGGGGAAAAAACCCCC";
	char* structure = new char[strlen(seq)];
	float energy = c_fold(seq, structure);
	std::cout << "Energy: " << energy << std::endl;
	std::cout << "MEA predicted struct: " << structure << std::endl;
}

void test_c_energy_of_structure(){
	char seq[] = "GGGGGAAAAAACCCCC";
	char structure[] = "(((((......)))))";
	float energy = c_energy_of_structure(seq, structure);
	std::cout << "Energy: " << energy << std::endl;
}

int main()
{
	test_pfunc();
	test_predict_struct();
	test_c_fold();
	test_c_energy_of_structure();
	return 0;
}
