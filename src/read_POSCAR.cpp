#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>

using namespace std;

/*
###*******************************DOCUMENTATION********************************
### AUTHOR: Asif Iqbal Bhatti
### DATE : 01/11/2020
### USAGE: icpc/icc/g++ defor.cpp (Intel/gcc compiler)
### PURPOSE: Reads only a POSCAR (VASP) file.
###****************************END OF DOCUMENTATION****************************
*/
template <class T>
bool from_string(T& t, const std::string& s, std::ios_base& (*f)(std::ios_base&))
  {
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
  }

int main()
{
	int array_size = 100000;
	string array[array_size];
	std::string filename, line;
	std::ifstream indata;
	float latvec[4][4] = {}, xx;
	int loop = 0, loop2;
	int p = 0;
	
	cout << "__|  Enter the POSCAR filename : ";
	getline(cin, filename);
	indata.open(filename.c_str(), ios::in);
	
	if (indata.is_open())
		{
		while (! indata.eof() )
			{
			getline (indata,line);
	//		cout << line <<  endl;
			array[loop]=line;
			loop++;
			}
			cout <<"# of lines:   " << loop << endl;
			indata.close();
		}
	else cout << "Unable to open file";
	
	cout << array[0] <<endl;
	cout << array[1] <<endl;
	
	for (loop2=2; loop2<5;loop2++)	
	{
		int i = loop2 - 1, j = 1;
		stringstream iss(array[loop2]);
		string buf;
		vector<string> tokens;
		while (iss >> buf)
		{
			tokens.push_back(buf);
			if(from_string<float>(xx, std::string(buf), std::dec))		
			{
				latvec[i][j] = xx;
				printf("%16.12f\t", latvec[i][j]* ::atof(array[1].c_str()));
	// 			cout << i << "\t" << j << endl;
				j++;
			}
		}	
		std::cout << endl;
	}
	cout << array[5] <<endl;
	cout << array[6] <<endl;
	cout << array[7] <<endl;
	
	for (loop2=8; loop2 < loop ;loop2++)
		{
		cout << array[loop2] <<endl;
		}
	return 0;
}

