#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>



/*
 * a class to read data from a csv file.
 */
class csvreader
{
	std::string filename;
	std::string delimeter;

public:
	csvreader(std::string filename, std::string delm = ",");
	~csvreader();
	// function to fetch data from a csv file
	std::vector<std::vector<std::string>> getdata();
};



