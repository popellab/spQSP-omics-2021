#include "csvreader.h"

/*
* parses through csv file line by line and returns the data
* in vector of vector of strings.
*/

csvreader::csvreader(std::string file, std::string delm)
	: filename{ file }, delimeter{ delm } {}

csvreader::~csvreader() {}

std::vector<std::vector<std::string> > csvreader::getdata()
{
	std::ifstream file(filename);

	if (!file) 
	{
		std::cerr << "Antigen File Not Found" << std::endl;
		exit(1);
	}

	std::vector<std::vector<std::string>> datalist;

	std::string line{ "" };
	// iterate through each line and split the content using delimeter

	for (int i = 0; i < 4; i++) // ignore lines that contains other information (gene name, index, gene id etc.) 
	{							// total four lines
		getline(file, line);
	}

	while (getline(file, line))
	{
		std::vector<std::string> vec;
		boost::algorithm::split(vec, line, boost::is_any_of(delimeter));
		datalist.push_back(vec);
	}
	// close the file
	datalist.shrink_to_fit();
	file.close();

	return datalist;
}
