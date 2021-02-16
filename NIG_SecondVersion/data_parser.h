#ifndef SCRIPT1_DATA_PARSER_H
#define SCRIPT1_DATA_PARSER_H

#include <string>
#include <fstream>
#include <vector>
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream

std::vector<double> read_data_csv(std::string filename){
    // Reads a data CSV file into a vector of <double>

    // Create a vector of <double> pairs to store the result
    std::vector<double> result;

    // Create an input filestream
    std::ifstream myFile(filename);

    // Make sure the file is open
    if(!myFile.is_open()) throw std::runtime_error("Could not open file");

    // Helper vars
    std::string line;
    double val;


    // Read data, line by line
    while(std::getline(myFile, line))
    {
        // Create a stringstream of the current line
        std::stringstream ss(line);

        // Extract each double
        while(ss >> val){

            // Add the current integer to the 'colIdx' column's values vector
            result.push_back(val);

            // If the next token is a comma, ignore it and move on
            if(ss.peek() == ',') ss.ignore();
        }
    }

    // Close file
    myFile.close();

    return result;
}

std::vector<int> read_dims_csv(std::string filename){
    // Reads dims CSV file into a vector of <int>

    // Create a vector of <double> pairs to store the result
    std::vector<int> result;

    // Create an input filestream
    std::ifstream myFile(filename);

    // Make sure the file is open
    if(!myFile.is_open()) throw std::runtime_error("Could not open file");

    // Helper vars
    std::string line;
    int val;


    // Read data, line by line
    while(std::getline(myFile, line))
    {
        // Create a stringstream of the current line
        std::stringstream ss(line);

        // Extract each double
        while(ss >> val){

            // Add the current integer to the 'colIdx' column's values vector
            result.push_back(val);

            // If the next token is a comma, ignore it and move on
            if(ss.peek() == ',') ss.ignore();
        }
    }

    // Close file
    myFile.close();

    return result;
}

#endif //SCRIPT1_DATA_PARSER_H
