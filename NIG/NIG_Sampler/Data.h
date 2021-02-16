#include <iostream>
#include <vector>

#ifndef SCRIPT1_DATA_H
#define SCRIPT1_DATA_H

using std::vector;

class Data {
private:
    vector<vector<double>> Dataset;
    std::size_t groups_num=0;
    vector<int> dims;
public:
    Data()= default;;
    void AddGroup(vector<double> &vec){
        Dataset.push_back(vec);
        dims.push_back(vec.size());
        groups_num++;
    }

    void PrintDataGroup(std::size_t i) const{
        std::cout << "Group size: "<< Dataset[i].size() <<'\n';
        for (double k : Dataset[i]) {
            std::cout << k;
        }
        std::cout<< '\n';
    }

    void PrintData() const{
        for(std::size_t i=0; i<groups_num; i++) {
            std::cout << "Data from group " << i << ": \n";
            PrintDataGroup(i);
        }
    };

    vector<double> GetData() const{
        vector<double> D;
        for(std::size_t i=0; i<groups_num; i++) {
             D.insert(D.end(),Dataset[i].begin(),Dataset[i].end());
        }
        return D;
    }

    vector<int> GetDims() const{
        return dims;
    }


};


#endif //SCRIPT1_DATA_H
