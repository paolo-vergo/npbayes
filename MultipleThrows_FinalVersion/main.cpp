#include "Data.h"
#include "Sampler.h"
#include <vector>
#include "data_parser.h"

using std::vector;
int main() {
    //Import data
    vector<double> Data=read_data_csv(R"(C:\Users\aughi\Desktop\Hierarchical Bayesian Nonparametric models to smooth functional data\Shotput Dataset\Pre processing\Data\Data.csv)");
    vector<int>    Dims=read_dims_csv(R"(C:\Users\aughi\Desktop\Hierarchical Bayesian Nonparametric models to smooth functional data\Shotput Dataset\Pre processing\Data\Dims.csv)");
    vector<int>    Per_athlete_dims=read_dims_csv(R"(C:\Users\aughi\Desktop\Hierarchical Bayesian Nonparametric models to smooth functional data\Shotput Dataset\Pre processing\Data\per_athlete_dims.csv)");

    //Instantiate sampler object
    Sampler S(Data,Dims,Per_athlete_dims);

    //Initialize: eventually set parameters using setters
    S.Initialize();

    //Print initial state
    S.Print();

    //MCMC
    S.sampling();

    //Final state
    S.Print();

    //Produce output files in the working directory
    S.save_estimates();

}
