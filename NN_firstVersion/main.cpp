#include "Data.h"
#include "Sampler.h"
#include <vector>
#include "data_parser.h"

using std::vector;
int main() {

    //Import Rdata
    vector<double> Data=read_data_csv(R"(C:\Users\aughi\Desktop\Hierarchical Bayesian Nonparametric models to smooth functional data\NN_DataGeneration\Data2.csv)");
    vector<int>    Dims=read_dims_csv(R"(C:\Users\aughi\Desktop\Hierarchical Bayesian Nonparametric models to smooth functional data\NN_DataGeneration\Dims2.csv)");

    //Instantiate sampler object
    Sampler S(Data,Dims);

    //Set parameters of the model. IMPORTANT NOTE: must be done before calling Initialize!
    S.set_burnin(0);
    S.set_alpha(100);
    S.set_gamma(1);
    S.set_sigma0(0.3);

    //Initialize
    S.Initialize();

    //Print initial state
    S.Print();

    //Define grid refining
    S.Initialize_grid(200);

    //MCMC
    S.sampling();

    //Save output files in the working directory
    S.save_estimates();

}
