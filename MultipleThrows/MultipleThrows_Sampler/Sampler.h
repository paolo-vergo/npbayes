#include <vector>
#include <iostream>
#include <numeric>
#include <random>
#include <algorithm>

using std::vector;
using std::cout;

class Sampler{
private:

  //DATA STRUCTURES
  vector<vector<double>> data;           // each point in the data structure admits a variable_size

  vector<int> dims;                      // dimensions per each restaurants

  vector<int> per_athlete_dims;          // number of trials per each athlete

  vector<double> data_mean;              //vector of means per athlete. Only internal use for k_means initialization.

  //SAMPLER PARAMETERS

  int n_iter=1000;

  int burnin=0;

  double mu0=17.3;                                // P0 parameters (set to this by default)
  double tau=4.0;

  double a=10.0;                                  // IG parameters. Note that sigma2 is out of the mixture
  double b=9.0;

  double sigma_start=3.0;                         // Starting value for the vector of sigma estimates
  
  double alpha=1;                                 // Group specific dirichlet distribution parameter
  double gamma=1;                                 //  Base dirichlet distribution parameter

  vector<int> clust={3,3,3,3,3,3,3,3,3,3};        // Number of clusters per restaurant. Only internal use.
  //CLUSTERING DATA STRUCTURES

  vector<size_t> T;                               // Vector of table indices per each customer

  vector<size_t> K;                           //Vector of dish indices per each customer

  vector<vector<int>> Njt;                    //Number of customers at table t in restaurant j

  vector<int> Nk;                             //Number of customers eating dish k across all restaurants

  vector<int> Mjp;                            //Number of tables in restaurant j

  vector<int> Mpk;                            //Number of tables where dish k is served

  int Mpp=0;                                  //Total number of tables across all restaurants

  // DATA STRUCTURES FOR PARAMETERS ESTIMATES

  vector<double> sigma;                       // Sigma estimates
  vector<vector<double>> theta;               // Each cell in theta contains the length-variable vector of the unique values for MU
  vector<vector<double>> mu;                  // Each cell contains actual mu estimate per each athlete
  vector<vector<size_t>> clustering;          // Each cell contains cluster label per each athlete


  // RANDOM NUMBER GENERATOR

  std::mt19937 gen;

  // Debug Mode

  size_t debug=data.size()+1;                  // Non active; if an existent data index is provided, show additional info during run. May cause increase of time.

  /*------------------------------------------PRIVATE METHODS-----------------------------------------------------*/

  // MCMC Algorithm tools

  void Refresh();                           // change table allocation per each customer, then dish allocation per each table

  void evaluate();                          // given actual allocation, compute all quantities of interest

  // Marginals

  double Fk(size_t index, size_t dish) const;

  double Ftnew(size_t index) const;

  double Fknew(size_t index) const;

  double Fdish(size_t restaurant,size_t table, size_t dish ) const;

  double Fdishnew(size_t table, size_t restaurant) const;

  //Utilities for marginals computation

  double sum_from_ids(const vector<size_t> & ids) const;

  vector<size_t> get_index_from_k(size_t k_noto, size_t id_noto) const;

  vector<size_t> get_index_from_table(size_t t_value,size_t restaurant) const;

  vector<size_t> get_index_from_k_vector(size_t k_chosen, size_t exclude_table, size_t restaurant) const;

  bool out_of_table(size_t temp_index, size_t exclude_table, size_t restaurant) const;

  size_t rest_from_id (size_t index) const;

  //Utilities

  size_t from_table_get_dish(size_t table_value,size_t restaurant_value) const;

  bool check_table_empty(size_t table_value, size_t restaurant_value, size_t index_value) const;

  int acc_dims(size_t rest) const {
        return std::accumulate(dims.begin(),dims.begin()+rest,0);
    }

  double Ndf(double value, double MU, double SIGMA) const;

  //Weights generators for discrete probability distribution

  vector<double> generate_weights_changetable(size_t index, size_t table, size_t restaurant) const;

  vector<double> generate_weights_changedish(size_t table, size_t restaurant,size_t dish) const;

  vector<double> generate_weights_assignK(size_t index) const;

  // New dish setter

  size_t set_new_K(size_t index);


    /*------------------------------------------PUBLIC METHODS-----------------------------------------------------*/

public:

  //CONSTRUCTOR: Need vector<double> of data (all in one line) and vector<int> specifying dimensions of the groups (restaurants) in data.
  //             per_athlete_dims specifies number of observations related to each athlete

  // IMPORTANT NOTE: Data to constructor must be provided as a one_line vector

  Sampler(const vector<double>& Data,const vector<int>& dims, const vector<int>& per_athlete_dims):
          dims(dims), per_athlete_dims(per_athlete_dims) {
                               for (size_t j=0; j<per_athlete_dims.size(); ++j){
                                   int length=per_athlete_dims[j];
                                   vector<double> vec(length,0.0);
                                   data.push_back(vec);
                                   int res=0;
                                   for (size_t k=0;k<j;++k) {
                                       res += per_athlete_dims[k];
                                   }
                                   for (size_t i=0; i<length;++i) {
                                           data[j][i] = Data[res + i];
                                   }
                               }
                               T.reserve(data.size());
                               K.reserve(data.size());
                               std::random_device rd;
                               gen.seed(rd());
   };

  Sampler(const vector<double>& Data,const vector<int>& dims, const vector<int>& per_athlete_dims, unsigned int seed):
            dims(dims), per_athlete_dims(per_athlete_dims) {
        for (size_t j=0; j<per_athlete_dims.size(); ++j){
            int length=per_athlete_dims[j];
            vector<double> vec(length,0.0);
            data.push_back(vec);
            int res=0;
            for (size_t k=0;k<j;++k) {
                res += per_athlete_dims[k];
            }
            for (size_t i=0; i<length;++i) {
                data[j][i] = Data[res + i];
            }
        }
        T.reserve(data.size());
        K.reserve(data.size());
        gen.seed(seed);
  };

  //INITIALIZER: custom initialization with k-means. With different data need to be modified. General version should be implemented

  void Initialize();

  //Initialization tools

  void Initialize_mu(){
        vector<double> zero_vec( data.size(), 0.0);
        for (size_t i=0; i< (n_iter - burnin);++i)
            mu.push_back(zero_vec);
    }

  vector <size_t> kmeans_run(int start,int end,int num_groups) const;


  //SETTERS: Allow non-admin users to set parameters. Default initialization is made upon base case.
  void set_burnin(int n){
      burnin=n;
  }
  void set_n_iter(int n) {
    n_iter=n;
  };

  void set_mu0(double mu) {
    mu0=mu;
  };

  void set_tau(double tau) {
    tau=tau;
  };

  void set_alpha(double al) {
    alpha=al;
  };

  void set_gamma(double gam) {
        gamma=gam;
  };

  void set_a(double a_) {
    a=a_;
  };

  void set_b(double b_) {
        b=b_;
  };

  void set_sigma_start(double sig) {
        sigma_start=sig;
  };

  void set_clust(const vector<int>& clust_) {
        if( clust_.size()!= dims.size())
            std::cerr << "Incorrect dimensions\n";
        clust=clust_;
  };

  void set_debug(size_t deb){
      if(deb>data.size())
          std::cerr<< "Debug Mode already set to Off\n";
      else
          debug=deb;
  }

  //ACTUAL SAMPLING:

  void sampling() ;

  void save_estimates() const;

  // PRINT

  void Print() const;

  void Print_data() const{
      cout << " Data:\n";
      for (auto el : data){
          for (auto el1:el){
              cout << el1 << " ";
          }
          cout << '\n';
      }

      cout << " Dims:\n";
      for(auto el : dims)
          cout << el << " ";
      cout << '\n';

      cout << " Per_athlete_dims:\n";
      for(auto el : per_athlete_dims)
          cout << el << " ";
      cout << '\n';

  }

};


