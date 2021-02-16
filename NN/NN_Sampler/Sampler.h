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
  vector<double> data;

  vector<int> dims;

  //SAMPLER PARAMETERS
  int n_iter=100;

  int burnin=10;

  double sigma0=0.5;                           // variance of the sampling model. HERE IS FIXED! (N_N conjugated model)

  double mu0=0;                                // P0 parameters (set to this by default)

  double tau=10;
  
  double alpha=1;                              // group specific dirichlet distribution parameter

  double gamma=1  ;                            // Base Dirichlet distribution parameter

  vector<int> clust={3,2,2};                   // Initial number of clusters in each season

  //CLUSTERING DATA STRUCTURES

  vector<size_t> T;

  vector<size_t> K;

  vector<vector<int>> Njt;                    //Number of customers at table t in restaurant j

  vector<int> Nk;                             //Number of customers eating dish k across all restaurants

  vector<int> Mjp;                            //Number of tables in restaurant j

  vector<int> Mpk;                            //Number of tables where dish k is served

  int Mpp=0;                                  //Total number of tables across all restaurants

  // DATA STRUCTURES FOR DENSITY ESTIMATION
  vector<double> grid;                        //Initialization through specific method

  vector<vector<double>> density_estimates;   //Dimensions set accordingly to grid initialization

  // DATA STRUCTURES FOR PARAMETERS ESTIMATES

  vector<vector<double>> theta;               // Each cell in theta contains the length-variable vector of the unique values for MU
  vector<vector<double>> mu;                  // Each cell contains actual mu estimate per each athlete
  vector<vector<size_t>> clustering;          // Each cell contains cluster label per each athlete

  // RANDOM NUMBER GENERATOR

  std::mt19937 gen;

  /*------------------------------------------PRIVATE METHODS-----------------------------------------------------*/

  // MCMC Algorithm tools

  void Refresh();

  void evaluate();

  //Initialization tools

  vector <size_t> kmeans_run(int start,int end,int num_groups) const;

  //Densities

  double Fk(size_t index, size_t dish) const;

  double Ftnew(size_t index) const;

  double Fknew(size_t index) const;

  double Fdish(size_t restaurant,size_t table, size_t dish ) const;

  double Fdishnew(size_t table, size_t restaurant) const;

  //Utilities for densities

  double sum_from_ids(const vector<size_t> & ids) const;
  vector<size_t> get_index_from_k(size_t k_noto, size_t id_noto) const;

  vector<size_t> get_index_from_table(size_t t_value,size_t restaurant) const;
  vector<size_t> get_index_from_k_vector(size_t k_chosen, size_t exclude_table, size_t restaurant) const;

  bool out_of_table(size_t temp_index, size_t exclude_table, size_t restaurant) const;

  size_t rest_from_id (size_t index) const;

  //Utilities

  size_t from_table_get_dish(size_t table_value,size_t restaurant_value) const;

  bool check_table_empty(size_t table_value, size_t restaurant_value) const;

  int acc_dims(size_t rest) const {
        return std::accumulate(dims.begin(),dims.begin()+rest,0);
    }


  //Weights generators for discrete probability distribution

  vector<double> generate_weights_changetable(size_t index, size_t table, size_t restaurant) const;

  vector<double> generate_weights_changedish(size_t table, size_t restaurant,size_t dish) const;

  vector<double> generate_weights_assignK(size_t index) const;

  // New K setter

  size_t set_new_K (size_t index);

  // DENSITY ESTIMATION FUNCTIONS

  void normalize_grid_estimates(){
      for(auto& el:density_estimates)
          for(auto& el1:el)
              el1/=static_cast<double>(n_iter-burnin);
  };

  vector<double> generate_weights_table1( size_t restaurant) const;

  double Fknew1(size_t index_grid) const;

  vector<double> generate_weights_assignK1() const;

  double Ndf(double value, double mu, double sigma) const;

  double predictive(size_t grid_index, size_t restaurant, const vector<double>& MU) const;


    /*------------------------------------------PUBLIC METHODS-----------------------------------------------------*/

public:

  //CONSTRUCTOR: Need vector<double> of data (all in one line) and vector<int> specifying dimensions of the groups (restaurants) in data
  Sampler(const vector<double>& data,const vector<int>& dims):
         data(data), dims(dims){
                               T.reserve(data.size());
                               K.reserve(data.size());
                               std::random_device rd;
                               gen.seed(rd());
   };

  //INITIALIZER: custom initialization with k-means. With different data need to be modified. General version should be implemented
  void Initialize();

  void Initialize_mu(){
        vector<double> zero_vec( data.size(), 0.0);
        for (size_t i=0; i< (n_iter - burnin);++i)
            mu.push_back(zero_vec);
  }

  //GRID INITIALIZER: Note that it MUST be invoked before Sampling method. "fining refers to total number of points in the grid
  void Initialize_grid(int fining){
      std::cout << "Created grid:\n";
      double inf_mio=*(std::min_element(data.begin(),data.end()));
      std::cout<<"inf: "<< inf_mio << '\n';
      double sup=*(std::max_element(data.begin(),data.end()));
      std::cout<<"sup: "<< sup << '\n';
      double h=(sup-inf_mio)/ static_cast<double>(fining-1);
      std::cout<< "h: " << h << '\n';
      grid.push_back(inf_mio);
      for (int i=1; i<fining;++i){
          double val=inf_mio+h*static_cast<double>(i);
          grid.push_back(val);
      }

      /*---------------------*/
      //once initialized the grid, must reserve space of data structure for density estimation
      for (size_t i=0;i<dims.size();++i){
          vector<double> vec(fining,0.0);
          density_estimates.push_back(vec);
      }
  }

  //SETTERS: Allow non-admin users to set parameters. Default initialization is made upon basic case.
  void set_burnin(int n){
      burnin=n;
  }
  void set_n_iter(int n) {
    n_iter=n;
  };

  void set_mu0(double m) {
    mu0=m;
  };

  void set_tau(double t) {
    tau=t;
  };

  void set_alpha(double al) {
    alpha=al;
  };

  void set_gamma(double gam) {
    gamma=gam;
  };

  void set_sigma0(double s0) {
        sigma0=s0;
  };

  void set_clust(const vector<int>& clust_) {
        if( clust_.size()!= dims.size())
            std::cerr << "Incorrect dimensions\n";
        clust=clust_;
  };

  //ACTUAL SAMPLING: iterate Refresh for n_iter times. Update MCMC density estimates, then evaluate grid.
  void sampling() ;

  void save_estimates() const;

  // PRINT
  void Print() const;

};


