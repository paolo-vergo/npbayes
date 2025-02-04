#include "Sampler.h"
#include <random>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>

using std::vector;
using std::size_t;
using std::cout;


void Sampler::sampling() {

    for (int i=0;i<n_iter;++i) {

        // Print progress
        if(i%2==0)
            cout << "iterazione: "<< i << "\n";

        // Refresh
        Refresh();

        // Print control
        if(i%200==0)
           Print();

        // Evaluate
        evaluate();        // update estimates for sigma and theta (Markov chains)


        // After burnin iterations (if set), start saving estimates

        if(i>=burnin) {
            for (size_t j = 0; j < K.size(); ++j)
                mu[i-burnin][j] = theta[i][K[j]];
            vector<std::size_t> temp=K;
            clustering.push_back(temp);
        }
    }
}

void Sampler::evaluate(){

    /* update parameters:
     * theta: vector<double> PUSH BACK*/
     vector<double> actual_theta;
     unsigned int M=0;
     for(size_t param=0;param<Nk.size();++param) {

        unsigned int n=0;
        double avg = 0.0;
        for (size_t i = 0; i < data.size(); ++i) {
            if (K[i] == param) {
                avg += std::accumulate(data[i].begin(), data[i].end(), 0.0);
                n += data[i].size();
            }
        }
        avg = avg / n;
        double mu_n = (n * tau * avg + sigma.back() * mu0) / (n * tau + sigma.back());
        double sigma_n = sigma.back() * tau / (n * tau + sigma.back());
        std::normal_distribution<double> d(mu_n, sigma_n);
        double value = d( gen);
        actual_theta.push_back(value);
        M += n ;
     }
     theta.push_back(actual_theta);

     double a_n = a + M;
     double ss = 0.0;
     for (size_t i=0; i<data.size(); ++i){
         double th = actual_theta[K[i]];
                 for (size_t j=0; j<data[i].size(); ++j){
                     ss += (data[i][j]-th)*(data[i][j]-th)/2.0;
                 }
     }
     double b_n= b + ss;
    std::gamma_distribution<double> s( a_n , 1.0/b_n);
    double val = 1.0/s(gen);
    sigma.push_back(val);

}


double Sampler::Ndf(double value, double MU, double SIGMA) const {
    double squared_term=sqrt(2*M_PI*SIGMA);
    double exp_term=exp(-1.0*(value-MU)*(value-MU)/(2*SIGMA));
    return exp_term/squared_term;
}


void Sampler::save_estimates() const {
    std::cout << "im printing\n";
    std::ofstream MyFile1("Mu_Estimates.csv");
    for (size_t j = 0; j < mu.size(); ++j) {
        for (size_t i = 0; i < data.size(); i++) {
            MyFile1 << mu[j][i];
            if (i != data.size() - 1)
                MyFile1 << ","; // No comma at end of line
        }
        MyFile1 << '\n';
    }

    std::ofstream MyFile2("Sigma_Estimates.csv");
    for (size_t j = 0; j < mu.size(); ++j) {
            MyFile2 << sigma[burnin + j];
            if (j != mu.size() - 1)
                MyFile2 << ","; // No comma at end of line
    }
    MyFile2 << '\n';

    std::ofstream MyFile3("Clust_Estimates.csv");
    for (size_t j = 0; j < clustering.size(); ++j) {
        for (size_t i = 0; i < data.size(); i++) {
            MyFile3 << clustering[j][i];
            if (i != data.size() - 1)
                MyFile3 << ","; // No comma at end of line
        }
        MyFile3 << '\n';
    }
}

void Sampler::Refresh() {
    //CHANGE TABLE
    for (size_t j = 0; j < dims.size(); ++j) {

        size_t restaurant=j;

        for (size_t i = acc_dims(restaurant); i < acc_dims(restaurant)+dims[restaurant]; ++i) {

            //EXTRACT FEATURES RELATIVE TO XIJ
            size_t index = i;                              // index of data item in vector DATA
            size_t table = T[index];                       // integer representing table in restaurant "restaurant"
            size_t dish  = K[index];                       // integer representing dish served at table "table"


            //COMPUTE CHANGE TABLE PROBABILITIES
            vector<double> weights = generate_weights_changetable(index, table,restaurant);
            std::discrete_distribution<int> d(weights.begin(), weights.end());

            //ASSIGN NEW TABLE
            size_t chosen_table = d(gen);

            //CASE 1: Table selected is the same as before
            if(chosen_table==table){
                //do nothing (?)
                if(i==debug)
                    std::cout <<" OK CASO 1"<<std::endl;
            }

            //CASE 2: New table is selected
            else if(chosen_table==Mjp[restaurant]){

                if(i==debug)
                    std::cout <<" start caso 2"<<std::endl;
                //SUB-CASE 2.1: Old table is empty when XIJ leaves it
                if(check_table_empty(table,restaurant,index)){
                    K[index]=set_new_K(index);                                   // Re-instantiate same table but change relative K;
                    if(Nk[dish]==1){
                        if(i==debug)
                            std::cout <<" start caso 2.1.1"<<std::endl;
                        //Update Nk[new_dish] and Mpk[new_dish], then erase Nk[dish] and Mpk[dish]
                        ++Nk[K[index]];
                        ++Mpk[K[index]];
                        Nk.erase(Nk.begin()+dish);
                        Mpk.erase(Mpk.begin()+dish);

                        //rescale K
                        for(size_t it=0;it<K.size();++it) {
                            if (K[it] > dish)
                                K[it] = K[it] - 1;
                        }
                        if(i==debug)
                            std::cout <<" ok caso 2.1.1"<<std::endl;
                    }
                    else {
                        if(i==debug)
                            std::cout <<" start caso 2.1.2"<<std::endl;
                        --Nk[dish];
                        ++Nk[K[index]];
                        --Mpk[dish];
                        ++Mpk[K[index]];
                        if(i==debug)
                            std::cout <<" ok caso 2.1.2"<<std::endl;
                    }
                }

                //SUB-CASE 2.2: Old table is not empty
                else{
                    if(i==debug)
                        std::cout <<" start caso 2.2"<<std::endl;
                    T[index]=chosen_table;
                    K[index]=set_new_K(index);
                    --Njt[restaurant][table];
                    Njt[restaurant].push_back(1);
                    ++Mjp[restaurant];
                    --Nk[dish];
                    ++Nk[K[index]];
                    ++Mpk[K[index]];
                    ++Mpp;
                    if(i==debug)
                        std::cout <<" OK caso 2.2"<<std::endl;
                }
                if(i==debug)
                    std::cout <<" OK caso 2"<<std::endl;
            }

            //CASE 3: Existent table, but not the old one, is selected
            else{
                if(i==debug)
                    std::cout <<" start caso 3"<<std::endl;
                //SUB-CASE 3.1: Old table is empty when XIJ leaves it
                if(check_table_empty(table,restaurant,index)){
                    if(i==debug)
                        std::cout <<" start caso 3.1"<<std::endl;
                    //Set new values
                    size_t chosen_dish=from_table_get_dish(chosen_table,restaurant);
                    T[index]=chosen_table;
                    K[index]=chosen_dish;

                    //Update Njt[restaurant][chosen_table], then erase Njt[restaurant][table]
                    ++Njt[restaurant][chosen_table];
                    Njt[restaurant].erase(Njt[restaurant].begin()+table);
                    --Mjp[restaurant];

                    //Rescale T vector
                    for(size_t it=acc_dims(restaurant);it<acc_dims(restaurant)+dims[restaurant];++it) {
                        if (T[it] > table)
                            T[it] = T[it] - 1;
                    }

                    //Update other quantities. Note: it is possible that a dish gets deleted
                    //SUB-SUB-CASE 3.1.1: no more people eat dish "dish"
                    if(Nk[dish]==1){
                        if(i==debug) {
                            std::cout << " start caso 3.1.1" << std::endl;
                        }
                        //Update Nk[new_dish] then erase Nk[dish] and Mpk[dish]
                        ++Nk[K[index]];
                        Nk.erase(Nk.begin()+dish);
                        Mpk.erase(Mpk.begin()+dish);
                        --Mpp;

                        //rescale K
                        for(size_t it=0;it<K.size();++it) {
                            if (K[it] > dish)
                                K[it] = K[it] - 1;
                        }
                        if(i==debug)
                            std::cout <<" OK caso 3.1.1"<<std::endl;
                    }

                    //SUB-SUB-CASE 3.1.2: table is empty, but there is another one (maybe in another restaurant) where people eat it
                    else{
                        if(i==debug)
                            std::cout <<" start caso 3.1.2"<<std::endl;
                        --Nk[dish];
                        ++Nk[K[index]];
                        --Mpk[dish];
                        --Mpp;
                        if(i==debug)
                            std::cout <<" OK caso 3.1.2"<<std::endl;

                    }
                    if(i==debug)
                        std::cout <<" OK caso 3.1"<<std::endl;
                }

                //SUB-CASE 3.2: Old table is not empty
                else{
                    if(i==debug)
                        std::cout <<" start caso 3.2"<<std::endl;
                    K[index]=from_table_get_dish(chosen_table,restaurant);
                    T[index]=chosen_table;

                    --Njt[restaurant][table];
                    ++Njt[restaurant][chosen_table];
                    --Nk[dish];
                    ++Nk[K[index]];
                    if(i==debug)
                        std::cout <<" OK caso 3.2"<<std::endl;
                }
                if(i==debug)
                    std::cout <<"OK caso 3"<<std::endl;
            }
        }
    }


    //CHANGE DISH

    for (size_t j = 0; j < dims.size(); ++j) {
        size_t restaurant=j;
        for (size_t i=0; i<Njt[restaurant].size();++i){

        // EXTRACT INFORMATIONS
        size_t table=i;
        size_t dish=from_table_get_dish(table,restaurant);

        // Assign new K to each table
        vector<double> weights = generate_weights_changedish(table,restaurant,dish);
        std::discrete_distribution<int> d(weights.begin(), weights.end());

        //New dish
        size_t chosen_dish= d(gen);

        //CASE 1: dish selected is the same as before
        if(chosen_dish==dish){
            //do nothing
        }

        //CASE 2: New dish selected
        else if(chosen_dish==Mpk.size()){

            //SUB-CASE 2.1: No more tables serve old dish
            if(Mpk[dish]==1){
                //do nothing
            }

            //SUB-CASE 2.2: Other tables keep serving old dish
            else{
                 // Add new dish
                 Nk[dish]-=Njt[restaurant][table];
                 Nk.push_back(Njt[restaurant][table]);
                 --Mpk[dish];
                 Mpk.push_back(1);

                 //Update k
                 for(size_t index=acc_dims(restaurant); index<acc_dims(restaurant)+dims[restaurant];++index){
                     if(K[index]==dish && T[index]==table)
                         K[index]=chosen_dish;
                 }
            }
        }

        //CASE 3: Chosen dish is already existent, but not the old one
        else{

            //SUB-CASE 3.1: No more tables serve old dish
            if(Mpk[dish]==1){
                Nk[chosen_dish]+=Njt[restaurant][table];
                ++Mpk[chosen_dish];
                Nk.erase(Nk.begin()+dish);
                Mpk.erase(Mpk.begin()+dish);

                //Update k
                for(size_t index=acc_dims(restaurant); index<acc_dims(restaurant)+dims[restaurant];++index){
                    if(K[index]==dish && T[index]==table)
                        K[index]=chosen_dish;
                }

                //rescale K
                for(size_t it=0;it<K.size();++it) {
                    if (K[it] > dish)
                        K[it] = K[it] - 1;
                }

            }

            //SUB-CASE 3.2: Other tables keep serving old dish
            else{
                // Add new dish
                Nk[dish]-=Njt[restaurant][table];
                Nk[chosen_dish]+=Njt[restaurant][table];
                --Mpk[dish];
                ++Mpk[chosen_dish];

                //Update k
                for(size_t index=acc_dims(restaurant); index<acc_dims(restaurant)+dims[restaurant];++index){
                    if(K[index]==dish && T[index]==table)
                        K[index]=chosen_dish;
                }
            }

        }
       }
    }
}


size_t Sampler::from_table_get_dish(size_t table_value,size_t restaurant_value) const {
    size_t ind=acc_dims(restaurant_value);
    while(T[ind]!=table_value)
        ++ind;
    return K[ind];
}

bool Sampler::check_table_empty(size_t table_value, size_t restaurant_value, size_t index) const {
    size_t value=Njt[restaurant_value][table_value];
    return (value==1);
}



size_t Sampler::set_new_K(size_t index) {
    vector<double> weights=generate_weights_assignK(index);
    std::discrete_distribution<int> d(weights.begin(),weights.end());
    size_t chosen_dish=d(gen);
    if(chosen_dish==Nk.size()) {
        Nk.push_back(0);
        Mpk.push_back(0);
    }
    return chosen_dish;
}



vector<double> Sampler::generate_weights_changetable(size_t index, size_t table, size_t restaurant) const{
        vector<double> weights;
        for (size_t i=0;i<Mjp[restaurant];++i){
            double value;
            size_t actual_dish=from_table_get_dish(i,restaurant);
            if(i==table)
                value=(Njt[restaurant][i]-1)*Fk(index,actual_dish);                 //Fk
            else
                value=(Njt[restaurant][i])*Fk(index,actual_dish);
            weights.push_back(value);
        }
        double value=alpha*Ftnew(index);                                            //Ftnew
        weights.push_back(value);

        return weights;
}

vector<double> Sampler::generate_weights_changedish(size_t table, size_t restaurant,size_t dish) const{
    vector<double> weights;                                             //Structure: implementation must be verified
    for (size_t k=0;k<Mpk.size();++k){
        double value;
        if(k==dish)
            value=(Mpk[k]-1)*Fdish(restaurant, table, k);               //Fdish
        else
            value=(Mpk[k])*Fdish(restaurant, table, k);
        weights.push_back(value);
    }
    double value=gamma*Fdishnew( table, restaurant);                    //Fdishnew
    weights.push_back(value);
    return weights;
}

vector<double> Sampler::generate_weights_assignK(size_t index) const{
    vector<double> weights;                                          //Structure: implementation must be verified
    for (size_t i=0;i<Mpk.size();++i){
        double value=(Mpk[i])*Fk(index,i);                               //Fk
        weights.push_back(value);
    }
    double value=gamma*Fknew(index);                                    //Fknew
    weights.push_back(value);
    return weights;
}


double Sampler::Fk(size_t index, size_t dish) const{

    vector<size_t> ids = get_index_from_k(dish, index);
    unsigned int count = 0;
    for (size_t i=0; i<ids.size(); ++i)
        count += data[ids[i]].size();
    double sum_metoo = sum_from_ids(ids) + std::accumulate(data[index].begin(), data[index].end(), 0.0);
    double sum = sum_from_ids(ids) ;
    double avg = sum / count;
    double avg_metoo = sum_metoo / (count + data[index].size());
    double ss=0.0;
    for (auto el:ids)
        for(auto el1:data[el])
           ss+=el1*el1;
    double ss_metoo= ss;
    for (auto el: data[index])
        ss_metoo += el*el;

    double sigma0=sigma.back();
    unsigned int card= data[index].size();

    double exp_11=-1.0*ss_metoo/(2*sigma0)-mu0*mu0/(2*tau);
    double exp_21=((tau*(count+card)*(count+card)*avg_metoo*avg_metoo/sigma0)+(sigma0*mu0*mu0/tau)+(2*(count+card)*avg_metoo*mu0))/(2*((count+card)*tau+sigma0));

    double exp_12=-1.0*ss/(2*sigma0)-mu0*mu0/(2*tau);
    double exp_22=((tau*(count)*(count)*avg*avg/sigma0)+(sigma0*mu0*mu0/tau)+(2*count*avg*mu0))/(2*((count)*tau+sigma0));

    double squared_term = sqrt((count*tau+sigma0)/((count+card)*tau+sigma0))/(pow(sqrt(2*M_PI*sigma0),card));
    double exp_term=exp_11+exp_21-exp_12-exp_22;

    return squared_term*exp(exp_term);
}

double  Sampler::sum_from_ids(const vector<size_t> & ids) const{
    double sum=0.0;
    for(auto el : ids){
        for (auto el1: data[el] )
         sum+=el1;
    }
    return sum;
}

vector<size_t> Sampler::get_index_from_k(size_t k_noto, size_t id_noto) const{    //escludo XJI
    vector <size_t> found_ind;
    for(size_t ii=0;ii<K.size();++ii){
        if(K[ii]==k_noto && id_noto!=ii)
            found_ind.push_back(ii);
    }
    return found_ind;
}


double Sampler::Ftnew(size_t index) const{
    double sum=0.0;
    for(size_t k=0; k<Nk.size();++k){
        sum+=(Mpk[k]/(Mpp+gamma))*Fk(index,k);
    }
    sum+=(gamma/(Mpp+gamma))*Fknew(index);
    return sum;
}

double Sampler::Fknew(size_t index) const {
    unsigned int card = data[index].size();
    double avg = 0.0;
    double ss = 0.0;
    for (auto el:data[index]) {
        avg += el;
        ss += el*el;
    }
    avg = avg / card;

    double sigma0=sigma.back();

    double squared_term = sqrt(sigma0) / (pow(sqrt(2 * M_PI * sigma0), card) * sqrt(card * tau + sigma0));

    double exp_1 = -1.0 * ss / (2 * sigma0) - mu0 * mu0 / (2 * tau);
    double exp_2 = ((tau * card * card * avg * avg / sigma0) + (sigma0 * mu0 * mu0 / tau) + (2 * card * avg * mu0)) /
                   (2 * (card * tau + sigma0));
    double exp_term = exp_1 + exp_2;

    return squared_term * exp(exp_term);
}

double Sampler::Fdish(size_t restaurant,size_t table, size_t dish ) const{

    vector<size_t> ids_table=get_index_from_table(table,restaurant);
    vector<size_t> ids=get_index_from_k_vector(dish,table, restaurant);

    unsigned count=0;
    for (auto el : ids)
        count+= data[el].size();

    unsigned count_all= count;
    for (auto el : ids_table)
        count_all += data[el].size();

    double sum = sum_from_ids(ids);
    double sum_metoo = sum +sum_from_ids(ids_table);

    double avg= sum/count;
    double avg_metoo=sum_metoo/(count_all);

    double ss = 0.0;
    for (auto el:ids)
        for(auto el1:data[el])
            ss += el1*el1;
    double ss_metoo = ss;
    for (auto el:ids_table)
        for(auto el1:data[el])
            ss_metoo += el1*el1;

    double sigma0= sigma.back();

    double exp_11=-1.0*ss_metoo/(2*sigma0)-mu0*mu0/(2*tau);
    double exp_21=((tau*(count_all)*(count_all)*avg_metoo*avg_metoo/sigma0)+(sigma0*mu0*mu0/tau)+(2*(count_all)*avg_metoo*mu0))/(2*((count_all)*tau+sigma0));

    double exp_12=-1.0*ss/(2*sigma0)-mu0*mu0/(2*tau);
    double exp_22=((tau*(count)*(count)*avg*avg/sigma0)+(sigma0*mu0*mu0/tau)+(2*count*avg*mu0))/(2*((count)*tau+sigma0));

    double squared_term = sqrt(((count)*tau+sigma0)/((count_all)*tau+sigma0))/pow(sqrt(2*M_PI*sigma0),(count_all-count));
    double exp_term=exp_11+exp_21-exp_12-exp_22;

    return squared_term*exp(exp_term);
}

vector<size_t> Sampler::get_index_from_table(size_t t_value,size_t restaurant) const{
    vector <size_t> found_ind;
    for(unsigned ii=acc_dims(restaurant);ii<acc_dims(restaurant)+dims[restaurant];++ii){
        if(T[ii]==t_value)
            found_ind.push_back(ii);
    }
    return found_ind;
}

vector<size_t> Sampler::get_index_from_k_vector(size_t k_chosen, size_t exclude_table, size_t restaurant) const{
    vector <size_t> found_ind;
    for(size_t ii=0;ii<K.size();++ii){
        if(K[ii]==k_chosen && out_of_table(ii,exclude_table,restaurant))
            found_ind.push_back(ii);
    }
    return found_ind;
}

bool Sampler::out_of_table(size_t temp_index, size_t exclude_table, size_t restaurant) const {

    return ((T[temp_index]!=exclude_table) ||  (rest_from_id(temp_index)!=restaurant));

}

size_t Sampler::rest_from_id (size_t index) const {
    size_t start=0;
    size_t it=0;
    while(!(start<=index && start+dims[it]>index)){
        start+=dims[it];
        ++it;
    }
    return it;
}



double Sampler::Fdishnew(size_t table, size_t restaurant) const{
    vector<size_t> ids_table=get_index_from_table(table,restaurant);
    unsigned int card=0;
    for (auto el : ids_table)
        card += data[el].size();
    double avg=0.0;
    double ss=0.0;
    for(auto el:ids_table) {
       for (auto el1 : data[el]) {
           avg += el1;
           ss += el1 * el1;
       }
    }
    avg=avg/card;
    double sigma0 = sigma.back();
    double squared_term=sqrt(sigma0)/(pow(sqrt(2*M_PI*sigma0),card)*sqrt(card*tau+sigma0));
    double exp_1=-1.0*ss/(2*sigma0)-mu0*mu0/(2*tau);
    double exp_2=((tau*card*card*avg*avg/sigma0)+(sigma0*mu0*mu0/tau)+(2*card*avg*mu0))/(2*(card*tau+sigma0));
    double exp_term= exp_1+exp_2;
    return squared_term*exp(exp_term);
}

void Sampler::Initialize() {
    for (size_t i=0; i<data.size();++i){
        data_mean.push_back( std::accumulate(data[i].begin(),data[i].end(), 0.0)/data[i].size());
    }

    vector<vector<std::size_t >> indices(dims.size());

    for (size_t i=0; i<dims.size();++i) {
        vector<size_t> ind = kmeans_run(acc_dims(i), acc_dims(i) + dims[i], clust[i]);
        indices [i] = ind;
    }

    for (size_t j = 0; j < dims.size(); ++j) {
        for (size_t i = 0; i < dims[j]; ++i) {
            T.push_back(indices[j][i]);
            K.push_back(T.back() + std::accumulate(clust.begin(),clust.begin()+j,0));
         }
    }

    Mpp=std::accumulate(clust.begin(),clust.end(),0);
    vector<int> temp(Mpp, 1);
    Mpk=temp;
    Mjp=clust;
    for ( size_t j=0; j<dims.size();++j){
        vector<int> temp1(clust[j],0);
        Njt.push_back(temp1);
    }
    for (size_t j=0; j<dims.size();++j)
        for(size_t i=0; i<dims[j];++i)
            Njt[j][T[acc_dims(j)+i]]++;

    for (size_t i=0; i<Mpp; ++i)
        Nk.push_back(0);
    for(auto el:K)
        Nk[el]++;

    // mu initialization

    Initialize_mu();

    // sigma initialization

    sigma.push_back(sigma_start);

}




vector<size_t> Sampler::kmeans_run(int start, int end, int num_groups) const {

        vector<size_t> label;
        vector <size_t> old_label(end-start);

        bool cond=true;
        int maxit=100;

        int iter=0;

        for(int i=start;i<end;++i){
            label.push_back(rand() % num_groups);
        }
        while(cond) {
            vector <int> count(num_groups,0);
            vector <double> centroids(num_groups,0.0);
            for (size_t j = 0; j < label.size(); ++j) {
                centroids[label[j]] += data_mean[start+j];
                ++count[label[j]];

            }
            for (size_t kk = 0; kk < centroids.size(); ++kk) {
                if(centroids[kk]==0 && count[kk]==0) {
                    centroids[kk] = 5;
                }
                else
                    centroids[kk] /= count[kk];
            }

            cond = (!equal(label.cbegin(), label.cend(), old_label.cbegin())) & (iter<maxit);

            old_label=label;

            for(size_t jj=0;jj<label.size();++jj) {
                vector<double> dist(num_groups,0.0);
                for (size_t ss = 0; ss < num_groups; ++ss) {
                    dist[ss] = std::abs(centroids[ss] - data_mean[start + jj]);
                }
                label[jj] = std::min_element(dist.cbegin(), dist.cend())-dist.cbegin();
            }
            ++iter;
        }
        return label;
}

void Sampler::Print() const {
    cout<<" TABLE ALLOCATION: \n";
    for(auto el:T)
        cout << el<< " ";
    cout <<"\n";
    cout<<" DISH ALLOCATION: \n";
    for(auto el:K)
        cout << el<< " ";
    cout <<"\n";
    cout<<" TOTAL NUMBER OF TABLES: \n";
    cout << Mpp;
    cout <<"\n";
    cout<<" NUMBER OF TABLES IN EACH RESTAURANT: \n";
    for(auto el:Mjp)
        cout << el <<" ";
    cout <<"\n";
    cout<<"NUMBER OF TABLES SERVING EACH DISH: \n";
    for(auto el:Mpk)
        cout << el <<" ";
    cout <<"\n";
    cout<<"NUMBER OF CUSTOMERS EATING EACH DISH: \n";
    for(auto el:Nk)
        cout << el <<" ";
    cout <<"\n";
    cout<<"NUMBER OF CUSTOMERS IN EACH TABLE IN EACH RESTAURANT: \n";
    for(auto el: Njt) {
        for (auto el1:el)
            cout << el1 << " ";
        cout << "\n";
    }
}



