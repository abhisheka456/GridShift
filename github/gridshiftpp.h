#include <map> 
#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;


void generate_offsets_cy(int d,
                         int base,
                         int * offsets) {
    /*
        Generate 3**d neighbors for any point.

        Parameters
        ----------
        d: Dimensions
        base: 3, corresponding to (-1, 0, 1)
        offsets: (3**d, d) array of offsets to be added to 
                 a bin to get neighbors

    */

    int tmp_i;

    for (int i = 0; i < pow(base, d); i++) {
        tmp_i = i;
        for (int j = 0; j < d; j++) {
            if (tmp_i == 0) break;
            offsets[i * d + j] = tmp_i % base - 1;
            tmp_i /= base;
        }
    }
}


void grid_cluster(int n,
                  int d,
                  int base,
                  int iterations,
                  float bandwidth,
                  int * offsets,
                  float * X_shifted,
                  int * membership, 
                  int * k_num) {
                  
    map< vector<int>, pair< vector<float>, int> > cluster_grid;
    map< vector<int>, int > map_cluster;
    map< int, int > clus;
    map< int, int > :: iterator it2;
    map< vector<int>, pair< vector<float>, int> >:: iterator it;
    map< vector<int>, pair< vector<float>, int> > means;
    

    int iter = 0;
    vector<int> current_bin(d);
    vector<int> bin(d);
    vector<int> membershipp(n);
    vector<int> membershipp_old(n);

    // new clustering at grids
    int temp = 0;
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < d; k++) {
            bin[k] = X_shifted[i * d + k] / bandwidth;
        }
        if (cluster_grid.find(bin) == cluster_grid.end()) {
            cluster_grid[bin] = make_pair(std::vector<float>(d, 0), 0);
        }
        for (int k = 0; k < d; k++){
            cluster_grid[bin].first[k] += X_shifted[i * d + k];
        }
        cluster_grid[bin].second++;
        if (map_cluster.find(bin) == map_cluster.end()){
            map_cluster[bin] = temp++;
        }
        membershipp[i] = map_cluster[bin] * 1.0;
    }
    
    
    //copy(membershipp.begin(), membershipp.end(), membership);              
    while (iter <= iterations){
        iter++;
        means.clear();
        for (it = cluster_grid.begin(); it != cluster_grid.end(); ++it ){
        
            for (int j = 0; j < pow(base, d); j++) {
            
                for (int k = 0; k < d; k++) {
                
                    current_bin[k] = it->first[k] + offsets[j * d + k];
                    
                    if (j == 0){
                    
                        bin[k] =  it->first[k] ;
                        
                    }
                     
                 
                }
            
                // If neighbor exists, add it to the mean
                if (cluster_grid.find(current_bin) != cluster_grid.end()) {
                
                    if (means.find(current_bin) == means.end()) {
                        means[current_bin] = make_pair(std::vector<float>(d, 0), 0);
                    }
                
                    for (int k = 0; k < d; k++) {
                        means[current_bin].first[k] += cluster_grid[bin].first[k] * 1.0;
                    }

                    means[current_bin].second += cluster_grid[bin].second;
                    }
                }
            
            }
        
         for (it = cluster_grid.begin(); it != cluster_grid.end(); ++it ){
            for (int k = 0; k < d; k++) {
                current_bin[k] = it->first[k];
            }
        
            for (int k = 0; k < d; k++) {
                cluster_grid[current_bin].first[k] = means[current_bin].first[k] * 1.0 / means[current_bin].second;
            }
        
        }
       
        // update cluster grid and membership 
        map< vector<int>, pair< vector<float>, int> > cluster_grid_old = cluster_grid;
        map< vector<int>, int > map_cluster_old = map_cluster;

        cluster_grid.clear();
        map_cluster.clear();
        clus.clear();
        

        int temp = 0;
        for (it = cluster_grid_old.begin(); it != cluster_grid_old.end(); ++it ){

            for (int k = 0; k < d; k++) {
                bin[k] = it->second.first[k] / bandwidth;
                current_bin[k] = it->first[k];
            }
            
            if (cluster_grid.find(bin) == cluster_grid.end()) {
                cluster_grid[bin] = make_pair(std::vector<float>(d, 0), 0);
            }   
            
            for (int k = 0; k < d; k++){
                cluster_grid[bin].first[k] += it->second.first[k] * 1.0 * it->second.second;
            }
            cluster_grid[bin].second += it->second.second;
            
            if (map_cluster.find(bin) == map_cluster.end()){
                map_cluster[bin] = temp++;
            }

            //replace (membershipp.begin(), membershipp.end(), map_cluster_old[current_bin], map_cluster[bin]);
            clus[map_cluster_old[current_bin]] = map_cluster[bin];
        }
        //membershipp_old = membershipp;
        int break_points = 0;
        for (it2 = clus.begin(); it2 != clus.end(); ++it2) {
        if (it2->first !=  it2->second){
            replace (membershipp.begin(), membershipp.end(), it2->first, it2->second);
            break_points += it2->first - it2->second;
            }
           
        }
        //if (membershipp_old == membershipp) {
        if (break_points == 0){
            break;
        }
    
    }
    copy(membershipp.begin(), membershipp.end(), membership);  
    
    vector<int> k_num2(1);
    k_num2[0] = cluster_grid.size();
    vector<float> bins(k_num2[0] * d);
    int itt = 0;
     for (it = cluster_grid.begin(); it != cluster_grid.end(); ++it ){
            
            for (int k = 0; k < d; k++) {
                bins[itt * d + k] = it->second.first[k] *1.0 / it->second.second;

            }
            itt++;
        }
    copy(bins.begin(), bins.end(),X_shifted);  
    copy(k_num2.begin(), k_num2.end(),k_num); 
}
