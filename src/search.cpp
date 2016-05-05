#include <Rcpp.h>
#include <stdio.h>

using namespace std;
using namespace Rcpp;

int max(SEXP g) {
  NumericVector v(g);
  //The minimum possible value is 1, so max = min-1
	int max = 0;
	
	for (int i = 0; i < v.size(); i++)
		if(v[i] > max)
			max = v[i];
  
	return max; 
}


//precondition: a list of positive numbers
//[[Rcpp::export]]
SEXP normalizeVector (SEXP v) 
{    
  NumericVector vec(v);
  
  if (vec.size() <= 0) 
    Rcpp::stop("Empty sample");
    
  int max_seen = 0;
  
  int max_vec = max(vec);
  
  int x = 0;
  
  for (int i = 0; i < vec.size(); i++) {
    x = vec[i];
    
    if (x > max_seen+1) {
      
      //the start should be the ith element
      std::replace (vec.begin(), vec.end(), max_seen+1, ++max_vec);
      
      std::replace (vec.begin(), vec.end(), x, max_seen+1);
    }
    
    if (vec[i] > max_seen)
        max_seen = vec[i];
  }
  
  return(vec);
}

double computeMQClusterFactor(double m, double n) {
    if (m > 0) return ((2 * m) / ((2 * m) + n));
    return(0);
}
  
double computeLQClusterFactor(double m, double n) {
  if (m > 0) return(m / (m + (0.5 *  n)));
    return(0);
}
  

//The mydata is a symmetric matrix
//[[Rcpp::export]]
SEXP evaluateLQ(SEXP d, SEXP g) {

	NumericMatrix mydata(d);
	NumericVector group(g);
	
	int k = max(group);
	//if (k <= 0)
	//	Rcpp::stop("wrong number of clusters");
	
  vector<vector<int> > group_indices(k);	
  
  for (int i=0; i < group.size(); i++)
    group_indices[group[i]-1].push_back(i);
		
	vector<double> intra_sim_dis(k);
	vector<double> inter_sim_dis(k);
  vector<double> all_sim_dis(k);
	
	for (int i =0; i<k;i++){
		intra_sim_dis[i] = 0;
		all_sim_dis[i] = 0;
	}
	
  int indices_size;
  
  int row_index;
  //the diagonal of submatrix = 0, so exclude
  for (int m=0; m < k; m++) {
    indices_size = group_indices[m].size();
    //start from 1
    for (int i=0; i < indices_size; i++){
      row_index = group_indices[m][i];
      for (int j = i+1; j < indices_size; j++)
        intra_sim_dis[m] += mydata(row_index, group_indices[m][j]);  
    }
  }
  
  for (int m =0; m < k; m++) {
    indices_size = group_indices[m].size();
    //start from 1
    for (int i=0; i < indices_size; i++)      
        all_sim_dis[m] += sum(mydata(group_indices[m][i], _));

  }
  
  for (int m =0; m < k; m++)    
    inter_sim_dis[m] = all_sim_dis[m] - 2 * intra_sim_dis[m];

  
	double sum = 0;
	
	for(int i=0; i < k; i++)
		sum += computeLQClusterFactor(intra_sim_dis[i], inter_sim_dis[i]);
	
   //cout << intra_sim_dis[2] << ' ' << inter_sim_dis[2];
  
	return wrap(sum); 

}

//[[Rcpp::export]]
SEXP evaluateMQ(SEXP d, SEXP g) {
  
  NumericMatrix mydata(d);
	NumericVector group(g);
	
	int k = max(group);
//	if (k <= 0)
//		Rcpp::stop("wrong number of clusters");
	
  vector<vector<int> > group_indices(k);	
  
  for (int i=0; i < group.size(); i++)
    group_indices[group[i]-1].push_back(i);
		
	vector<int> intra_edges(k);
	vector<int> inter_edges(k);
  vector<int> all_edges(k);
	
	for (int i =0; i<k;i++){
		intra_edges[i] = 0;
		all_edges[i] = 0;
	}
	
  int indices_size;
  
  int row_index;
  //the diagonal of submatrix = 0, so exclude
  for (int m=0; m < k; m++) {
    indices_size = group_indices[m].size();
    //start from 1
    for (int i=0; i < indices_size; i++){
      row_index = group_indices[m][i];
      for (int j = 0; j < indices_size; j++)
      intra_edges[m] += mydata(row_index, group_indices[m][j]);  
      
    }
  }
  
  for (int m =0; m < k; m++) {
    indices_size = group_indices[m].size();
    //start from 1
    for (int i=0; i < indices_size; i++)      {
      all_edges[m] += sum(mydata(group_indices[m][i], _));
      all_edges[m] += sum(mydata(_, group_indices[m][i]));
    }
  }
  
  for (int m =0; m < k; m++)    
    inter_edges[m] = all_edges[m] - 2 * intra_edges[m];

  
	double sum = 0;
	
	for(int i=0; i < k; i++)
		sum += computeMQClusterFactor(intra_edges[i], inter_edges[i]);
	
   //cout << intra_sim_dis[2] << ' ' << inter_sim_dis[2];
  
	return wrap(sum); 

}

/* parent: the current partitions, for which neighbors are to be computed
 * k: number of partitions in parent
 * percentage: the percentage of neighbors to compute
 * columns: number of elements in the parent
 */

/*
//[[Rcpp::export]]
SEXP computeNearestAscent(SEXP parent, int columns, int k, double percentage) {
  
  NumericVector group(parent);
  
  vector<NumericVector> neighbors;

  neighbors.attr("dimnames") = List::create(CharacterVector(n), group.attr("names"));
  
  int total = 0;
    
  for (int i=0; i < columns; i++) {
    
    int current_column = group[i];
    
    int iterations = ceil((k-1) * percentage)
    
    random_clusters = generateUniqueNumbers(k, iterations, current_column);
    
    for (int j =0; j < iterations; j++) {
      
      NumericVector group_clone(clone(group));
    
      group_clone[current_column] = random_clusters[j];
    
      neighbors[total++].push_back(group_clone);

    }
  
  } 
  
  return(wrap(neighbors));
  
}
*/

/* n: max
 * m: m random numbers
 * x: excluding this number
 * start from 1 to including n
*/
/*
int[] generateUniqueNumbers(int n, int m, int x){
  int arr[n];
  int res[m];
  
  int idx = 0;
  while ((idx < n) && (idx != (x-1))) {
    arr[idx] = idx +1;
    idx++;
  }

  for (int i =0; i < m; ++i) {
    int j = rand(idx-i); // Pick random number from 0 <= r < n-i.  Pick favorite method
    // j == 0 means don't swap, otherwise swap with the element j away
    if (j != 0) { 
        std::swap(arr[i], arr[i+j]);
    }
  }
  
  for (int i =0; i < m; ++i)
    res[i] = arr[idx-i-1];
  
  return(res);
  
}

*/

//[[Rcpp::export]]
NumericVector computeRandomNeighbor(NumericVector group, int columns, int k) {
   
   NumericVector group_clone(clone(group));
    
   int column_index = rand() % columns;
        
   int column_value = rand() % (k+1) + 1;
        
   group_clone[column_index] = column_value;
        
   return(normalizeVector(group_clone));
 
}


//[[Rcpp::export]]
SEXP computeRandomNeighbors(NumericVector group, int columns, int k, int n) {
    
    NumericMatrix neighbors(n, columns);

    neighbors.attr("dimnames") = List::create(CharacterVector(n), group.attr("names"));
    
    for (int i=0; i < n; i++) 
      neighbors(i, _) = computeRandomNeighbor(group, columns, k);
          
    return(neighbors);
    
  }
    


