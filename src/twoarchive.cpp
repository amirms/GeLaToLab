#include <Rcpp.h>
#include <stdio.h>

using namespace std;
using namespace Rcpp;

vector<List> getResults(List collection) {
  
  vector<List> results;
  
  for(int i =0; i < collection.size(); i++) 
    results.push_back(((List)collection[i])["result"]);
  
  return(results);
  
  
}


int all(vector<int> lelements) {
  for(vector<int>::size_type i = 0; i < lelements.size(); i++)
    if (lelements[i] == 0 )
      return(0);
    
  return(1);
  
}

int any(vector<int> lelements) {
  for(vector<int>::size_type i = 0; i < lelements.size(); i++)
      if (lelements[i] > 0)
        return(1);
      
  return(0);
  
}

int equal(List r1, List r2) {
 
  for (int i=0; i < r1.size(); i++){

    if (((NumericVector)r1[i])[0] != ((NumericVector)r2[i])[0])
      return(0);
      
  }
  return(1);    
  
}

//Domination happend iff there is at least one improvement, that is, 
// f_i(x) > f_i(y) 

//But let's get rid of equal elements
int dominates(List r1, List r2) {
  if (r1.size() != r2.size())
    Rcpp::stop("sizes do not match")  ;
  
  if(equal(r1, r2))
    return(0);

  for (int i=0; i < r1.size(); i++)
    if (((NumericVector)r1[i])[0] < ((NumericVector)r2[i])[0])
        return(0);      
        
  return(1);
  
}

//Input: a set of results, compared against r2 to check whether
// any element in r1 dominates r2
int checkWhetherDominates(vector<List> r1s, List r2) {
  
  for (int i = 0; i < r1s.size(); i++)
    if (dominates(r1s[i], r2) > 0 )
      return(1);
  
  return(0);
  
}

vector<int> giveDominatedMembers(vector<List> r1s, List r2) {
  
  vector<int> l_dominated;
  
  for (int i = 0; i < r1s.size(); i++) 
      l_dominated.push_back(dominates(r2, r1s[i]));

  
  return(l_dominated);
  
}

vector<int> findUnique(vector<int> k) {
  vector< int >::iterator r , w ;

  set< int > tmpset ;

	for( r = k.begin() , w = k.begin() ; r != k.end() ; ++r )
	{
		if( tmpset.insert( *r ).second )
		{
			*w++ = *r ;
		}
	}

	k.erase( w , k.end() );
  
  return(k);
}


bool myfunction (int i,int j) { return (i<j); }

List findNondominated(List population) {
  
  vector<int> indices;
  
  for (int i =0; i < population.size(); i++){
    List r1 = ((List)population[i])["result"];
    for (int j= 0; j < population.size(); j++) {
      if (i != j) { 
          List r2 = ((List)population[j])["result"];
          //if r1 dominates r2
          if (dominates(r1, r2) > 0)
           indices.push_back(j);
          
        }
      }
  }
  
  indices = findUnique(indices);
  sort (indices.begin(), indices.end(), myfunction);
  
  for(vector<int>::size_type i = indices.size() - 1; 
     i != (vector<int>::size_type) -1; i--)
        population.erase(population.begin()+indices[i]);

  /*  
  for (int i =0; i < population.size(); i++){
    
    List result = ((List)population[i])["result"];
    cout << endl;
    for (int j=0; j < result.size(); j++)
        cout << ((NumericVector)result[j])[0] << " ";
    
  }
  
  cout << endl<< "The nondominated pop size is: " << population.size() << endl;
  
  */
  return(population);
  
}
  

//[[Rcpp::export]]
SEXP collectNondominated(SEXP pop, SEXP arch) {
   
   List nonDominated_pop(pop);
   List archives(arch);

   
   List converging_archive = archives["converging"];
   List diversity_archive = archives["diversity"];
   
   nonDominated_pop = findNondominated(nonDominated_pop);

   for(int i =0; i < nonDominated_pop.size(); i++) {
     //if no member in both archives can dominate the indivisual(i)

     List result = ((List)nonDominated_pop[i])["result"];
     vector<List> ca_results = getResults(converging_archive);
     vector<List> da_results = getResults(diversity_archive);
     
    //enforce short-circuiting 
    if (checkWhetherDominates(ca_results, result) == 0) 
      if (checkWhetherDominates(da_results, result) == 0)  {
        //if individual(i) dominates any member in both archives then 
        int dFlag = 0;
        
        ca_results = getResults(converging_archive);
        da_results = getResults(diversity_archive);
     
        vector<int> ca_dominated = 
          giveDominatedMembers(ca_results, result);
        vector<int> da_dominated =
          giveDominatedMembers(da_results, result);
        if (any(ca_dominated)) {
           dFlag = 1; 
           //Delete elements
           for(vector<int>::size_type j = ca_dominated.size() - 1; 
              j != (vector<int>::size_type) -1; j--)
            if (ca_dominated[j] > 0)
                converging_archive.erase(converging_archive.begin()+j);
                
         }
         
        if (any(da_dominated)) {
           dFlag = 1;
           //Delete elements           
          for(vector<int>::size_type j = da_dominated.size() - 1; 
              j != (vector<int>::size_type) -1; j--)
            if (da_dominated[j] > 0) 
             diversity_archive.erase(diversity_archive.begin()+j);

         }

         if (dFlag > 0){
           char buffer [10];
           itoa (converging_archive.size(),buffer,10);
           converging_archive[string(buffer)] = nonDominated_pop[i];
           //cout << endl << "The size of converging archive: " <<
            //    converging_archive.size() << endl;
         
         }else{
           char buffer [10];
           itoa (diversity_archive.size(),buffer,10);
           diversity_archive[string(buffer)] = nonDominated_pop[i];
          // cout << endl << "The size of divesity archive: " <<
            //    diversity_archive.size() << endl;
         }
         
     
     }
     
     //Removing strategy to shrink the size of the archives below the limit
   // FIXME based on the euclidean distance
   //while(archives["converging"].size() + archives["converging"].size() > limit)
     //    || (length(archives$diverging) != 0)) {
     
     //index <- sample(1:length(archives$diverging), 1)    
     //archives$diverging <-  archives$diverging[-index]
     
   //}
     
   }
              
List x;
x["converging"] = converging_archive;
x["diversity"] = diversity_archive;

return(x);
//   return(List::create(Named("converging") = converging_archive,
//             Named("diversity") = diversity_archive));
}