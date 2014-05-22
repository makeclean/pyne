#include "material.h" 
#include "nucname.h"
#include "data.h"

int main (int argc, char* argv[])
{
   using namespace pyne;

   Material test_mat;

   // Define a comp_map by hand 
   comp_map nucvec;
   nucvec.insert(std::pair<int,double>(10010000, 1.0)); 
   nucvec.insert(std::pair<int,double>(80160000, 1.0));   
   nucvec.insert(std::pair<int,double>(691690000, 1.0));   
   nucvec.insert(std::pair<int,double>(922350000, 1.0));   
   nucvec.insert(std::pair<int,double>(922380000, 1.0));   
   nucvec.insert(std::pair<int,double>(942390000, 1.0));   
   nucvec.insert(std::pair<int,double>(942410000, 1.0));   
   nucvec.insert(std::pair<int,double>(952420000, 1.0));   
   nucvec.insert(std::pair<int,double>(962440000, 1.0));   

   // Create a basic filled Material object
   test_mat = Material(nucvec);

   // Test if can call added functionality in pyne library
   std::string mcnp_ret = test_mat.mcnp();
   std::cout << mcnp_ret << std::endl;
}