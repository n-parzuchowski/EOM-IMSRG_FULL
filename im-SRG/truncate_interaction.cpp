#include <iostream> 
#include <string>
#include <fstream>

using namespace std;


int main( int argc, char* argv[] ) {

  string infile,outfile,Rstr;

  if (argc < 3 ) {
    cout << "Incorrect number of arguments.\n";
    return 0;
  }
  
  // read input parameters 
  if ( argc == 4 ) {
    infile = argv[1];
    outfile = argv[2];
    Rstr = argv[3]; 
  }
  else { 
    infile = argv[1];
    outfile = "temp.int";
    Rstr = argv[2];     
  }
  
  int R = stoi(Rstr); 
  int numstates = (R+1)*(R+2);     
  ifstream inputfile;
  ofstream outputfile;
  
  // open files check that they exist
  inputfile.open("../../TBME_input/"+infile);
  outputfile.open("../../TBME_input/"+outfile);
  
  if (! inputfile.is_open()) {
    cout << "File failed to open\n";
    return 0;
  }
  
  // read files and write
  string linein;
  int a,b,c,d;
  
 
  while ( getline (inputfile ,linein) ) {
 
    if (! isdigit(linein[3]) ) {
      
      outputfile << linein << "\n";
    }
    else { 
      a = stoi(linein.substr(12,4));
      if (a > numstates) continue;
      b = stoi(linein.substr(16,4));
      if (b > numstates) continue;
      c = stoi(linein.substr(20,4));
      if (c > numstates) continue;
      d = stoi(linein.substr(24,4));
      if (d > numstates) continue;
      
      outputfile << linein << "\n";
      
    }
  }
  
    inputfile.close();
    outputfile.close();
  
  return 0;
  
}
    
