#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <iterator>
#include <map>

using namespace std;

int main() {
 
  string badChannels = "noisyChannels"; 
  string line;
  map<unsigned short, vector<unsigned char> > chmap;
  map<unsigned short, vector<unsigned char> >::iterator it;
  unsigned short a, b;
  unsigned char c;
  //istringstream strm;
  ifstream f;
  f.open(badChannels.c_str());
  if (f.is_open()) 
  {
    cout << "opened well" << endl;
    
    while (getline(f, line))
    {
      istringstream iss(line);
      vector<string> pieces;
 
      copy(istream_iterator<string>(iss), istream_iterator<string>(), back_inserter<vector<string> >(pieces));
      istringstream strm(pieces[5]); strm >> a;
      istringstream strm2(pieces[7]); strm2 >> b; c = (unsigned char)b;

      it = chmap.find(a);
      if (it != chmap.end())
      {
        it->second.push_back(b);
      } else {
        vector<unsigned char> vec;
        vec.push_back(b);
        chmap[a] = vec;
      }
    }
    f.close();
  
  string outputfile = "analysisMaskFileName.xml";
  ofstream f;
  f.open(outputfile.c_str());

  int arm, half, det, vfat;
  arm = half = det = vfat = -1;
  bool arm1, half1, det1;
  arm1 = half1 = det1 = true;
  
  if (f.is_open()) 

    f << "<top>" << endl << endl;    

    for (it = chmap.begin(); it != chmap.end(); it++)
    {
      if (((it->first - 2000*arm - 1000*half)/100 != det) && (det1 == false)) {
        f << "    </t2_detector>" << endl << endl;
      }

      if (((it->first - 2000*arm)/1000 != half) && (half1 == false)) {
        f << "  </t2_half>" << endl << endl;
      }
  
      if (((it->first)/2000 != arm) && (arm1 == false)) {
        f << "</arm>" << endl << endl;
      }
      
      if ((it->first)/2000 != arm) {
        arm1 = false;
        arm = (it->first)/2000;
        f << "<arm id=\"" << arm << "\">" << endl << endl;
      }
      
      if ((it->first - 2000*arm)/1000 != half) {
        half1 = false;
        half = (it->first - 2000*arm)/1000;
        f << "  <t2_half id=\"" << half << "\">" << endl << endl;
      }
      
      if ((it->first - 2000*arm - 1000*half)/100 != det) {
        det = (it->first - 2000*arm - 1000*half)/100;
        det1 = false;
        f << "    <t2_detector position=\"" << det << "\">" << endl << endl;
      }
      
      vfat = (it->first)%100;
      det = (it->first - 2000*arm - 1000*half)/100;
      f << "      <vfat iid=\"" << vfat << "\" fullmask=\"no\">" << endl;

      for (vector<unsigned char>::iterator ti = it->second.begin(); ti != it->second.end(); ti++)
      {
        f << "        <channel ch_num=\"" << (int)*ti << "\"/>" << endl;
      }
      f << "      </vfat>" << endl << endl;

    }
     
    f <<"    </t2_detector>"<<endl<<endl<<"  </t2_half>"<<endl<<endl<<"</arm>"<<endl<<endl<<"</top>"<<endl<<endl;
  } 
  else 
  {
    cout << "didn't open" << endl;
  }
  

}
