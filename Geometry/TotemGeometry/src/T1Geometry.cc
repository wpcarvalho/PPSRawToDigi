#include <fstream>
#include<stdlib.h>
#include <string>
#include <iostream>
#include <cassert>
#include "Geometry/TotemGeometry/interface/T1Geometry.h"

//#define _DEBUG_
using namespace std;

void T1Geometry::buildGeometry(std::string file) {
  std::cout <<" buildGeometry: Leggo il file. Costruisco la geometria ****************"<<std::endl;
  int plane=0;
  std::string type;
  float bmin=0.0;
  float bminOld=0.0;
  float bmax=0.0;
  float h=0.0;
  std::ifstream f;
  // std::string dummiest="/home/capurro/WD/CMSSW_1_0_1/src/Geometry/TotemGeometry/data/T1_data_geometry.dat";

  std::string cmssw_base = getenv("CMSSW_BASE");
  
  cmssw_base.append("/src/");
  cmssw_base.append(file);

  std::cout << " OPENING GEOMETRY FILE : "<< cmssw_base << std::endl;

  f.open(cmssw_base.c_str());
  //f.open(dummiest.c_str());
  std::string dummy;
  // if(!(f.is_open())) std::cout <<" ERROR OPEN FILE !!!!!!"<<std::endl;
  if(!(f.is_open()))
    {
      throw " ERROR OPENING GEOMETRY FILE !!!!!!";

    }
  if(f.is_open())
    {

      for (int i=0; i<8; i++){
	f >> dummy;
	//   std::cout << dummy <<std::endl;
      }
      for(int g=0; g<10 && !f.eof(); g++){


	f >> plane;
	f>> type;

	if(type =="P")
	  {
	    f>>bmin;
	    bminP[plane]=bmin;
	    f>>bminOld;
	    bminOldP[plane]=bminOld;
	    
	    f>>bmax;
	    
	    bmaxP[plane]=bmax;
	    f>>h;
	    _yP[plane]=h/2.;   
	    _xP[plane]=-(3./8.)*bminOld+(bmax/8.);
	    f>> widthP[plane];
	    offsetP[plane]=bmin;
	    f>>WireNP[plane];
	    //	  stripNP[plane]=(int)(widthP[plane]/parametri.Pitch());
	    f>>stripNP[plane];
	    // da verificare!!!


//	    offsetStripP[plane]=(widthP[plane]-parametri.Pitch()*(stripNP[plane]-1))/2.;
	    offsetStripP[plane]=0.5+2.25;
//	    std::cout << "offsetStripP["<<plane<<"]"<<offsetStripP[plane] << std::endl;
	    if(offsetStripP[plane]<0){
	      throw " ERROR : offsetStripP < 0 ! ";
	    }


	    // f>>offsetP[plane];
	  }
	if(type =="G")
	  {
	    f>>bmin;
	    bminG[plane]=bmin;
	    f>>bminOld;
	    bminOldG[plane]=bminOld;
	    f>>bmax;
	    bmaxG[plane]=bmax;
	    f>>h;
	    _yG[plane]=h/2.;   
	    _xG[plane]=-(3./8.)*bminOld+(bmax/8.);
	    f>> widthG[plane];
	    offsetG[plane]=bmin;
	    f>>WireNG[plane];
	    //	  stripNG[plane]=(int)(widthG[plane]/parametri.Pitch());
	    f>>stripNG[plane];
	    // da verificare!!!
	  //    offsetStripG[plane]=(widthG[plane]-parametri.Pitch()*(stripNG[plane]-1))/2.;
	    offsetStripG[plane]=0.5+2.25;
//	    std::cout << "offsetStripG["<<plane<<"]"<<offsetStripG[plane] << std::endl;
	    if(offsetStripG[plane]<0){
	      throw " ERROR : offsetStripG < 0 ! ";
	    }
	
	    //f>>stripNG[plane];
	    // f>>offsetG[plane];
	  }
	
      }
    }
  f.close();
  return; 
}

void T1Geometry::buildGeometry() {
  std::cout <<" buildGeometry: Leggo il file. Costruisco la geometria ****************"<<std::endl;

  int plane=0;
  std::string type;
  float bmin=0.0;
  float bminOld=0.0;
  float bmax=0.0;
  float h=0.0;
  std::ifstream f;
  //  std::string dummiest="/mnt/babar0_mnt/ferro/WD/CMSSW_1_2_3/src/Geometry/TotemGeometry/data/T1_data_geometry.dat";
  std::string 
    dummiest="/home/ferro/workdir2/CMSSW_1_2_3/src/Geometry/TotemGeometry/data/T1_data_geometry.dat";
  std::cout <<" buildGeometry: Leggo "<< dummiest <<" Costruisco la geometria ****************"<<std::endl;
  // f.open(file.c_str());
  f.open(dummiest.c_str());
  std::string dummy;
  if(!(f.is_open()))
    {
      throw " ERROR OPENING GEOMETRY FILE !!!!!!";

    }
  if(f.is_open())
    {

      for (int i=0; i<8; i++){
	f >> dummy;
	//   std::cout << dummy <<std::endl;
      }
      for(int g=0; g<10 && !f.eof(); g++){


	f >> plane;
	f>> type;

	if(type =="P")
	  {
	    f>>bmin;
	    bminP[plane]=bmin;
	    f>>bminOld;
	    bminOldP[plane]=bmin;
	    f>>bmax;
	    
	    bmaxP[plane]=bmax;
	    f>>h;
	    _yP[plane]=h/2.;   
	    _xP[plane]=-(3./8.)*bminOld+(bmax/8.);
	    f>> widthP[plane];
	    offsetP[plane]=bmin;
	    f>>WireNP[plane];
	    //	  stripNP[plane]=(int)(widthP[plane]/parametri.Pitch());
	    f>>stripNP[plane];
	    // da verificare!!!


//	    offsetStripP[plane]=(widthP[plane]-parametri.Pitch()*(stripNP[plane]-1))/2.;
offsetStripP[plane]=0.5+2.25;
// std::cout << "offsetStripP["<<plane<<"]"<<offsetStripP[plane] << std::endl;

	    if(offsetStripP[plane]<0){
	      throw " ERROR : offsetStripP < 0 ! ";
	    }


	    // f>>offsetP[plane];
	  }
	if(type =="G")
	  {
	    f>>bmin;
	    bminG[plane]=bmin;
	    f>>bminOld;
	    bminOldG[plane]=bminOld;
	    f>>bmax;
	    bmaxG[plane]=bmax;
	    f>>h;
	    _yG[plane]=h/2.;   
	    _xG[plane]=-(3./8.)*bminOld+(bmax/8.);
	    f>> widthG[plane];
	    offsetG[plane]=bmin;
	    f>>WireNG[plane];
	    //	  stripNG[plane]=(int)(widthG[plane]/parametri.Pitch());
	    f>>stripNG[plane];
	    // da verificare!!!
//	    offsetStripG[plane]=(widthG[plane]-parametri.Pitch()*(stripNG[plane]-1))/2.;
	    
offsetStripG[plane]=0.5+2.25;
//std::cout << "offsetStripG["<<plane<<"]"<<offsetStripG[plane] << std::endl;
	    if(offsetStripG[plane]<0){
	      throw " ERROR : offsetStripG < 0 ! ";
	    }
	
	    //f>>stripNG[plane];
	    // f>>offsetG[plane];
	  }
	
      }
    }



  f.close();
  return; 
}

float T1Geometry::getGlobalZ(uint32_t detid){
    T1DetId object(detid);
    float arm_factor = 1;
    if(object.Arm()==1)arm_factor=-1;
    float globalZ = Z[object.Plane()][object.CSC()];
    assert(globalZ>7000 && globalZ<11000);
    return arm_factor*globalZ;
      
}

void T1Geometry::buildAlignment(){


  cout << "Reading T1 alignment parameters (if available)" << endl;
  cout << endl;
  cout << " ***** Alignment applied ONLY if T1Geometry::Align method called ***** " << endl;
  cout << endl;


  for(int i=0; i<2; i++){
    Arm_Align_x[i]=0;
    Arm_Align_y[i]=0;
    Arm_Align_z[i]=0;
    for(int j=0; j<5; j++)
      for(int k=0; k<6; k++){
	Align_x[i][j][k]=0;
	Align_y[i][j][k]=0;
	Align_th[i][j][k]=0;
	Sex_Align_x[i][k]=0;
	Sex_Align_y[i][k]=0;
      }
  }





  unsigned int a,p,c;
  float x_,y_,z_,th_;

  bool theAlignment = true;
  bool theSexAlignment = true;
  bool theArmAlignment = true;
  string theAlignmentFile = "TotemAlignment/T1Alignment/data/T1Align.dat";
  string theSexAlignmentFile = "TotemAlignment/T1Alignment/data/T1SexAlign.dat";
  string theArmAlignmentFile = "TotemAlignment/T1Alignment/data/T1ArmAlign.dat";
  FILE * file;

  if(theAlignment){
    
    char *cmsswPath = getenv("CMSSW_BASE");
    if (cmsswPath && theAlignmentFile[0] != '/')      
      theAlignmentFile= string(cmsswPath) + string("/src/") + theAlignmentFile;
    file = fopen(theAlignmentFile.c_str(),"r");
    if(file!=NULL){
    while(!feof(file)){
  


      fscanf(file,"%u%u%u%f%f%f",&a,&p,&c,&x_,&y_,&th_);


      Align_x[a][p][c]=x_;
      Align_y[a][p][c]=y_;
      Align_th[a][p][c]=th_;
   }
     
    fclose(file);
    cout << " Alignment parameters " <<endl;
    for(int i=0; i<2; i++)
    for(int j=0; j<5; j++)
      for(int k=0; k<6; k++){
	cout <<i << " "<< j << " "<< k << " "<<  Align_x[i][j][k]<< " "
	     << Align_y[i][j][k] << " "<< Align_th[i][j][k] << endl;
      }
    }else{
      cout << "WARNING: Alignment file " << theAlignmentFile << " not found. Running with no external alignment " << endl;
    }


  }else{
    cout << " No external alignment " << endl;
  }


  if(theSexAlignment){
    char *cmsswPath = getenv("CMSSW_BASE");
    if (cmsswPath && theSexAlignmentFile[0] != '/')      
      theSexAlignmentFile= string(cmsswPath) + string("/src/") + theSexAlignmentFile;
    file = fopen(theSexAlignmentFile.c_str(),"r");
    if(file!=NULL){
    while(!feof(file)){
  


      fscanf(file,"%u%u%f%f",&a,&c,&x_,&y_);


      Sex_Align_x[a][c]=x_;
      Sex_Align_y[a][c]=y_;
  
   }
     
    fclose(file);
    cout << " Sextant Alignment parameters " <<endl;
    for(int i=0; i<2; i++)
      for(int k=0; k<6; k++){
	cout <<i << " "<< k << " "<<  Sex_Align_x[i][k]<< " "
	     << Sex_Align_y[i][k] <<  endl;
      }
    }else{
      cout << "WARNING: Sextant Alignment file " << theSexAlignmentFile << " not found. Running with no external SEXTANT alignment " << endl;
    }
    
    
  }else{
    cout << " No external SEXTANT alignment " << endl;
  }
  
  






  if(theArmAlignment){
    char *cmsswPath = getenv("CMSSW_BASE");
    if (cmsswPath && theArmAlignmentFile[0] != '/')      
      theArmAlignmentFile= string(cmsswPath) + string("/src/") + theArmAlignmentFile;
    file = fopen(theArmAlignmentFile.c_str(),"r");
    if(file!=NULL){
  
  


      fscanf(file,"%f%f%f",&x_,&y_,&z_);


      Arm_Align_x[0]=x_;
      Arm_Align_y[0]=y_;
      Arm_Align_z[0]=z_;

      fscanf(file,"%f%f%f",&x_,&y_,&z_);


      Arm_Align_x[1]=x_;
      Arm_Align_y[1]=y_;
      Arm_Align_z[1]=z_;

  
   
     
    fclose(file);
    cout << " Arm Alignment parameters " <<endl;
    for(int i=0; i<2; i++)
      {
	cout <<i << " "<<  Arm_Align_x[i] << " "
	     << Arm_Align_y[i] << " "
	     << Arm_Align_z[i] <<endl;
      }
    }else{
      cout << "WARNING: Arm Alignment file " << theArmAlignmentFile << " not found. Running with no global arm alignment " << endl;
    }
    
    
  }else{
    cout << " No external arm alignment " << endl;
  }
  
  








}
