#ifndef _T1_GEO_GEO_
#define _T1_GEO_GEO_

#include <iostream>
#include <string>
#include "assert.h"
#include "DataFormats/T1DetId/interface/T1DetId.h"
#include "Geometry/TotemGeometry/interface/T1ChamberSpecs.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"

#ifndef RADQ3
#define RADQ3 1.7320508
#endif
#ifndef PI
#define PI 3.14159265
#endif


//#define _DEBUG_

class T1Geometry {

  public:
   T1Geometry(){

//    float Rmin_temp[5]={148,158,167,178,187.5};

  float Rmin_temp[5]={144.5, 155.5, 164.5, 176.5, 185.5}; 
/*     float Z_temp[5][6]={{7521.4,7569.6,7526.4,7569.6,7521.4,7564.6},
                         {8191.4,8239.4,8196.4,8239.6,8191.4,8234.6},
                         {8798.4,8846.6,8803.4,8846.6,8798.4,8841.6},
                         {9405.4,9453.6,9410.4,9453.6,9405.4,9448.6},
                         {10168.4,10216.6,10173.4,10216.6,10168.4,10211.6}};
*/
float Z_temp[5][6]={{ 7569.6, 7521.4, 7564.6, 7521.4, 7569.6,  7526.4,},
                     { 8239.4, 8191.4, 8234.6, 8191.4, 8239.6,  8196.4,},
                     { 8846.6, 8798.4, 8841.6, 8798.4, 8846.6,  8803.4,},
                     { 9453.6, 9405.4, 9448.6, 9405.4, 9453.6,  9410.4,},
                     {10216.6,10168.4,10211.6,10168.4,10216.6, 10173.4,}};


     float ROT0_temp[5]={66,63,54,57,60};
     float disassamentoX_temp[6]={3.035,-3.035,-3.035,3.035,-3.035,3.035};

     for(int i=0; i<5; i++){
       _xP[i]=0;
       _xG[i]=0;
       _yP[i]=0;
       _yG[i]=0;
       widthP[i]=0;
       offsetStripG[i]=0;
       offsetStripP[i]=0;
       stripNG[i]=0;
       widthG[i]=0;
       WireNP[i]=0;
       WireNG[i]=0;
       stripNP[i]=0;
       bmaxP[i]=0;
       bminP[i]=0;
       bminOldP[i]=0;
       bmaxG[i]=0;
       bminG[i]=0;
       bminOldG[i]=0;
       offsetP[i]=0;
       offsetG[i]=0;

       RminSens[i]=Rmin_temp[i];
       Rot0[i]=ROT0_temp[i];

       for(int g=0;g<6;g++){
         Z[i][g]=Z_temp[i][g];
         if(i==0)disassamentoX[g] = disassamentoX_temp[g];
       }
     }

     


     std::cout << " Calling buildGeometry " <<std::endl;
     buildGeometry();
     std::cout << " Calling buildAlignment " <<std::endl;
     buildAlignment();

   }
   T1Geometry(std::string filename){
     //   float Rmin_temp[5]={148,158,167,178,187.5};
  float Rmin_temp[5]={144.5, 155.5, 164.5, 176.5, 185.5}; 
/*     float Z_temp[5][6]={{7521.4,7569.6,7526.4,7569.6,7521.4,7564.6},
                         {8191.4,8239.4,8196.4,8239.6,8191.4,8234.6},
                         {8798.4,8846.6,8803.4,8846.6,8798.4,8841.6},
                         {9405.4,9453.6,9410.4,9453.6,9405.4,9448.6},
                         {10168.4,10216.6,10173.4,10216.6,10168.4,10211.6}};
*/

float Z_temp[5][6]={{ 7569.6, 7521.4, 7564.6, 7521.4, 7569.6,  7526.4,},
                     { 8239.4, 8191.4, 8234.6, 8191.4, 8239.6,  8196.4,},
                     { 8846.6, 8798.4, 8841.6, 8798.4, 8846.6,  8803.4,},
                     { 9453.6, 9405.4, 9448.6, 9405.4, 9453.6,  9410.4,},
                     {10216.6,10168.4,10211.6,10168.4,10216.6, 10173.4,}};

     float ROT0_temp[5]={66,63,54,57,60};
     float disassamentoX_temp[6]={3.035,-3.035,-3.035,3.035,-3.035,3.035};
     for(int i=0; i<5; i++){
       _xP[i]=0;
       _xG[i]=0;
       _yP[i]=0;
       _yG[i]=0;
       widthP[i]=0;
       offsetStripG[i]=0;
       offsetStripP[i]=0;
       stripNG[i]=0;
       widthG[i]=0;
       WireNP[i]=0;
       WireNG[i]=0;
       stripNP[i]=0;
       bmaxP[i]=0;
       bminP[i]=0;
       bminOldP[i]=0;
       bmaxG[i]=0;
       bminG[i]=0;
       bminOldG[i]=0;
       offsetP[i]=0;
       offsetG[i]=0;

       RminSens[i]=Rmin_temp[i];
       Rot0[i]=ROT0_temp[i];

       for(int g=0;g<6;g++){
         Z[i][g]=Z_temp[i][g];
         if(i==0)disassamentoX[g] = disassamentoX_temp[g];
}
     }

     std::cout << " Calling buildGeometry " <<std::endl;
     buildGeometry(filename);
     std::cout << " Calling buildAlignment " <<std::endl;
     buildAlignment();

   }

   ~T1Geometry(){}

   void buildGeometry(std::string );
   void buildGeometry();
   void buildAlignment();
   float getGlobalZ(uint32_t);
   float getX(int plane, std::string  type  )  {
   
     if( (type != "P" && type != "G") || plane < 0 || plane > 4){
       std::cout <<"Error: Ivalid Plane or CSC type " <<std::endl;
       std::cout <<"Plane: 0 to 4   CSC type: P or G"<<std::endl;
     }
     if(type=="P")
       return _xP[plane];
     if(type=="G")
       return _xG[plane];    
   }
   float getY(int plane, std::string type)  {

     if((type != "P" && type != "G") || plane < 0 || plane > 4){
       std::cout <<"Error: Ivalid Plane or CSC type " <<std::endl;
       std::cout <<"Plane: 0 to 4   CSC type: P or G"<<std::endl;
     }
     if(type=="P")
       return _yP[plane];
     if(type=="G")
       return _yG[plane];    
   }
   //metodo che restituisce la y del filo nel sistema di riferimento di Geant4 (xml)
   float yOfWire(uint32_t detid, int wire )  {
     T1DetId object(detid);
     float yWire=0.0;
#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif 
     if(object.CSC() == 5 || object.CSC() == 2){

       //    yWire = wire*parametri.wireSpacing()-_yP[object.Plane()] - parametri.wireSpacing()/2.;
       yWire = wire*parametri.wireSpacing()-_yP[object.Plane()];
       return yWire;
     }
     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){     
//     yWire = wire*parametri.wireSpacing()-_yG[object.Plane()] - parametri.wireSpacing()/2.;
       yWire = wire*parametri.wireSpacing()-_yG[object.Plane()];
       return yWire;
     }
     return 0.0;
   }

   //y del filo nel sistema di riferimento centrato nell'angolo in basso a sinistra della camera signal side
   float yOfWireCornerRS(uint32_t detid, int wire)  {
     T1DetId object(detid);
     float yWire=0.0;
#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif 
     if(object.CSC() == 5 || object.CSC() == 2){

       yWire = wire*parametri.wireSpacing() ;
    
       return yWire;
     }
     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){          yWire = wire*parametri.wireSpacing() ;

     return yWire;
     }

   }

   //x of center of wire
   float xOfWire(uint32_t detid){//, int wire)  {
     T1DetId object(detid);
     float xWire=0.0;
#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif 
     if(object.CSC() == 5 || object.CSC() == 2){

       xWire = -(bminOldP[object.Plane()]/2.-_xP[object.Plane()]);
    
       return xWire;
     }
     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){ 
       xWire = -(bminOldG[object.Plane()]/2.-_xG[object.Plane()]);

       return xWire;
     }
     return 0.0;
   }
   //Geant4 reference system
   float xOfStripA(uint32_t detid, int strip, float ywire)  {
     T1DetId object(detid);
     float xStrip=0.0;
     float ywirecorner=0.;
     float xstripcorner=0.;

#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif 
     if(object.CSC() == 5 || object.CSC() == 2){

       ywirecorner = ywire + _yP[object.Plane()];
       //offsetStripP al momento e nullo 
       xstripcorner=(ywirecorner-2*(offsetStripP[object.Plane()]+((float)strip-1)*parametri.Pitch()))/RADQ3;
       xStrip=xstripcorner-_xP[object.Plane()];
       return xStrip;
     }
     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){ 
      
       ywirecorner = ywire + _yG[object.Plane()];
       //offsetStripP al momento e nullo 

       xstripcorner=(ywirecorner-2*(offsetStripG[object.Plane()]+((float)strip-1)*parametri.Pitch()))/RADQ3;


       xStrip=xstripcorner-_xG[object.Plane()];
       return xStrip;
     }
     return -1;
   }



   float xOfStripBWireCrossing(uint32_t detid, float strip, int wire)  {
     T1DetId object(detid);
     float xStrip=0.0;
     float ywirecorner=0.;
     float xstripcorner=0.;

     ywirecorner = ( (float)wire )* parametri.wireSpacing();

     
#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif 
     if(object.CSC() == 5 || object.CSC() == 2){

       //      ywirecorner = ywire + _yP[object.Plane()];
       //offsetStripP al momento e nullo 
       xstripcorner=-ywirecorner/RADQ3 + 2.*(offsetStripP[object.Plane()]+((float)strip-1)*parametri.Pitch())/RADQ3 - offsetP[object.Plane()];
       xStrip=xstripcorner-_xP[object.Plane()];
       return xStrip;
     }
     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){ 
      
       // ywirecorner = ywire + _yG[object.Plane()];
       //offsetStripP al momento e nullo 
       xstripcorner=-ywirecorner/RADQ3 + 2.*(offsetStripG[object.Plane()]+((float)strip-1)*parametri.Pitch())/RADQ3 - offsetG[object.Plane()];
       xStrip=xstripcorner-_xG[object.Plane()];
       return xStrip;
     }
     return -1;
   }



   float xOfStripAWireCrossing(uint32_t detid, float strip, int wire)  {
     T1DetId object(detid);
     float xStrip=0.0;
     float ywirecorner=0.;
     float xstripcorner=0.;

     ywirecorner = ( (float)wire )* parametri.wireSpacing();



#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif 
     if(object.CSC() == 5 || object.CSC() == 2){

       //       ywirecorner = ywire + _yP[object.Plane()];
       //offsetStripP al momento e nullo 
       xstripcorner=(ywirecorner-2*(offsetStripP[object.Plane()]+((float)strip-1)*parametri.Pitch()))/RADQ3;
       xStrip=xstripcorner-_xP[object.Plane()];
       return xStrip;
     }
     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){ 
      
       //  ywirecorner = ywire + _yG[object.Plane()];
       //offsetStripP al momento e nullo 

       xstripcorner=(ywirecorner-2*(offsetStripG[object.Plane()]+((float)strip-1)*parametri.Pitch()))/RADQ3;


       xStrip=xstripcorner-_xG[object.Plane()];
       return xStrip;
     }
     return -1;
   }






   //Geant4 reference system
   float xOfStripB(uint32_t detid, int strip, float ywire)  {
     T1DetId object(detid);
     float xStrip=0.0;
     float ywirecorner=0.;
     float xstripcorner=0.;

#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif 
     if(object.CSC() == 5 || object.CSC() == 2){

       ywirecorner = ywire + _yP[object.Plane()];
       //offsetStripP al momento e nullo 
       xstripcorner=-ywirecorner/RADQ3 + 2.*(offsetStripP[object.Plane()]+((float)strip-1)*parametri.Pitch())/RADQ3 - offsetP[object.Plane()];
       xStrip=xstripcorner-_xP[object.Plane()];
       return xStrip;
     }
     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){ 
      
       ywirecorner = ywire + _yG[object.Plane()];
       //offsetStripP al momento e nullo 
       xstripcorner=-ywirecorner/RADQ3 + 2.*(offsetStripG[object.Plane()]+((float)strip-1)*parametri.Pitch())/RADQ3 - offsetG[object.Plane()];
       xStrip=xstripcorner-_xG[object.Plane()];
       return xStrip;
     }
     return -1;
   }


   float yOfStripCrossing(uint32_t detid, float stripA, float stripB){
     T1DetId object(detid);
     float ystripcorner=0.0;
     float yStrip=0.0;
     if(object.CSC() == 5 || object.CSC() == 2){
       ystripcorner=((float)stripB+(float)stripA-2.0)*parametri.Pitch() - offsetP[object.Plane()]*RADQ3/2. + 2*offsetStripP[object.Plane()];
       yStrip=ystripcorner-_yP[object.Plane()];
       return yStrip;
     }

     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){ 

       //      std::cout << stripA << " " <<stripB  << " " <<offsetG[object.Plane()] << " " <<offsetStripG[object.Plane()] << " " <<_yG[object.Plane()] << " " <<std::endl;
       ystripcorner=((float)stripB+(float)stripA-2.0)*parametri.Pitch() - offsetG[object.Plane()]*RADQ3/2. + 2*offsetStripG[object.Plane()];
       yStrip=ystripcorner-_yG[object.Plane()];
       return yStrip;
     }
     return 0.0;
   }

   float xOfStripCrossing(uint32_t detid, float stripA, float stripB){
     T1DetId object(detid);
     float xstripcorner=0.;
     float xStrip=0.0;
     if(object.CSC() == 5 || object.CSC() == 2){
       xstripcorner=((float)stripB-(float)stripA)*parametri.Pitch()/RADQ3 - offsetP[object.Plane()]/2.;
       xStrip=xstripcorner-_xP[object.Plane()];
       return xStrip;
     }

     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){ 
       xstripcorner=((float)stripB-(float)stripA)*parametri.Pitch()/RADQ3 - offsetG[object.Plane()]/2.;
       xStrip=xstripcorner-_xG[object.Plane()];
       return xStrip;
     }
     return 0.0;
   }


   // from G4 system to a system centered in the beam axis and with y axis going through the middle of the chamber
   float xFromLocal2BeamSystem(uint32_t detid, float localx){
     T1DetId object(detid);
     float xstripcorner=0.;
     float x =0.0;
     if(object.CSC() == 5 || object.CSC() == 2){
       xstripcorner=localx+_xP[object.Plane()];
       x = xstripcorner + bminOldP[object.Plane()]/2.-disassamentoX[object.CSC()] ;
       return x;
     }

     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){ 
       xstripcorner=localx+_xG[object.Plane()];
       x = xstripcorner + bminOldG[object.Plane()]/2.-disassamentoX[object.CSC()];
       return x;
     }

     return 0.0;
   }

   // from G4 system to a system centered in the beam axis and with y axis going through the middle of the chamber
   float yFromLocal2BeamSystem(uint32_t detid,float localy){
     T1DetId object(detid);
  
     float y =0.0;
     if(object.CSC() == 5 || object.CSC() == 2){
       y=localy+_yP[object.Plane()]+RminSens[object.Plane()];
    
       return y;
     }

     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){ 
       y=localy+_yG[object.Plane()]+RminSens[object.Plane()];
   
       return y;
     }

     return 0.0;
   }


void Local2Global(uint32_t detid,float lx, float ly, float & x, float & y, float & z){

     T1DetId object(detid);

     float globalX, globalY;

     float tempX = xFromLocal2BeamSystem(detid,lx);
     float tempY = yFromLocal2BeamSystem(detid,ly);
     RotationLocal2Global(detid,tempX,tempY,globalX,globalY);
     
     object.Arm();
     object.Plane();
     object.CSC();
     
     x=globalX;
     y=globalY;
     z=Zeta(detid);

}

void Align(uint32_t detid,float globalX,float globalY,float & x,float & y,float & z){
   T1DetId object(detid);
     int lArm=object.Arm();
     int ilPlane=object.Plane();
     int laCSC=object.CSC();
     
   x = globalX*cos( Align_th[lArm][ilPlane][laCSC] ) - globalY*sin( Align_th[lArm][ilPlane][laCSC] ) + Align_x[lArm][ilPlane][laCSC] + Sex_Align_x[lArm][laCSC] + Arm_Align_x[lArm];
   y = globalX*sin( Align_th[lArm][ilPlane][laCSC] ) + globalY*cos( Align_th[lArm][ilPlane][laCSC] ) + Align_y[lArm][ilPlane][laCSC] + Sex_Align_y[lArm][laCSC] + Arm_Align_y[lArm];
   z = Zeta(detid) + Arm_Align_z[lArm];
}
void InvAlign(uint32_t detid,float globalX,float globalY,float & x,float & y,float & z){
   T1DetId object(detid);
     int lArm=object.Arm();
     int ilPlane=object.Plane();
     int laCSC=object.CSC();
  
     x = (globalX - Align_x[lArm][ilPlane][laCSC] - Sex_Align_x[lArm][laCSC] - Arm_Align_x[lArm])*cos( Align_th[lArm][ilPlane][laCSC] )  + (globalY - Align_y[lArm][ilPlane][laCSC] - Sex_Align_y[lArm][laCSC] - Arm_Align_y[lArm])*sin( Align_th[lArm][ilPlane][laCSC] );
     y = -((globalX - Align_x[lArm][ilPlane][laCSC] - Sex_Align_x[lArm][laCSC] - Arm_Align_x[lArm])*sin( Align_th[lArm][ilPlane][laCSC] ) )  + (globalY - Align_y[lArm][ilPlane][laCSC] - Sex_Align_y[lArm][laCSC] - Arm_Align_y[lArm])*cos( Align_th[lArm][ilPlane][laCSC] );
     z = Zeta(detid);

}




   void RotationLocal2Global(uint32_t detid,float lx, float ly, float & x, float & y){
     T1DetId object(detid);

     float phi = 0;
     phi = (Rot0[object.Plane()]+object.CSC()*60-90)*(PI/180.);
     // ATTENZIONE all' ARM negativo

     //   std::cout <<Rot0[object.Plane()]<< " " << object.CSC()<<" " << std::endl;
     //   std::cout << " PHI " << phi << std::endl;

     x = lx*cos(phi)-ly*sin(phi);
     if(object.Arm()==1)
       x=-x;

     y = lx*sin(phi)+ly*cos(phi);
     return;
   }


   void RotationGlobal2Local(uint32_t detid,float gx, float gy, float & x, float & y){
     T1DetId object(detid);
     if(object.Arm()==1)
       gx=-gx;
     float phi = 0;
     phi = -(Rot0[object.Plane()]+object.CSC()*60-90)*(PI/180.);
     // ATTENZIONE all' ARM negativo

     //   std::cout <<Rot0[object.Plane()]<< " " << object.CSC()<<" " << std::endl;
     //   std::cout << " PHI " << phi << std::endl;

     x = gx*cos(phi)-gy*sin(phi);
    
     y = gx*sin(phi)+gy*cos(phi);
     return;
   }

   float xFromBeamSystem2Local(uint32_t detid, float globalx){
     T1DetId object(detid);
     float xstripcorner=0.;
     float x =0.0;
     if(object.CSC() == 5 || object.CSC() == 2){
      xstripcorner = globalx- bminOldP[object.Plane()]/2.+disassamentoX[object.CSC()];
       x =  xstripcorner - _xP[object.Plane()];
       return x;
     }

     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){ 



       xstripcorner = globalx- bminOldG[object.Plane()]/2.+disassamentoX[object.CSC()];
       x =  xstripcorner - _xG[object.Plane()];
    

       return x;
     }
     assert(object.CSC()>=0 && object.CSC()<=5);
     return 0.0;

   }

   float yFromBeamSystem2Local(uint32_t detid,float globaly){
     T1DetId object(detid);
  
     float y =0.0;
     if(object.CSC() == 5 || object.CSC() == 2){
       y=globaly-_yP[object.Plane()]-RminSens[object.Plane()];
    
       return y;
     }

     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){ 
         y=globaly-_yG[object.Plane()]-RminSens[object.Plane()];
   
       return y;
     }

     assert(object.CSC()>=0 && object.CSC()<=5);
     return 0.0;
   }


   void RotationLocalError2Global(uint32_t detid,float lex, float ley, float & ex, float & ey){
     T1DetId object(detid);

     float phi = 0;
     phi = (Rot0[object.Plane()]+object.CSC()*60-90)*(PI/180.);
     // ATTENZIONE all' ARM negativo

    
/*
     ex = sqrt(lex*lex*fabs(cos(phi))*fabs(cos(phi))+ley*fabs(sin(phi))*ley*fabs(sin(phi)));
     

     ey = sqrt(lex*fabs(sin(phi))*lex*fabs(sin(phi))+ley*fabs(cos(phi))*ley*fabs(cos(phi)));
 */
     ex = lex*cos(phi)*cos(phi)+sin(phi)*ley*sin(phi);
     

     ey = lex*sin(phi)*sin(phi)+cos(phi)*ley*cos(phi);
     return;
   }



   //x of center of wire in the bottom left corner ref. sys. (bmin/2)
   float xOfWireCornerRS(uint32_t detid, int wire)  {
     T1DetId object(detid);
     float xWire=0.0;
#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif 
     if(object.CSC() == 5 || object.CSC() == 2){

       xWire = bminP[object.Plane()]/2.;
    
       return xWire;
     }
     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){ 
       xWire = bminG[object.Plane()]/2.;

       return xWire;
     }

   }

   //metodo che restituisce il numero di Wire per camera
   int numberOfwires(uint32_t detid)  {

     T1DetId object(detid);
#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif 
     if(object.CSC() == 5 || object.CSC() == 2)
       {
         return WireNP[object.Plane()];
       }
     
     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 )
       {
         return WireNG[object.Plane()];
       }
   }


   //metodo che restituisce il numero di Strips per camera
   int numberOfStrips(uint32_t detid)  {

     T1DetId object(detid);
#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif 
     if(object.CSC() == 5 || object.CSC() == 2)
       {
         return stripNP[object.Plane()];
       }
     
     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 )
       {
         return stripNG[object.Plane()];
       }
     return -1;
   }


   //metodo che restituisce il numero di Strips per camera
   int numberOfStrips(T1DetId object)  {

     // T1DetId object(detid);
#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif 
     if(object.CSC() == 5 || object.CSC() == 2)
       {
         return stripNP[object.Plane()];
       }
     
     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 )
       {
         return stripNG[object.Plane()];
       }
     return -1;
   }

   float Zeta(uint32_t detid){
     T1DetId object(detid);
     if(object.Arm()==1) return -Z[object.Plane()][object.CSC()];
     return Z[object.Plane()][object.CSC()];
   }

   float eta(uint32_t detid, float localx, float localy){
     T1DetId object(detid);
     float x=0;
     float y=0;
     float ETA=0;
     
#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif 
     if(object.CSC() == 5 || object.CSC() == 2)
       {
         //using reference system centered on the bottom left corner of the chamber SIGNAL side
         x = localx+_xP[object.Plane()]+bminOldP[object.Plane()]/2.;
         y = localy+_yP[object.Plane()];
         y = y + RminSens[object.Plane()];

         float c = atan(sqrt(x*x+y*y)/Z[object.Plane()][object.CSC()]);

         ETA = -log(tan(c/2.));

         if(object.Arm()==1)ETA = -ETA;

       } else{
         //using reference system centered on the bottom left corner of the chamber SIGNAL side
         x = localx+_xG[object.Plane()]+bminOldG[object.Plane()]/2.;
         y = localy+_yG[object.Plane()];
         y = y + RminSens[object.Plane()];

         float c = atan(sqrt(x*x+y*y)/Z[object.Plane()][object.CSC()]);

         ETA = -log(tan(c/2.));

         if(object.Arm()==1)ETA = -ETA;


       }


     return ETA;

   }



   
   bool inside(uint32_t detid, LocalPoint newPoint)  {
     T1DetId object(detid);
#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif 
     double x = newPoint.x();
     double y = newPoint.y();
     double z = newPoint.z();
#ifdef _DEBUG_
     std::cout << " x " <<  x << " y " << y<<" z " << z<<std::endl;
#endif 

     double y0=0.0;
     double y1=0.0;
     double h=0.0;

     if(object.CSC() == 5 || object.CSC() == 2)
       {
         //using reference system centered on the bottom left corner of the chamber SIGNAL side
         x = x+_xP[object.Plane()];
         y = y+_yP[object.Plane()];

         h=WireNP[object.Plane()]*parametri.wireSpacing();
         y0=((bmaxP[object.Plane()]-bminOldP[object.Plane()])/2.)*RADQ3 + h;
         y1=((bmaxP[object.Plane()]+bminOldP[object.Plane()])/2.)*RADQ3 +h;
         if( z>-5 && z<5 &&  y<h && y>0 && y > (x*RADQ3) && y < -RADQ3*x + y0 && y < RADQ3*x + y1 && y>-RADQ3*(x+bminOldP[object.Plane()])) {return true;}
         else {
           //     std::cout << " ERROR " << x << " " << y << " " << z << std::endl; 
           return false;
         } 
   
       }
     
     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 )
       {
         x = x+_xG[object.Plane()];
         y = y+_yG[object.Plane()];
         h=WireNG[object.Plane()]*parametri.wireSpacing();
         y0=((bmaxG[object.Plane()]-bminOldG[object.Plane()])/2.)*RADQ3 + h;
         y1=((bmaxG[object.Plane()]+bminOldG[object.Plane()])/2.)*RADQ3 +h;
         if( z>-5 && z<5 &&  y<h && y>0 && y > (x*RADQ3) && y < -RADQ3*x + y0 && y < RADQ3*x + y1 && y>-RADQ3*(x+bminOldG[object.Plane()])) {return true;}
         else {
           //     std::cout << x << " " << y << " "<<z<<std::endl;
           //     std::cout << h << " " <<  (x*RADQ3)  << " " << -RADQ3*x + y0 << " " << RADQ3*x + y1 << " " <<-RADQ3*(x+bminOldG[object.Plane()])<<std::endl;

           return false;} 
       }
     return false;
   }
   //////////
   bool insideRECO(T1DetId object, LocalPoint newPoint)  {


     //T1DetId object(detid);
#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif 
     double x = newPoint.x();
     double y = newPoint.y();
     double z = newPoint.z();
#ifdef _DEBUG_
     std::cout << " x " <<  x << " y " << y<<" z " << z<<std::endl;
#endif 

     double y0=0.0;
     double y1=0.0;
     double h=0.0;

     if(object.CSC() == 5 || object.CSC() == 2)
       {
         //using reference system centered on the bottom left corner of the chamber SIGNAL side
         x = x+_xP[object.Plane()];
         y = y+_yP[object.Plane()];

         float DELTA = fabs((bminOldP[object.Plane()]-bminP[object.Plane()])/2.);
//      std::cout << " DELTA " <<  DELTA << std::endl;

         h=WireNP[object.Plane()]*parametri.wireSpacing();
         y0=((bmaxP[object.Plane()]-bminOldP[object.Plane()])/2.)*RADQ3 + h;
         y1=((bmaxP[object.Plane()]+bminOldP[object.Plane()])/2.)*RADQ3 +h;
         if( z>-5 && z<5 &&  y<h && y>0 && y > ((x+DELTA)*RADQ3) && y < -RADQ3*x + y0 && y < RADQ3*x + y1 && y>-RADQ3*(x+bminOldP[object.Plane()]-DELTA)) {return true;}
         else {
           //  std::cout << " ERROR in Plane "<< object.Plane()<< " CSC " << object.CSC()<<"  (" << x << " " << y << " " << z <<")"<< std::endl; 
           return false;
         } 
   
       }
     
     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 )
       {
         x = x+_xG[object.Plane()];
         y = y+_yG[object.Plane()];
         float DELTA = fabs((bminOldG[object.Plane()]-bminG[object.Plane()])/2.);//      std::cout << " DELTA " <<  DELTA << std::endl;
         h=WireNG[object.Plane()]*parametri.wireSpacing();
         y0=((bmaxG[object.Plane()]-bminOldG[object.Plane()])/2.)*RADQ3 + h;
         y1=((bmaxG[object.Plane()]+bminOldG[object.Plane()])/2.)*RADQ3 +h;
         if( z>-5 && z<5 &&  y<h && y>0 && y > ((x+DELTA)*RADQ3) && y < -RADQ3*x + y0 && y < RADQ3*x + y1 && y>-RADQ3*(x+bminOldG[object.Plane()]-DELTA)) {return true;}
         else {
           //     std::cout << x << " " << y << " "<<z<<std::endl;
           //     std::cout << h << " " <<  (x*RADQ3)  << " " << -RADQ3*x + y0 << " " << RADQ3*x + y1 << " " <<-RADQ3*(x+bminOldG[object.Plane()])<<std::endl;
           // std::cout << " ERROR in Plane "<< object.Plane()<< " CSC " << object.CSC()<<"  (" << x << " " << y << " " << z <<")"<< std::endl; 
           return false;} 
       }

     std::cout << " SOMETHING WRONG IN insideReco! " << std::endl;
     return false;
   }


   ///////////

   int tellWireTB(uint32_t detid, float x, float y)  {
     T1DetId object(detid);
#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif 
     if(object.CSC() == 5 || object.CSC() == 2){

       //     float xx = x+_xP[object.Plane()];
       float yy = y+_yP[object.Plane()];
       int index = int(yy/3.+0.5);

#ifdef _DEBUG_
       std::cout << " yy traslato  " <<  yy << " index " << index <<std::endl;
#endif
       if(index<1 && index >= 0)index=1;
       if(index < 0 || index >WireNP[object.Plane()])return -1;
       return index;
     }
     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){
       //      float xx = x+_xG[object.Plane()];
       float yy = y+_yG[object.Plane()];
       int index = int(yy/3.+0.5);
#ifdef _DEBUG_
       std::cout << " yy traslato  " <<  yy << " index " << index <<std::endl;
#endif
       if(index<1 && index >= 0)index=1;
       if(index < 0 || index >WireNG[object.Plane()])return -1;
       return index;
     }
     return -1;
   }
   //implementare il controllo sul numero totale di strip per camera e gli offset dal tendifilo
   float tellStripA_TB(uint32_t detid, float x, float y)  {
     T1DetId object(detid);
#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif    
     if(object.CSC() == 5 || object.CSC() == 2){

       float xx = x+_xP[object.Plane()];
       float yy = y+_yP[object.Plane()];

#ifdef _DEBUG_
       std::cout << " x " << x << "   xx " << xx <<std::endl;
       std::cout << " y " << y << "   yy " << yy <<std::endl;
#endif

       //   float pitch = widthP[4]/stripNP[4];
       float pitch = parametri.Pitch();
       float index = ((yy-xx*sqrt(3.))/(2*pitch)+1 - offsetStripP[object.Plane()]/pitch);

       // -xx oppure +xx ??????????????????????????????????

       // if(index < 0 || index >WireNP[object.Plane()])return -1;
       return index;

     }

     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){
       float xx = x+_xG[object.Plane()];
       float yy = y+_yG[object.Plane()];
#ifdef _DEBUG_
       std::cout << " x " << x << "   xx " << xx <<std::endl;
       std::cout << " y " << y << "   yy " << yy <<std::endl;
#endif

       //   float pitch = widthP[4]/stripNP[4];
       float pitch = parametri.Pitch();
       float index = ((yy-xx*sqrt(3.))/(2*pitch)+1 - offsetStripG[object.Plane()]/pitch);

       // -xx oppure +xx ??????????????????????????????????

       // if(index < 0 || index >WireNG[object.Plane()])return -1;

       return index;

     }

     return -1;    
   }

   float tellStripB_TB(uint32_t detid, float x, float y)  {
     T1DetId object(detid); 
#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif
     if(object.CSC() == 5 || object.CSC() == 2){
       float xx = x+_xP[object.Plane()];
       float yy = y+_yP[object.Plane()];

       //float pitch = widthP[4]/stripNP[4];
       float pitch = parametri.Pitch();
       float index = ((yy+sqrt(3.)*(offsetP[object.Plane()]+xx))/2./pitch+1  - offsetStripP[object.Plane()]/pitch);

       // -xx oppure +xx ??????????????????????????????????

       // if(index < 0 || index >stripNP[object.Plane()])return -1;

       return index;
     }
     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){
       float xx = x+_xG[object.Plane()];
       float yy = y+_yG[object.Plane()];
       float pitch = parametri.Pitch();
       //float pitch = widthG[4]/stripNP[4];

       //    float index = ((yy+sqrt(3.)*(offsetG[object.Plane()]+xx))/(2.*pitch)+1);
       float index = ((yy+sqrt(3.)*(offsetG[object.Plane()]+xx))/2./pitch+1 - offsetStripG[object.Plane()]/pitch);
       // -xx oppure +xx ??????????????????????????????????

       // if(index < 0 || index >stripNG[object.Plane()])return -1;

       return index;
     }
     return -1;
   }


   float tellStripA_TB2(uint32_t detid, float x, float y, float discrepancy)  {
     T1DetId object(detid);
#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif    
     if(object.CSC() == 5 || object.CSC() == 2){

       float xx = x+_xP[object.Plane()];
       float yy = y+_yP[object.Plane()];

#ifdef _DEBUG_
       std::cout << " x " << x << "   xx " << xx <<std::endl;
       std::cout << " y " << y << "   yy " << yy <<std::endl;
#endif

       //   float pitch = widthP[4]/stripNP[4];
       float pitch = parametri.Pitch();
       float index = ((yy-xx*sqrt(3.))/(2*pitch)+1 - offsetStripP[object.Plane()]/pitch + discrepancy/pitch);

       // -xx oppure +xx ??????????????????????????????????

       // if(index < 0 || index >WireNP[object.Plane()])return -1;
       return index;

     }

     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){
       float xx = x+_xG[object.Plane()];
       float yy = y+_yG[object.Plane()];
#ifdef _DEBUG_
       std::cout << " x " << x << "   xx " << xx <<std::endl;
       std::cout << " y " << y << "   yy " << yy <<std::endl;
#endif

       //   float pitch = widthP[4]/stripNP[4];
       float pitch = parametri.Pitch();
       float index = ((yy-xx*sqrt(3.))/(2*pitch)+1 - offsetStripG[object.Plane()]/pitch + discrepancy/pitch);

       // -xx oppure +xx ??????????????????????????????????

       // if(index < 0 || index >WireNG[object.Plane()])return -1;

       return index;

     }

     return -1;    
   }

   float tellStripB_TB2(uint32_t detid, float x, float y, float discrepancy)  {
     T1DetId object(detid); 
#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif
     if(object.CSC() == 5 || object.CSC() == 2){
       float xx = x+_xP[object.Plane()];
       float yy = y+_yP[object.Plane()];

       //float pitch = widthP[4]/stripNP[4];
       float pitch = parametri.Pitch();
       float index = ((yy+sqrt(3.)*(offsetP[object.Plane()]+xx))/2./pitch+1  - offsetStripP[object.Plane()]/pitch + discrepancy/pitch);

       // -xx oppure +xx ??????????????????????????????????

       // if(index < 0 || index >stripNP[object.Plane()])return -1;

       return index;
     }
     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){
       float xx = x+_xG[object.Plane()];
       float yy = y+_yG[object.Plane()];
       float pitch = parametri.Pitch();
       //float pitch = widthG[4]/stripNP[4];

       //    float index = ((yy+sqrt(3.)*(offsetG[object.Plane()]+xx))/(2.*pitch)+1);
       float index = ((yy+sqrt(3.)*(offsetG[object.Plane()]+xx))/2./pitch+1 - offsetStripG[object.Plane()]/pitch + discrepancy/pitch);
       // -xx oppure +xx ??????????????????????????????????

       // if(index < 0 || index >stripNG[object.Plane()])return -1;

       return index;
     }
     return -1;
   }



   float getPitch(int plane, std::string  type)  {
   
     if( (type != "P" && type != "G") || plane < 0 || plane > 4){
       std::cout <<"Error: Ivalid Plane or CSC type " <<std::endl;
       std::cout <<"Plane: 0 to 4   CSC type: P or G"<<std::endl;
     }
     if(type=="P"){
       float pitch = parametri.Pitch();
       return pitch;
       //   return widthP[plane]/stripNP[plane];
     }
     if(type=="G"){
       float pitch = parametri.Pitch();
       return pitch;
       //  return widthG[plane]/stripNG[plane];
     }
     return -1;
   }
   
   //funzione che ritorna il gruppo inserendo il filo
   
   int GetGroup(T1DetId object, int nWire){
   
     //   T1DetId object(detid);
#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif 

     if(object.Plane() == 0){
       if(object.CSC() == 5 || object.CSC() == 2)
         {
           if( nWire <=6) return 0;
           if( nWire > 6 && nWire <=13) return 1;   
           if( nWire > 13 && nWire <=22) return 2;  
           if( nWire > 22 && nWire <=30) return 3;  
           if( nWire > 30 && nWire <=39) return 4; 
           if( nWire > 39 && nWire <=50) return 5; 
           if( nWire > 50 && nWire <=60) return 6;
           if( nWire > 60 && nWire <=72) return 7;  
           if( nWire > 72 && nWire <=84) return 8; 
           if( nWire > 84 && nWire <=96) return 9;
           if( nWire > 96 && nWire <=108) return 10;
           if( nWire > 108 && nWire <=119) return 11;
           if( nWire > 119 && nWire <=128) return 12;  //dovrebbe essere 127 ???   
         }
     
       if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 )
         {
           if( nWire <=6) return 0;
           if( nWire > 6 && nWire <=13) return 1;   
           if( nWire > 13 && nWire <=22) return 2;  
           if( nWire > 22 && nWire <=30) return 3;  
           if( nWire > 30 && nWire <=39) return 4; 
           if( nWire > 39 && nWire <=50) return 5; 
           if( nWire > 50 && nWire <=60) return 6;
           if( nWire > 60 && nWire <=72) return 7;  
           if( nWire > 72 && nWire <=84) return 8; 
           if( nWire > 84 && nWire <=96) return 9;
           if( nWire > 96 && nWire <=108) return 10;
           if( nWire > 108 && nWire <=119) return 11;
           if( nWire > 119 && nWire <=131) return 12;  
           if( nWire > 131 && nWire <=143) return 13; 
           if( nWire > 143 && nWire <=155) return 14;
           if( nWire > 155 && nWire <=166) return 15;
         }

     }
     if(object.Plane() == 1){
       if(object.CSC() == 5 || object.CSC() == 2)
         {
           if( nWire <=7) return 0;
           if( nWire > 7 && nWire <=15) return 1;   
           if( nWire > 15 && nWire <=25) return 2;  
           if( nWire > 25 && nWire <=33) return 3;  
           if( nWire > 33 && nWire <=44) return 4; 
           if( nWire > 44 && nWire <=55) return 5; 
           if( nWire > 55 && nWire <=67) return 6;
           if( nWire > 67 && nWire <=80) return 7;  
           if( nWire > 80 && nWire <=92) return 8; 
           if( nWire > 92 && nWire <=105) return 9;
           if( nWire > 105 && nWire <=118) return 10;
           if( nWire > 118 && nWire <=131) return 11;
           if( nWire > 131 && nWire <=144) return 12;
           if( nWire > 144 && nWire <=151) return 13;    
     
         }
     
       if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 )
         {
           if( nWire <=7) return 0;
           if( nWire > 7 && nWire <=15) return 1;   
           if( nWire > 15 && nWire <=25) return 2;  
           if( nWire > 25 && nWire <=33) return 3;  
           if( nWire > 33 && nWire <=44) return 4; 
           if( nWire > 44 && nWire <=55) return 5; 
           if( nWire > 55 && nWire <=67) return 6;
           if( nWire > 67 && nWire <=80) return 7;  
           if( nWire > 80 && nWire <=92) return 8; 
           if( nWire > 92 && nWire <=105) return 9;
           if( nWire > 105 && nWire <=118) return 10;
           if( nWire > 118 && nWire <=131) return 11;
           if( nWire > 131 && nWire <=144) return 12;
           if( nWire > 144 && nWire <=157) return 13;    
           if( nWire > 157 && nWire <=170) return 14;
           if( nWire > 170 && nWire <=181) return 15;
         }

     }
     if(object.Plane() == 2){
       if(object.CSC() == 5 || object.CSC() == 2)
         {
           if( nWire <=9) return 0;
           if( nWire > 9 && nWire <=17) return 1;   
           if( nWire > 17 && nWire <=27) return 2;  
           if( nWire > 27 && nWire <=37) return 3;  
           if( nWire > 37 && nWire <=48) return 4; 
           if( nWire > 48 && nWire <=60) return 5; 
           if( nWire > 60 && nWire <=72) return 6;
           if( nWire > 72 && nWire <=86) return 7;  
           if( nWire > 86 && nWire <=100) return 8; 
           if( nWire > 100 && nWire <=114) return 9;
           if( nWire > 114 && nWire <=128) return 10;
           if( nWire > 128 && nWire <=142) return 11;
           if( nWire > 142 && nWire <=156) return 12;
           if( nWire > 156 && nWire <=159) return 13;  //dovrebbero essere 158    
     
         }
     
       if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 )
         {
           if( nWire <=9) return 0;
           if( nWire > 9 && nWire <=17) return 1;   
           if( nWire > 17 && nWire <=27) return 2;  
           if( nWire > 27 && nWire <=37) return 3;  
           if( nWire > 37 && nWire <=48) return 4; 
           if( nWire > 48 && nWire <=60) return 5; 
           if( nWire > 60 && nWire <=72) return 6;
           if( nWire > 72 && nWire <=86) return 7;  
           if( nWire > 86 && nWire <=100) return 8; 
           if( nWire > 100 && nWire <=114) return 9;
           if( nWire > 114 && nWire <=128) return 10;
           if( nWire > 128 && nWire <=142) return 11;
           if( nWire > 142 && nWire <=156) return 12;
           if( nWire > 156 && nWire <=169) return 13;    
           if( nWire > 169 && nWire <=184) return 14;
           if( nWire > 184 && nWire <=197) return 15;
         }

     }
     if(object.Plane() == 3){
       if(object.CSC() == 5 || object.CSC() == 2)
         {
           if( nWire <=9) return 0;
           if( nWire > 9 && nWire <=19) return 1;   
           if( nWire > 19 && nWire <=29) return 2;  
           if( nWire > 29 && nWire <=40) return 3;  
           if( nWire > 40 && nWire <=51) return 4; 
           if( nWire > 51 && nWire <=64) return 5; 
           if( nWire > 64 && nWire <=78) return 6;
           if( nWire > 78 && nWire <=92) return 7;  
           if( nWire > 92 && nWire <=107) return 8; 
           if( nWire > 107 && nWire <=122) return 9;
           if( nWire > 122 && nWire <=137) return 10;
           if( nWire > 137 && nWire <=152) return 11;
           if( nWire > 152 && nWire <=166) return 12;
           if( nWire > 166 && nWire <=181) return 13;    
     
         }
     
       if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 )
         {
           if( nWire <=9) return 0;
           if( nWire > 9 && nWire <=19) return 1;   
           if( nWire > 19 && nWire <=29) return 2;  
           if( nWire > 29 && nWire <=40) return 3;  
           if( nWire > 40 && nWire <=51) return 4; 
           if( nWire > 51 && nWire <=64) return 5; 
           if( nWire > 64 && nWire <=78) return 6;
           if( nWire > 78 && nWire <=92) return 7;  
           if( nWire > 92 && nWire <=107) return 8; 
           if( nWire > 107 && nWire <=122) return 9;
           if( nWire > 122 && nWire <=137) return 10;
           if( nWire > 137 && nWire <=152) return 11;
           if( nWire > 152 && nWire <=166) return 12;
           if( nWire > 166 && nWire <=181) return 13;    
           if( nWire > 181 && nWire <=197) return 14;
           if( nWire > 197 && nWire <=214) return 15;
         }

     }
     if(object.Plane() == 4){
       if(object.CSC() == 5 || object.CSC() == 2)
         {
           if( nWire <=12) return 0;
           if( nWire > 12 && nWire <=22) return 1;   
           if( nWire > 22 && nWire <=34) return 2;  
           if( nWire > 34 && nWire <=44) return 3;  
           if( nWire > 44 && nWire <=57) return 4; 
           if( nWire > 57 && nWire <=71) return 5; 
           if( nWire > 71 && nWire <=86) return 6;
           if( nWire > 86 && nWire <=102) return 7;  
           if( nWire > 102 && nWire <=118) return 8; 
           if( nWire > 118 && nWire <=134) return 9;
           if( nWire > 134 && nWire <=150) return 10;
           if( nWire > 150 && nWire <=166) return 11;
           if( nWire > 166 && nWire <=182) return 12;
           if( nWire > 182 && nWire <=198) return 13;    
           if( nWire > 198 && nWire <=205) return 14;    
     
         }
     
       if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 )
         {
           if( nWire <=12) return 0;
           if( nWire > 12 && nWire <=22) return 1;   
           if( nWire > 22 && nWire <=34) return 2;  
           if( nWire > 34 && nWire <=44) return 3;  
           if( nWire > 44 && nWire <=57) return 4; 
           if( nWire > 57 && nWire <=71) return 5; 
           if( nWire > 71 && nWire <=86) return 6;
           if( nWire > 86 && nWire <=102) return 7;  
           if( nWire > 102 && nWire <=118) return 8; 
           if( nWire > 118 && nWire <=134) return 9;
           if( nWire > 134 && nWire <=150) return 10;
           if( nWire > 150 && nWire <=166) return 11;
           if( nWire > 166 && nWire <=182) return 12;
           if( nWire > 182 && nWire <=198) return 13;    
           if( nWire > 198 && nWire <=214) return 14;
           if( nWire > 214 && nWire <=227) return 15;
         }

     }

     return -1;
   }


   int GetGroup(uint32_t ID, int nWire){
   
     T1DetId object(ID);
#ifdef _DEBUG_
     std::cout << " Camera " <<  object.Plane() << " " << object.CSC()<<std::endl;
#endif 

     if(object.Plane() == 0){
       if(object.CSC() == 5 || object.CSC() == 2)
         {
           if( nWire <=6) return 0;
           if( nWire > 6 && nWire <=13) return 1;   
           if( nWire > 13 && nWire <=22) return 2;  
           if( nWire > 22 && nWire <=30) return 3;  
           if( nWire > 30 && nWire <=39) return 4; 
           if( nWire > 39 && nWire <=50) return 5; 
           if( nWire > 50 && nWire <=60) return 6;
           if( nWire > 60 && nWire <=72) return 7;  
           if( nWire > 72 && nWire <=84) return 8; 
           if( nWire > 84 && nWire <=96) return 9;
           if( nWire > 96 && nWire <=108) return 10;
           if( nWire > 108 && nWire <=119) return 11;
           if( nWire > 119 && nWire <=128) return 12;  //dovrebbe essere 127 ???   
         }
     
       if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 )
         {
           if( nWire <=6) return 0;
           if( nWire > 6 && nWire <=13) return 1;   
           if( nWire > 13 && nWire <=22) return 2;  
           if( nWire > 22 && nWire <=30) return 3;  
           if( nWire > 30 && nWire <=39) return 4; 
           if( nWire > 39 && nWire <=50) return 5; 
           if( nWire > 50 && nWire <=60) return 6;
           if( nWire > 60 && nWire <=72) return 7;  
           if( nWire > 72 && nWire <=84) return 8; 
           if( nWire > 84 && nWire <=96) return 9;
           if( nWire > 96 && nWire <=108) return 10;
           if( nWire > 108 && nWire <=119) return 11;
           if( nWire > 119 && nWire <=131) return 12;  
           if( nWire > 131 && nWire <=143) return 13; 
           if( nWire > 143 && nWire <=155) return 14;
           if( nWire > 155 && nWire <=166) return 15;
         }

     }
     if(object.Plane() == 1){
       if(object.CSC() == 5 || object.CSC() == 2)
         {
           if( nWire <=7) return 0;
           if( nWire > 7 && nWire <=15) return 1;   
           if( nWire > 15 && nWire <=25) return 2;  
           if( nWire > 25 && nWire <=33) return 3;  
           if( nWire > 33 && nWire <=44) return 4; 
           if( nWire > 44 && nWire <=55) return 5; 
           if( nWire > 55 && nWire <=67) return 6;
           if( nWire > 67 && nWire <=80) return 7;  
           if( nWire > 80 && nWire <=92) return 8; 
           if( nWire > 92 && nWire <=105) return 9;
           if( nWire > 105 && nWire <=118) return 10;
           if( nWire > 118 && nWire <=131) return 11;
           if( nWire > 131 && nWire <=144) return 12;
           if( nWire > 144 && nWire <=151) return 13;    
     
         }
     
       if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 )
         {
           if( nWire <=7) return 0;
           if( nWire > 7 && nWire <=15) return 1;   
           if( nWire > 15 && nWire <=25) return 2;  
           if( nWire > 25 && nWire <=33) return 3;  
           if( nWire > 33 && nWire <=44) return 4; 
           if( nWire > 44 && nWire <=55) return 5; 
           if( nWire > 55 && nWire <=67) return 6;
           if( nWire > 67 && nWire <=80) return 7;  
           if( nWire > 80 && nWire <=92) return 8; 
           if( nWire > 92 && nWire <=105) return 9;
           if( nWire > 105 && nWire <=118) return 10;
           if( nWire > 118 && nWire <=131) return 11;
           if( nWire > 131 && nWire <=144) return 12;
           if( nWire > 144 && nWire <=157) return 13;    
           if( nWire > 157 && nWire <=170) return 14;
           if( nWire > 170 && nWire <=181) return 15;
         }

     }
     if(object.Plane() == 2){
       if(object.CSC() == 5 || object.CSC() == 2)
         {
           if( nWire <=9) return 0;
           if( nWire > 9 && nWire <=17) return 1;   
           if( nWire > 17 && nWire <=27) return 2;  
           if( nWire > 27 && nWire <=37) return 3;  
           if( nWire > 37 && nWire <=48) return 4; 
           if( nWire > 48 && nWire <=60) return 5; 
           if( nWire > 60 && nWire <=72) return 6;
           if( nWire > 72 && nWire <=86) return 7;  
           if( nWire > 86 && nWire <=100) return 8; 
           if( nWire > 100 && nWire <=114) return 9;
           if( nWire > 114 && nWire <=128) return 10;
           if( nWire > 128 && nWire <=142) return 11;
           if( nWire > 142 && nWire <=156) return 12;
           if( nWire > 156 && nWire <=159) return 13;  //dovrebbero essere 158    
     
         }
     
       if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 )
         {
           if( nWire <=9) return 0;
           if( nWire > 9 && nWire <=17) return 1;   
           if( nWire > 17 && nWire <=27) return 2;  
           if( nWire > 27 && nWire <=37) return 3;  
           if( nWire > 37 && nWire <=48) return 4; 
           if( nWire > 48 && nWire <=60) return 5; 
           if( nWire > 60 && nWire <=72) return 6;
           if( nWire > 72 && nWire <=86) return 7;  
           if( nWire > 86 && nWire <=100) return 8; 
           if( nWire > 100 && nWire <=114) return 9;
           if( nWire > 114 && nWire <=128) return 10;
           if( nWire > 128 && nWire <=142) return 11;
           if( nWire > 142 && nWire <=156) return 12;
           if( nWire > 156 && nWire <=169) return 13;    
           if( nWire > 169 && nWire <=184) return 14;
           if( nWire > 184 && nWire <=197) return 15;
         }

     }
     if(object.Plane() == 3){
       if(object.CSC() == 5 || object.CSC() == 2)
         {
           if( nWire <=9) return 0;
           if( nWire > 9 && nWire <=19) return 1;   
           if( nWire > 19 && nWire <=29) return 2;  
           if( nWire > 29 && nWire <=40) return 3;  
           if( nWire > 40 && nWire <=51) return 4; 
           if( nWire > 51 && nWire <=64) return 5; 
           if( nWire > 64 && nWire <=78) return 6;
           if( nWire > 78 && nWire <=92) return 7;  
           if( nWire > 92 && nWire <=107) return 8; 
           if( nWire > 107 && nWire <=122) return 9;
           if( nWire > 122 && nWire <=137) return 10;
           if( nWire > 137 && nWire <=152) return 11;
           if( nWire > 152 && nWire <=166) return 12;
           if( nWire > 166 && nWire <=181) return 13;    
     
         }
     
       if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 )
         {
           if( nWire <=9) return 0;
           if( nWire > 9 && nWire <=19) return 1;   
           if( nWire > 19 && nWire <=29) return 2;  
           if( nWire > 29 && nWire <=40) return 3;  
           if( nWire > 40 && nWire <=51) return 4; 
           if( nWire > 51 && nWire <=64) return 5; 
           if( nWire > 64 && nWire <=78) return 6;
           if( nWire > 78 && nWire <=92) return 7;  
           if( nWire > 92 && nWire <=107) return 8; 
           if( nWire > 107 && nWire <=122) return 9;
           if( nWire > 122 && nWire <=137) return 10;
           if( nWire > 137 && nWire <=152) return 11;
           if( nWire > 152 && nWire <=166) return 12;
           if( nWire > 166 && nWire <=181) return 13;    
           if( nWire > 181 && nWire <=197) return 14;
           if( nWire > 197 && nWire <=214) return 15;
         }

     }
     if(object.Plane() == 4){
       if(object.CSC() == 5 || object.CSC() == 2)
         {
           if( nWire <=12) return 0;
           if( nWire > 12 && nWire <=22) return 1;   
           if( nWire > 22 && nWire <=34) return 2;  
           if( nWire > 34 && nWire <=44) return 3;  
           if( nWire > 44 && nWire <=57) return 4; 
           if( nWire > 57 && nWire <=71) return 5; 
           if( nWire > 71 && nWire <=86) return 6;
           if( nWire > 86 && nWire <=102) return 7;  
           if( nWire > 102 && nWire <=118) return 8; 
           if( nWire > 118 && nWire <=134) return 9;
           if( nWire > 134 && nWire <=150) return 10;
           if( nWire > 150 && nWire <=166) return 11;
           if( nWire > 166 && nWire <=182) return 12;
           if( nWire > 182 && nWire <=198) return 13;    
           if( nWire > 198 && nWire <=205) return 14;    
     
         }
     
       if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 )
         {
           if( nWire <=12) return 0;
           if( nWire > 12 && nWire <=22) return 1;   
           if( nWire > 22 && nWire <=34) return 2;  
           if( nWire > 34 && nWire <=44) return 3;  
           if( nWire > 44 && nWire <=57) return 4; 
           if( nWire > 57 && nWire <=71) return 5; 
           if( nWire > 71 && nWire <=86) return 6;
           if( nWire > 86 && nWire <=102) return 7;  
           if( nWire > 102 && nWire <=118) return 8; 
           if( nWire > 118 && nWire <=134) return 9;
           if( nWire > 134 && nWire <=150) return 10;
           if( nWire > 150 && nWire <=166) return 11;
           if( nWire > 166 && nWire <=182) return 12;
           if( nWire > 182 && nWire <=198) return 13;    
           if( nWire > 198 && nWire <=214) return 14;
           if( nWire > 214 && nWire <=227) return 15;
         }

     }

     return 0;
   }

   float get_bmin(int layer,char size){
     if (size=='G') return bminG[layer]; else return bminP[layer];
   }
         
   float get_bmax(int layer,char size){
     if (size=='G') return bmaxG[layer]; else return bmaxP[layer];
   }

   float get_Rmin(int layer){
     return RminSens[layer];
   }

   float get_Rot(int layer){
     return Rot0[layer];
   }

   float get_Z(int layer, int part){
     return Z[layer][part];
   }

   float get_disX(int part){
     return disassamentoX[part];
   }

   int get_wire(int layer,char size){
     if (size=='G') return WireNG[layer]; else return WireNP[layer];
   }

   float get_h(int layer,char size){
     if (size=='G') return 2*_yG[layer]; else return 2*_yP[layer];
   }


   float get_offsetGP(uint32_t detid){
     T1DetId object(detid);

     if(object.CSC() == 5 || object.CSC() == 2){
      
       return offsetP[object.Plane()];

     }

     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){ 
       
       return offsetG[object.Plane()];
       
     }

   }
   
   float get_offsetstripGP(uint32_t detid){
     T1DetId object(detid);

     if(object.CSC() == 5 || object.CSC() == 2){
  
       return offsetStripP[object.Plane()];

     }

     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){ 
       return offsetStripG[object.Plane()];
  
     }

   }


   float get_xGP(uint32_t detid){
     T1DetId object(detid);

     if(object.CSC() == 5 || object.CSC() == 2){
  
       return _xP[object.Plane()];
     }

     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){ 

       return _xG[object.Plane()];
     }

   }

   float get_yGP(uint32_t detid){
     T1DetId object(detid);

     if(object.CSC() == 5 || object.CSC() == 2){
  
       return _yP[object.Plane()];
     }

     if(object.CSC() == 0 || object.CSC() == 1|| object.CSC() == 4|| object.CSC() == 3 ){ 

       return _yG[object.Plane()];
     }

   }

   
  private:

   T1ChamberSpecs parametri; 

   float _xP[5];
   float _xG[5];
   float _yP[5];
   float _yG[5];
   float widthP[5];
   float offsetStripG[5];
   float offsetStripP[5];
   int stripNG[5];
   float widthG[5];
   int WireNP[5];
   int WireNG[5];
   int stripNP[5];
   float bmaxP[5];
   float bminP[5];
   float bminOldP[5];
   float bmaxG[5];
   float bminG[5];
   float bminOldG[5];
   float offsetP[5];
   float offsetG[5];


   float RminSens[5];
   float Z[5][6];

   float Rot0[5];
   float disassamentoX[6];



   double Align_x[2][5][6];
   double Align_y[2][5][6];
   double Sex_Align_x[2][6];
   double Sex_Align_y[2][6];
   double Align_th[2][5][6];
   double Arm_Align_x[2];
   double Arm_Align_y[2];
   double Arm_Align_z[2];

};




#endif
