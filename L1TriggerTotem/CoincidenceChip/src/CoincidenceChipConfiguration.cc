#include "L1TriggerTotem/CoincidenceChip/interface/CoincidenceChipConfiguration.h"


CoincidenceChipConfiguration::CoincidenceChipConfiguration() {
  verbose_ = false; 
  controlRegisters_.reset();
}

CoincidenceChipConfiguration::~CoincidenceChipConfiguration() {
}

void CoincidenceChipConfiguration::configure(const edm::ParameterSet& iConfig){

  //int controlRegister1=iConfig.getParameter<int>("controlRegister1"); 
  //int controlRegister2=iConfig.getParameter<int>("controlRegister2"); 
  //int controlRegister3=iConfig.getParameter<int>("controlRegister3"); 
  bool useControlRegisters=iConfig.getParameter<bool>("useControlRegisters"); 
  //int test=iConfig.getParameter<int>("controlRegister1",-1); 
  if(useControlRegisters){
      //if((test1>=0) && (test2>=0) && (test3>=1)){
      if(verbose_) edm::LogInfo("CChipConfiguration") << "Configuretion of CC via controlRegisters{1,2,3} ..." << std::endl;
      setControlRegister1(iConfig.getParameter<unsigned int>("controlRegister1")); 
      setControlRegister2(iConfig.getParameter<unsigned int>("controlRegister2")); 
      setControlRegister3(iConfig.getParameter<unsigned int>("controlRegister3")); 
      //}else if((test1>=0) || (test2>=0) || (test3>=0)){
      //}else{
      //std::cout << "\nCoincidenceChipConfiguration::configure" << std::endl;
      //std::cout << "    for configuring CoincidenceChip use ALL control registers controlRegister{1,2,3} XOR numbers V,NP,OV,W...!!! \n" << std::endl;
      //assert(0);
      //}
  }else{
    if(verbose_) edm::LogInfo("CChipConfiguration") << "Configuretion of CC via numbers V,NP,OV,W ..." << std::endl;
    setV(iConfig.getParameter<unsigned int>("V")); 
    setNP(iConfig.getParameter<unsigned int>("NP")); 
    setOV(iConfig.getParameter<unsigned int>("OV")); 
    setW(iConfig.getParameter<unsigned int>("W")); 
    setZ(iConfig.getParameter<unsigned int>("Z")); 
    setO2(iConfig.getParameter<unsigned int>("O2")); 
    setLI(iConfig.getParameter<unsigned int>("LI")); 
    setLO(iConfig.getParameter<unsigned int>("LO")); 
    setAO(iConfig.getParameter<unsigned int>("AO")); 
  }

  unsigned int _useLogicWithWrongNP = iConfig.getParameter<unsigned int>("useLogicWithWrongNP");
  if( _useLogicWithWrongNP==1)
	  useLogicWithWrongNP();
  else
	  useLogicWithCorrectNP();

#if 0
  std::cout << configSummary();

  setControlRegister1(83);
  std::cout << configSummary();

  setControlRegister2(242);
  std::cout << configSummary();

  setControlRegister3(200);
  std::cout << configSummary();
#endif
}

void CoincidenceChipConfiguration::setControlReg(std::bitset<24> controlRegisters) {
  controlRegisters_ = controlRegisters;
}

std::bitset<24> CoincidenceChipConfiguration::getControlReg() const{
  return controlRegisters_;
}

/********* logic configuration **********/

void CoincidenceChipConfiguration::useLogicWithWrongNP(){
  useLogicWithWrongNP_ = 1;
}

void CoincidenceChipConfiguration::useLogicWithCorrectNP(){
  useLogicWithWrongNP_ = 0;
}

bool CoincidenceChipConfiguration::getLogicWithWrongNPFlag()const{
  if( useLogicWithWrongNP_ )
	  return true;
  else
	  return false;
}

/************  V **************/

void CoincidenceChipConfiguration::setV(unsigned short v_) {
  controlRegisters_.set(0, (v_ & 0x1));
  controlRegisters_.set(1, (v_ & 0x2));
  controlRegisters_.set(2, (v_ & 0x4));
  controlRegisters_.set(3, (v_ & 0x8));
}

unsigned short CoincidenceChipConfiguration::getV() const{
  return (unsigned short) (controlRegisters_ & std::bitset<24>(0x0000000FUL)).to_ulong();
}

/************  NP **************/

void CoincidenceChipConfiguration::setNP(unsigned short np_) {
  controlRegisters_.set(4, (np_ & 0x1));
}

unsigned short CoincidenceChipConfiguration::getNP() const {
  return (controlRegisters_[4] ? 1 : 0);
}

/************  OV **************/

void CoincidenceChipConfiguration::setOV(unsigned short ov_) {
  controlRegisters_.set(5, (ov_ & 0x1));
  controlRegisters_.set(6, (ov_ & 0x2));
  controlRegisters_.set(7, (ov_ & 0x4));
}

unsigned short CoincidenceChipConfiguration::getOV() const{
  return (unsigned short) ((controlRegisters_ & std::bitset<24>(0x000000E0UL))>>=5).to_ulong();
}

/************  W **************/

void CoincidenceChipConfiguration::setW(unsigned short w_) {
  controlRegisters_.set(8, (w_ & 0x1));
  controlRegisters_.set(9, (w_ & 0x2));
  controlRegisters_.set(10, (w_ & 0x4));
  controlRegisters_.set(11, (w_ & 0x8));
}

unsigned short CoincidenceChipConfiguration::getW() const{
  return (unsigned short) ((controlRegisters_ & std::bitset<24>(0x00000F00UL))>>=8).to_ulong();
}

/************  Z **************/

void CoincidenceChipConfiguration::setZ(unsigned short z_) {
  controlRegisters_.set(12, (z_ & 0x1));
  controlRegisters_.set(13, (z_ & 0x2));
  controlRegisters_.set(14, (z_ & 0x4));
  controlRegisters_.set(15, (z_ & 0x8));
}

unsigned short CoincidenceChipConfiguration::getZ() const{
  return (unsigned short) ((controlRegisters_ & std::bitset<24>(0x0000F000UL))>>12).to_ulong();
}

/************ O2 **************/

void CoincidenceChipConfiguration::setO2(unsigned short o2_) {
  controlRegisters_.set(16, (o2_ & 0x1));
  controlRegisters_.set(17, (o2_ & 0x2));
  controlRegisters_.set(18, (o2_ & 0x4));
}

unsigned short CoincidenceChipConfiguration::getO2() const{
  return (unsigned short) ((controlRegisters_ & std::bitset<24>(0x00070000UL))>>16).to_ulong();
}

/************ LI **************/

void CoincidenceChipConfiguration::setLI(unsigned short li_) {
  controlRegisters_.set(19, (li_ & 0x1));
}

unsigned short CoincidenceChipConfiguration::getLI() const{
  return (controlRegisters_[19] ? 1 : 0);
}

/************ LO **************/

void CoincidenceChipConfiguration::setLO(unsigned short lo_) {
  controlRegisters_.set(20, (lo_ & 0x1));
  controlRegisters_.set(21, (lo_ & 0x2));
}

unsigned short CoincidenceChipConfiguration::getLO() const{
  return (unsigned short) ((controlRegisters_ & std::bitset<24>(0x00300000UL))>>20).to_ulong();
}

/************ AO **************/

void CoincidenceChipConfiguration::setAO(unsigned short ao_) {
  controlRegisters_.set(22, (ao_ & 0x1));
  controlRegisters_.set(23, (ao_ & 0x2));
}

unsigned short CoincidenceChipConfiguration::getAO() const{
  return (unsigned short) ((controlRegisters_ & std::bitset<24>(0x00C00000UL))>>22).to_ulong();
}

void CoincidenceChipConfiguration::CheckControlRegister(unsigned short regNumber, unsigned int n) const{
    if(n>=256){
      TString msg = "CoincidenceChipConfiguration::setControlRegister";
      msg+=regNumber;
      throw cms::Exception((std::string)msg) << "Control register " << regNumber << " must be <=" << 255 << " (not" <<  n <<")!!!" << std::endl 
                                << "Check configuration of CC..." << std::endl;
    }
}

/************ control register 1 **************/
unsigned short CoincidenceChipConfiguration::getControlRegister1() const{
   return (unsigned short) ((controlRegisters_ & std::bitset<24>(0x000000FFUL))).to_ulong();
}

void CoincidenceChipConfiguration::setControlRegister1(unsigned long n){
    CheckControlRegister(1,n);
    
    std::bitset<8> reg(n);
	setV( (reg & std::bitset<8>(0x0FUL)).to_ulong());
	setNP(((reg & std::bitset<8>(0x10UL))>>4).to_ulong());
	setOV(((reg & std::bitset<8>(0xE0UL))>>5).to_ulong());
}

/************ control register 2 **************/
unsigned short CoincidenceChipConfiguration::getControlRegister2() const{
   return (unsigned short) ((controlRegisters_ & std::bitset<24>(0x0000FF00UL))>>8).to_ulong();
}

void CoincidenceChipConfiguration::setControlRegister2(unsigned long n){
    CheckControlRegister(2,n);

    std::bitset<8> reg(n);
	setW((reg & std::bitset<8>(0x0FUL)).to_ulong());
	setZ(((reg & std::bitset<8>(0xF0UL))>>4).to_ulong());
}

/************ control register 3 **************/
unsigned short CoincidenceChipConfiguration::getControlRegister3() const{
   return (unsigned short) ((controlRegisters_ & std::bitset<24>(0x00FF0000UL))>>16).to_ulong();
}

void CoincidenceChipConfiguration::setControlRegister3(unsigned long n){
    CheckControlRegister(3,n);

    std::bitset<8> reg(n);
	setO2(( reg & std::bitset<8>(0x07UL)).to_ulong());
	setLI(((reg & std::bitset<8>(0x08UL))>>3).to_ulong());
	setLO(((reg & std::bitset<8>(0x30UL))>>4).to_ulong());
	setAO(((reg & std::bitset<8>(0xC0UL))>>6).to_ulong());
}

/************ Print configuration of CC **************/
std::string CoincidenceChipConfiguration::configSummary() const {
	using namespace std;
	char res2[200];
	stringstream out;
	out << "CC configuration:" << endl;
	// char res1[200];
	// sprintf(res1, "   V=0x%x NP=0x%x OV=0x%x W=0x%x Z=0x%x O2=0x%x LI=0x%x LO=0x%x AO=0x%x", getV(), getNP(), getOV(), getW(), getZ(), getO2(), getLI(), getLO(), getAO());
    // out << res1 << endl;
	sprintf(res2, "   V=%u NP=%u OV=%u W=%u Z=%u O2=%u LI=%u LO=%u AO=%u", getV(), getNP(), getOV(), getW(), getZ(), getO2(), getLI(), getLO(), getAO());
    out << res2 << endl;
	out << "   controlRegister1    = " << getControlRegister1()     << endl;
	out << "   controlRegister2    = " << getControlRegister2()     << endl; 
	out << "   controlRegister3    = " << getControlRegister3()     << endl;
	out << "   useLogicWithWrongNP = " << getLogicWithWrongNPFlag() << endl;
    //out << "useControlRegisters =" << useControlRegisters << endl;
	return  out.str();
}

//************ summary **************/
//std::string CoincidenceChipConfiguration::summary() const{
//  char res[200];
//  sprintf(res, "V=0x%x NP=0x%x OV=0x%x W=0x%x Z=0x%x O2=0x%x LI=0x%x LO=0x%x AO=0x%x \n", getV(), getNP(), getOV(), getW(), getZ(), getO2(), getLI(), getLO(),
//      getAO());
//  return std::string(res);
//}

