/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *  Leszek Grzanka (Leszek.Grzanka@cern.ch)
 *
 * TODO replace cout and cerr by logger interface
 * TODO implement a method of providing user defined process description
 *
 ****************************************************************************/

#include "IOMC/Phojet/interface/Phojet.h"

// list files in given directory
void show_files( const path & directory, bool recurse_into_subdirs = true ){
  if( exists( directory ) ){
    directory_iterator end ;

    for( directory_iterator iter(directory) ; iter != end ; ++iter ){
      if ( is_directory( *iter ) ){
        cout << iter->path() << " (directory)" << endl;
        if( recurse_into_subdirs ) show_files(*iter) ;
      } else {
        cout << iter->path() << " (file)" << endl;
      } // endif

    } // end for

  } // end if
}

//----------------------------------------------------------------------------------------------------
int Phojet::generate_input_card( const path & file, double cms_energy_GeV, long number_of_events, string phojet_process_description){

  /** NA1,NA2,NA3,NB1  - values for initializing the generator
                            NA? must be in 1..178 and not all 1;
                            12,34,56  are the standard values
                            NB1 must be in 1..168;
                            78  is the standard value */

  int NA1 = 1 + (int)(177.0 * (rndEng->flat()));
  int NA2 = 1 + (int)(177.0 * (rndEng->flat()));
  int NA3 = 2 + (int)(176.0 * (rndEng->flat()));
  int NB1 = 1 + (int)(167.0 * (rndEng->flat()));

  boost::filesystem::ofstream phojet_cfg(file);

  phojet_cfg << "INIT-RNDM " << NA1 << " " << NA2 << " " << NA3 << " " << NB1 << " \n";
  phojet_cfg << "PROCESS     " << phojet_process_description << " \n";
  phojet_cfg << "PARTICLE1   2212    0.0 \n";
  phojet_cfg << "PARTICLE2   2212    0.0 \n";
  phojet_cfg << "SETMODEL   22  2 \n";
  phojet_cfg << "EVENT-CMS   " << cms_energy_GeV << "   " << number_of_events <<  " \n";
  phojet_cfg << "ENDINPUT \n";

  phojet_cfg.close();
  return 0;
}


//----------------------------------------------------------------------------------------------------
void Phojet::run_phojet( const path & workDirectory, const path & phojetExeFile, const path & phojetCfgFile, const bool printCrossSection){

  string phojetExeBasename = boost::filesystem::basename(phojetExeFile);

  string command = "cd " + workDirectory.string() + " && pwd && ./" + phojetExeBasename + " < " + phojetCfgFile.string();

  // go to the tmp directory, and run phojet executable
  FILE * fp = popen(command.c_str(), "r");
  if (fp == NULL){
    cerr << " Phojet::run_phojet - problem running command [ " << command << " ]" << endl;
  }

  // read what phojet prints to the stdout
  char outputStream[LINE_MAX];
  while (fgets(outputStream, LINE_MAX, fp) != NULL){
    string outputLine(outputStream);
    if( printCrossSection ){
      if( outputLine.find("max. cross section (mb)") != string::npos){
        cout << outputLine << endl;
      }
    }
  }

  int status = pclose(fp);
  if (status == -1) {
      /* Error reported by pclose() */
      cerr << "Phojet::run_phojet , phojet finished with error\n" << endl;
  } else {
      /* Use macros described under wait() to inspect `status' in order
         to determine success/failure of command executed by popen() */
  }
}

//----------------------------------------------------------------------------------------------------
Phojet::Phojet(const ParameterSet& pSet) :
  verbosity(pSet.getUntrackedParameter<unsigned int>("verbosity", 1))
{
  produces<HepMCProduct>();
  if( verbosity > 10 )
    cout << "Phojet::Phojet()" << endl;

  // initialize random engine
  Service<RandomNumberGenerator> rng;
  rndEng = &(rng->getEngine());
  if( verbosity > 10 )
    cout << ">> Phojet > seed = " << rndEng->getSeed() << endl;

  buffer_size = pSet.getUntrackedParameter<unsigned int>("bufferSize", 10);

  cms_energy = pSet.getParameter<double>("cmsEnergy");

  phojet_executable = pSet.getParameter<string>("phojetExecutable");

  reader_ = HepMCFileReader::instance();

  string process_str = pSet.getParameter<string>("process");

  if( process_str == "DPE" ){
    phojet_process_description = "0 0 0 1 0 0 0 0";
  } else if ( process_str == "DD"){
    phojet_process_description = "0 0 0 0 0 0 1 0";
  } else if ( process_str == "MB"){
    phojet_process_description = "1 0 0 0 0 0 0 0";
  } else if ( process_str == "SD"){
    phojet_process_description = "0 0 0 0 1 1 0 0";
  } else {
    phojet_process_description = process_str;
    if( process_str.length() != 15 ){
      cout << "Wrong process description " << process_str << endl;
      cout << "Correct process description follows pattern \"[0|1] [0|1] [0|1] [0|1] [0|1] [0|1] [0|1] [0|1]\"" << endl;
      phojet_process_description = "";
      return;
    } else {
      int i;
      for( i = 1 ; i < 14 ; i+= 2){
        if( process_str[i] != ' ' ){
          cout << "Wrong process description " << process_str << endl;
          cout << "Correct process description follows pattern \"[0|1] [0|1] [0|1] [0|1] [0|1] [0|1] [0|1] [0|1]\"" << endl;
          phojet_process_description = "";
          return;
        }
      }
      for( i = 0 ; i < 15 ; i+= 2){
        if( (process_str[i] != '0') && (process_str[i] != '1')){
          cout << "Wrong process description " << process_str << endl;
          cout << "Correct process description follows pattern \"[0|1] [0|1] [0|1] [0|1] [0|1] [0|1] [0|1] [0|1]\"" << endl;
          phojet_process_description = "";
          return;
        }
      }
    }
    cout << "User defined process " << process_str << endl;
  }

  phojetStatus = 0;
}

//----------------------------------------------------------------------------------------------------
Phojet::~Phojet()
{
  if( verbosity > 10 )
    cout << "Phojet::~Phojet()" << endl;
}

//----------------------------------------------------------------------------------------------------
void Phojet::beginJob()
{
  if( verbosity > 10 )
    cout << "Phojet::beginJob()" << endl;

  // generate tmp directory name
  tmpDirectoryPath = path( tmpnam (NULL) );

  // create tmp directory
  create_directory(tmpDirectoryPath);
  if( verbosity > 10 )
    cout << "Directory " << tmpDirectoryPath << " created" << endl;

  // copy phojet executable to tmp dir
  path phojetExecutableFileName = path(getenv("CMSSW_BASE")) / "src" / phojet_executable;
  path tmpExeFile = tmpDirectoryPath / boost::filesystem::basename(phojetExecutableFileName);

  copy_file( phojetExecutableFileName,  tmpExeFile);

  if( verbosity > 10 )
    cout << "Copy " << phojetExecutableFileName << " to " << tmpExeFile << endl;

  if( verbosity > 10 )
    show_files( tmpDirectoryPath );
}

//----------------------------------------------------------------------------------------------------
void Phojet::produce(edm::Event &e, const edm::EventSetup &es)
{  
  // create event structure
  GenEvent* gEv = new GenEvent();
  gEv->set_event_number(e.id().event());

  if( verbosity > 10 )
    cout << "Phojet::produce() , event id: " << e.id().event() << endl;

  path outputFilePath = tmpDirectoryPath / "fredrik.txt";

  path phojetExecutableFileName = path(getenv("CMSSW_BASE")) / "src" / phojet_executable;
  path tmpExeFile = tmpDirectoryPath / boost::filesystem::basename(phojetExecutableFileName);

  // at the beginning of buffer loop (once every "buffer_size" events)
  if( (e.id().event() - 1) % buffer_size == 0 ){
    // create phojet configuration in tmp dir
    path tmpCfgFile = tmpDirectoryPath / "phojet_cfg.inp";
    phojetStatus = generate_input_card( tmpCfgFile, cms_energy, buffer_size, phojet_process_description );

    if( phojetStatus == 0){
      // remove output ASCII file (if exists)
      remove(outputFilePath);

      // run phojet executable to generate file with "buffer_size" events
      // initialization might be time consuming
      if( e.id().event() == 1){
        run_phojet( tmpDirectoryPath, tmpExeFile, tmpCfgFile, true ); // first event, print cross section
      } else {
        run_phojet( tmpDirectoryPath, tmpExeFile, tmpCfgFile, false );
      }

      // initialize ASCII file reader
      reader_->initialize(outputFilePath.string());
    } else {
      return;
    }
  }

  if( phojetStatus == 0 ){
    // read one event from generated ASCII file (frederik.txt)
    gEv = reader_->fillCurrentEventData();
  }

  // store generator event to the FW event
  auto_ptr<HepMCProduct> output(new HepMCProduct()) ;
  output->addHepMCData(gEv);
  e.put(output);
}

//----------------------------------------------------------------------------------------------------
void Phojet::endJob()
{
  if (verbosity > 10 ){
    cout << "Phojet::endJob()" << endl;
    show_files( tmpDirectoryPath );
  }

  // 1. delete tmp directory and all files inside
  remove_all(tmpDirectoryPath);

  if (verbosity > 10 )
    cout << "Directory " << tmpDirectoryPath << " removed" << endl;
}


DEFINE_FWK_MODULE(Phojet);

