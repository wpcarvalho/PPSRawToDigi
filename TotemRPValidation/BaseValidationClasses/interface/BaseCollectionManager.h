#ifndef TotemRPValidation_BaseValidationClasses_BaseCollectionManager_h
#define TotemRPValidation_BaseValidationClasses_BaseCollectionManager_h

#include <map>
#include <boost/shared_ptr.hpp>
#include "TFile.h"
#include <sstream>
#include <iostream>


template <class DEB, class ID, class CONF>
class BaseCollectionManager
{
  public:
    typedef std::map<ID, boost::shared_ptr<DEB> > container;
    typedef class container::iterator container_it;

    BaseCollectionManager(std::string base_path, const CONF& conf) : conf_(conf) {base_path_=base_path;}

    boost::shared_ptr<DEB> GetObj(ID Id)
    {
      container_it it = coll_cont.find(Id);
      if(it!=coll_cont.end())
      {
        return it->second;
      }
      else 
      {
        //std::ostringstream osstr;
        //osstr<<Id<<'/';
        char temp[128];
        sprintf(temp, "%04d", Id);
        boost::shared_ptr<DEB> ptr = boost::shared_ptr<DEB>(new DEB(base_path_ + temp + "/", Id, conf_));
        coll_cont.insert(std::pair<ID, boost::shared_ptr<DEB> >(Id, ptr) );
        return ptr;
      }
    }
    void Write(TFile *f)
    {
      if(f)
      {
        for(container_it it=coll_cont.begin(); it!=coll_cont.end(); ++it)
        {
          //std::cout<<"object id="<<it->first<<" written to the file"<<std::endl;
          it->second->WriteRootFile(f);
        }
      }
    }
    
  private:
    container coll_cont;
    std::string base_path_;
    const CONF conf_;
};

#endif
