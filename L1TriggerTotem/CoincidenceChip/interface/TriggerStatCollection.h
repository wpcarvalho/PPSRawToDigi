// -*- C++ -*-
//
// Package:    L1TriggerTotem
// Class:      RPTriggerAnalyzer
//
// Original Author:  Jiri Prochazka
//         Created:  Mon Mar  1 15:34:46 CET 2010
// $Id$

#ifndef _L1TriggerTotemTriggerStatCollection_H_
#define _L1TriggerTotemTriggerStatCollection_H_

// system
#include <vector>
#include <iostream>
#include <cassert>

// ROOT
#include "TObject.h"
#include "TString.h"

template <class T> 
// class TriggerStatCollection: public TObject, public edm::OwnVector<T>{
// class TriggerStatCollection: public std::vector<T*>{
class TriggerStatCollection: public TObject, public std::vector<T*>{
// class TriggerStatCollection: public boost::ptr_vector<T>{
    public:
        TriggerStatCollection(){}
        ~TriggerStatCollection(){}
        void AddIfNotExists(const TString& name){
            if(Find(name)){
                //  cout << "ERROR: Trigger stats \"" + name+"\" already exists" << endl;
                //  assert(!FindTriggerStats(name));
            }else{
               this->push_back(new T(name));
            }
        }
        void AddIfNotExists( T* object){
            if(Find(object->GetName())){
                //  cout << "ERROR: Trigger stats \"" + name+"\" already exists" << endl;
                //  assert(!FindTriggerStats(name));
            }else{
               this->push_back(object);
            }
        }

#if 1        
        T* Find(const TString& iname) {
            for(unsigned int i=0; i < this->size(); ++i){
                // if(!strcmp((*this)[i].GetName(),iname)) return &(*this)[i];
                if(!strcmp((*this)[i]->GetName(),iname)) return (*this)[i];
            }
            //edm::LogInfo("CClogic") <<  "inputBits" << nevents;
            // cout << "Can not find TriggerStat for condition: \"" << name << "\"" << endl;
            //throw cms::Exception("PPBckgAnalyzer2::FindHist2") << "Can not find histogram named " << name << endl;
            //  assert(0);
            return 0;
        }
        T* AssertFind(const TString& iname) {
            T* t = Find(iname);
            if(!t){
            //edm::LogInfo("CClogic") <<  "inputBits" << nevents;
                std::cout << "Can not find TriggerStat for condition: \"" << iname << "\"" << std::endl;
            //throw cms::Exception("PPBckgAnalyzer2::FindHist2") << "Can not find histogram named " << name << endl;
              assert(0);
            }
            return t;
        }
       const T* AssertFind(const TString& iname) const{
            const T* t = Find(iname);
            if(!t){
            //edm::LogInfo("CClogic") <<  "inputBits" << nevents;
                std::cout << "Can not find TriggerStat for condition: \"" << iname << "\"" << std::endl;
            //throw cms::Exception("PPBckgAnalyzer2::FindHist2") << "Can not find histogram named " << name << endl;
              assert(0);
            }
            return t;
        }


        const T* Find(const TString& iname) const {
            for(unsigned int i=0; i < this->size(); ++i){
                // if(!strcmp((*this)[i].GetName(),iname)) return &(*this)[i];
                if(!strcmp((*this)[i]->GetName(),iname)) return (*this)[i];
            }
            //edm::LogInfo("CClogic") <<  "inputBits" << nevents;
            // cout << "Can not find TriggerStat for condition: \"" << name << "\"" << endl;
            //throw cms::Exception("PPBckgAnalyzer2::FindHist2") << "Can not find histogram named " << name << endl;
            //  assert(0);
            return 0;
        }
#endif  
      //  typedef TriggerStatCollection<T> tscol;
        ClassDef(TriggerStatCollection,1)
};
#endif
