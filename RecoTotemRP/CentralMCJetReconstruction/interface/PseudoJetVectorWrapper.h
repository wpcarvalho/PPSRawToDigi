/*
 * PseudoJetVectorWrapper.h
 *
 *  Created on: 17 Sep 2013
 *      Author: psikora
 */
#ifndef PSEUDOJETVECTORWRAPPER_H_
#define PSEUDOJETVECTORWRAPPER_H_

#include <vector>
#include "fastjet/PseudoJet.hh"

/**
 * This wrapper is needed because dictionary is not compiling
 * with direct reference to fastjet::PseudoJet. (Probably because template constructor.)
 *
 */
class PseudoJetVectorWrapper{
public:
	std::vector<fastjet::PseudoJet> pseudoJetsVector;
};


#endif /* PSEUDOJETVECTORWRAPPER_H_ */
