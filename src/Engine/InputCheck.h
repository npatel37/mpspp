/*
Copyright (c) 2012, UT-Battelle, LLC
All rights reserved

[MPS++, Version 0.1]
[by K. Al-Hassanieh, Oak Ridge National Laboratory]
[by J. Rincon, Oak Ridge National Laboratory]
[by G.A., Oak Ridge National Laboratory]

See full open source LICENSE under file LICENSE
in the root directory of this distribution.

*********************************************************

DISCLAIMER

UT-Battelle, LLC AND THE GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL
WARRANTIES, BOTH EXPRESSED AND IMPLIED.  THERE ARE NO EXPRESS OR IMPLIED
WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE
USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT, COPYRIGHT, TRADEMARK, OR
OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL ACCOMPLISH THE INTENDED
RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE.
THE USER ASSUMES RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS,
CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR ARISING
OUT OF, IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.

*********************************************************

*/
/** \ingroup MPSPP */
/*@{*/

/*! \file InputCheck.h
 *
 *  InputChecking functions
 */
#ifndef INPUT_CHECK_H
#define INPUT_CHECK_H
#include <vector>
#include <string>
#include <stdexcept>
#include "TypeToString.h"
#include "Options.h"

namespace Mpspp {

class InputCheck {

	typedef PsimagLite::Options::Readable OptionsReadableType;

public:

	InputCheck() : optsReadable_(0) {}

	~InputCheck()
	{
		if (optsReadable_!=0) delete optsReadable_;
	}

	bool check(const std::string& label,const std::vector<std::string>& vec,size_t line) const
	{
		//			if (label=="JMVALUES") {
		//				if (vec.size()!=2) return error1("JMVALUES",line);
		//				return true;
		//			} else if (label=="RAW_MATRIX") {
		//				size_t row = atoi(vec[0].c_str());
		//				size_t col = atoi(vec[1].c_str());
		//				size_t n = row*col;
		//				if (vec.size()!=n+2) return error1("RAW_MATRIX",line);
		//				return true;
		//			} else if (label=="Connectors") {
		//				return true;
		//			} else if (label=="MagneticField") {
		//				return true;
		//			} else if (label=="FiniteLoops") {
		//				size_t n = atoi(vec[0].c_str());
		//				if (vec.size()!=3*n+1)  return error1("FiniteLoops",line);
		//				return true;
		//			}
		return false;
	}

	void check(const std::string& label,const std::string& val,size_t line)
	{
		if (label!="SolverOptions") return;
		std::vector<std::string> registerOpts;

		//			registerOpts.push_back("restart");
		//			registerOpts.push_back("debugmatrix");
		//			registerOpts.push_back("test");
		//			registerOpts.push_back("useDavidson");
		//			registerOpts.push_back("verbose");
		//			registerOpts.push_back("nofiniteloops");
		//			registerOpts.push_back("nowft");
		//			registerOpts.push_back("inflate");
		//			registerOpts.push_back("none");
		//			registerOpts.push_back("ChebyshevSolver");
		//			registerOpts.push_back("InternalProductStored");
		//			registerOpts.push_back("InternalProductKron");
		//			registerOpts.push_back("useSu2Symmetry");
		//			registerOpts.push_back("TimeStepTargetting");
		//			registerOpts.push_back("DynamicTargetting");
		//			registerOpts.push_back("AdaptiveDynamicTargetting");
		//			registerOpts.push_back("CorrectionVectorTargetting");
		//			registerOpts.push_back("CorrectionTargetting");
		//			registerOpts.push_back("MettsTargetting");

		PsimagLite::Options::Writeable optWriteable(registerOpts,PsimagLite::Options::Writeable::DISABLED);
		optsReadable_ = new  OptionsReadableType(optWriteable,val);
	}

	bool isSet(const std::string& thisOption) const
	{
		return optsReadable_->isSet(thisOption);
	}

	//		void checkForThreads(size_t nthreads) const
	//		{
	//			if (nthreads==1) return;

	//			std::string message1(__FILE__);
	//			message1 += " FATAL: You are requesting nthreads>0 but you did not compile with USE_PTHREADS enabled\n";
	//			message1 += " Either set Threads=1 in the input file (you won't have threads though) or\n";
	//			message1 += " add -DUSE_PTHREADS to the CPP_FLAGS in your Makefile and recompile\n";
	//			throw std::runtime_error(message1.c_str());
	//		}

	void usageMain(const std::string& name) const
	{
		std::cerr<<"USAGE is "<<name<<"\n";
	}

private:

	bool error1(const std::string& message,size_t line) const
	{
		std::string s(__FILE__);
		s += " : Input error for label " + message + " near line " + ttos(line) + "\n";
		throw std::runtime_error(s.c_str());

	}

	OptionsReadableType* optsReadable_;

}; // class InputCheck
} // namespace Mpspp

/*@}*/
#endif

