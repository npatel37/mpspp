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

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************

*/
/** \ingroup MPSPP */
/*@{*/

#ifndef CONTRACTED_PART_H
#define CONTRACTED_PART_H

#include "ProgramGlobals.h"
#include "ContractedFactor.h"

namespace Mpspp {

template<typename MatrixProductOperatorType>
class ContractedPart {

	typedef typename MatrixProductOperatorType::MatrixProductStateType MatrixProductStateType;
	typedef typename MatrixProductOperatorType::MpoFactorType MpoFactorType;
	typedef typename MatrixProductStateType::ComplexOrRealType ComplexOrRealType;

	enum {TO_THE_RIGHT = ProgramGlobals::TO_THE_RIGHT, TO_THE_LEFT = ProgramGlobals::TO_THE_LEFT};

	enum {PART_LEFT = ProgramGlobals::PART_LEFT, PART_RIGHT = ProgramGlobals::PART_RIGHT};

public:

	typedef typename ProgramGlobals::CrsMatrix<ComplexOrRealType>::Type SparseMatrixType;
	typedef ContractedFactor<MatrixProductOperatorType> ContractedFactorType;

	ContractedPart(const MatrixProductStateType& abState,const MatrixProductOperatorType& h)
		: abState_(abState),h_(h)
	  //: dataLeft_(abState.center()+1),dataRight_(abState.sites()-dataLeft_.size()+1)
	{
		// page 62, equations 192 and 193
//		assert(dataLeft_.size()>0);
//		dataLeft_[0].init(abState(0),h(0),0,PART_LEFT,0);
//		for (size_t i=1;i<dataLeft_.size();i++) {
//			dataLeft_[i].init(abState(i-1),h(i-1),i,PART_LEFT,&dataLeft_[i-1]);
//			std::cerr<<"Testing: ContracedPart (left) i="<<i<<" out of "<<dataLeft_.size()<<"\n";
//		}
//		size_t center = abState.center();
//		std::cerr<<"ContractedPart center="<<center<<"\n";
//		dataRight_[center].initRight(abState(center),h(center),center,0);
//		for (size_t i=1;i<dataRight_.size();i++) {
//			dataRight_[i].init(abState(i+center-1),h(i+center-1),i+center,PART_RIGHT,&dataRight_[i-1]);
//			std::cerr<<"Testing: ContracedPart (right) i="<<i<<" out of "<<dataRight_.size()<<"\n";
//		}
	}

	void growRight(size_t currentSite)
	{
		ContractedFactorType cf(abState_.B(currentSite),h_(currentSite),currentSite,ProgramGlobals::PART_RIGHT,0);
		R_.push_back(cf);
//		if (L_.size()==0) {
//			ContractedFactorType cfL;
//			L_.push_back(cfL);
//		}
	}

	//! From As (or Bs) and Ws reconstruct *this
	void update(size_t currentSite,const MatrixProductStateType& abState,size_t direction)
	{
		if (direction==TO_THE_RIGHT) {
			updateLeft(currentSite,abState);
		} else {
			updateRight(currentSite,abState);
		}
	}

	const ContractedFactorType& operator()(size_t currentSite,size_t leftOrRight) const
	{
		if (leftOrRight == PART_LEFT) {
			assert(currentSite<L_.size());
			std::cout<<L_;
		} else {
			assert(currentSite<R_.size());
		}

		return (leftOrRight == PART_LEFT) ? L_[currentSite] : R_[currentSite];
	}

	template<typename MatrixProductOperatorType_>
	friend std::ostream& operator<<(std::ostream& os,const ContractedPart<MatrixProductOperatorType_>& contractedPart);

private:

	void updateLeft(size_t currentSite,const MatrixProductStateType& abState)
	{
		if (currentSite>=L_.size()) {
			ContractedFactorType* dataPrev = 0;
			if (L_.size()>0) dataPrev = &(L_[L_.size()-1]);
			ContractedFactorType cf(abState.A(currentSite),h_(currentSite),currentSite,ProgramGlobals::PART_LEFT,dataPrev);
			L_.push_back(cf);
			return;
		}
		L_[currentSite].update(abState.A(currentSite));
	}

	void updateRight(size_t currentSite,const MatrixProductStateType& abState)
	{
		std::string str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Need to updateRight(...) here. I cannot go further until this is implemented\n";
		throw std::runtime_error(str.c_str());
	}

	const MatrixProductStateType& abState_;
	const MatrixProductOperatorType& h_;
	typename ProgramGlobals::Vector<ContractedFactorType>::Type R_;
	typename ProgramGlobals::Vector<ContractedFactorType>::Type L_;

}; // ContractedPart

template<typename MatrixProductOperatorType>
std::ostream& operator<<(std::ostream& os,const ContractedPart<MatrixProductOperatorType>& contractedPart)
{
	os<<"ContractedPart: right size="<<contractedPart.R_.size()<<"\n";
	for (size_t i=0;i<contractedPart.R_.size();i++)
		os<<contractedPart.R_[i];
	os<<"ContractedPart: left size="<<contractedPart.L_.size()<<"\n";
	for (size_t i=0;i<contractedPart.L_.size();i++)
		os<<contractedPart.L_[i];
	return os;
}

} // namespace Mpspp

/*@}*/
#endif // CONTRACTED_PART_H

