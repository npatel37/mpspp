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
	  : dataLeft_(abState.center()),dataRight_(abState.sites()-dataLeft_.size()-1)
	{
		// page 62, equations 192 and 193
		for (size_t i=0;i<dataLeft_.size();i++) {
			dataLeft_[i].init(abState(i),h(i),i,PART_LEFT,(i>0) ? &dataLeft_[i-1] : 0);
			std::cerr<<"Testing: ContracedPart (left) i="<<i<<" out of "<<dataLeft_.size()<<"\n";
		}
		size_t center = dataLeft_.size();
		for (size_t i=0;i<dataRight_.size();i++) {
			dataRight_[i].init(abState(i+center),h(i+center),i+center,PART_RIGHT,(i>0) ? &dataRight_[i-1] : 0);
			std::cerr<<"Testing: ContracedPart (right) i="<<i<<" out of "<<dataRight_.size()<<"\n";
		}
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
		return (leftOrRight == PART_LEFT) ? dataLeft_[currentSite] : dataRight_[currentSite];
	}

private:

	void updateLeft(size_t currentSite,const MatrixProductStateType& abState)
	{
		std::string str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Need to update(...) here. I cannot go further until this is implemented\n";
		throw std::runtime_error(str.c_str());
	}

	void updateRight(size_t currentSite,const MatrixProductStateType& abState)
	{
		std::string str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Need to update(...) here. I cannot go further until this is implemented\n";
		throw std::runtime_error(str.c_str());
	}

	typename ProgramGlobals::Vector<ContractedFactorType>::Type dataLeft_;
	typename ProgramGlobals::Vector<ContractedFactorType>::Type dataRight_;

}; // ContractedPart

} // namespace Mpspp

/*@}*/
#endif // CONTRACTED_PART_H

