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

#ifndef LEFT_RIGHT_SUPER_H
#define LEFT_RIGHT_SUPER_H

#include "ContractedPart.h"
#include "ProgressIndicator.h"
#include <vector>

namespace Mpspp {

template<typename MatrixProductOperatorType_,typename VectorType,typename RealType_>
class LeftRightSuper {

	enum {TO_THE_RIGHT = ProgramGlobals::TO_THE_RIGHT, TO_THE_LEFT = ProgramGlobals::TO_THE_LEFT};

public:

	typedef MatrixProductOperatorType_ MatrixProductOperatorType;
	typedef typename MatrixProductOperatorType::MatrixProductStateType MatrixProductStateType;
	typedef typename MatrixProductStateType::SymmetryLocalType SymmetryLocalType;
	typedef typename MatrixProductStateType::ComplexOrRealType ComplexOrRealType;
	typedef RealType_ RealType;
	typedef ContractedPart<MatrixProductOperatorType> ContractedPartType;

	LeftRightSuper(MatrixProductStateType& ab,
				   ContractedPartType& contracted)
	: progress_("LeftRightSuper",0),
	  ab_(ab),
	  contracted_(contracted)
	{}

	const MatrixProductStateType& abState() const { return ab_; }

	const ContractedPartType& contracted() const { return contracted_; }

	const SymmetryLocalType& symmetry() const { return ab_.symmetry(); }

	void updateContracted(size_t currentSite,const MatrixProductStateType& abState,size_t direction)
	{
		contracted_.update(currentSite,abState,direction);
	}

	void updateMps(size_t currentSite,const VectorType& v,size_t direction)
	{
		ab_.update(currentSite,v,direction);
	}

private:

	PsimagLite::ProgressIndicator progress_;
	MatrixProductStateType& ab_;
	ContractedPartType& contracted_;
}; // LeftRightSuper

} // namespace Mpspp

/*@}*/
#endif // LEFT_RIGHT_SUPER_H

