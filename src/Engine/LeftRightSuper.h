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

#include "ContractedLeftPart.h"
#include "ContractedRightPart.h"
#include "ProgressIndicator.h"

namespace Mpspp {

template<typename MatrixProductOperatorType>
class LeftRightSuper {

public:

	typedef ContractedLeftPart<MatrixProductOperatorType> ContractedLeftPartType;
	typedef ContractedRightPart<MatrixProductOperatorType> ContractedRightPartType;

	LeftRightSuper(MatrixProductState& A,
	               ContractedLeftPartType& cL,
	               MatrixProductState& B,
	               ContractedRightPartType& cR)
	: progress_("LeftRightSuper",0),
	  A_(A),
	  cL_(cL),
	  B_(B),
	  cR_(cR)
	{}

	//! Moves the center of orthogonality by one to the right
	void moveRight()
	{
		updateA();
		cL_.update(A_);
	}

	//! Moves the center of orthogonality by one to the left
	void moveLeft()
	{
		updateB();
		cR_.update(B_);
	}

	void printReport(std::ostream& os) const
	{
		os<<"Nothing to report so far, except that I need a progress indicator\n";
	}

private:

	//! From cL and cR construct a new A, only A changes here
	void updateA()
	{
		std::string str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Need to updateA(...) here. I cannot go further until this is implemented\n";
		throw std::runtime_error(str.c_str());
	}

	//! From cL and cR construct a new B, only B changes here
	void updateB()
	{
		std::string str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Need to updateB(...) here. I cannot go further until this is implemented\n";
		throw std::runtime_error(str.c_str());
	}

	PsimagLite::ProgressIndicator progress_;
	MatrixProductState& A_;
	ContractedLeftPartType& cL_;
	MatrixProductState& B_;
	ContractedRightPartType& cR_;
}; // LeftRightSuper

} // namespace Mpspp

/*@}*/
#endif // LEFT_RIGHT_SUPER_H

