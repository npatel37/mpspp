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

namespace Mpspp {

template<typename MatrixProductOperatorType>
class LeftRightSuper {

public:

	typedef ContractedLeftPart<MatrixProductOperatorType> ContractedLeftPartType;
	typedef ContractedRightPart<MatrixProductOperatorType> ContractedRightPartType;

	LeftRightSuper(MatrixProductState& A,
				   MatrixProductState& B,
				   ContractedLeftPartType& cL,
				   ContractedRightPartType& cR)
	{}

	//! Moves the center of orthogonality by one to the right
	void moveRight() const
	{
		std::string str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Need to implement moveRight here. I cannot go further until this is implemented\n";
		throw std::runtime_error(str.c_str());
	}


	//! Moves the center of orthogonality by one to the left
	void moveLeft() const
	{
		std::string str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += "Need to implement moveLeft here. I cannot go further until this is implemented\n";
		throw std::runtime_error(str.c_str());
	}

	void printReport(std::ostream& os) const
	{
		os<<"Nothing to report so far, except that I need a progress indicator\n";
	}

}; // LeftRightSuper

} // namespace Mpspp

/*@}*/
#endif // LEFT_RIGHT_SUPER_H

