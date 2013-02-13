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

#ifndef MODEL_SELECTOR_H
#define MODEL_SELECTOR_H

#include <string>
#include <stdexcept>
#include "HubbardOneOrbital.h"
#include "HeisenbergSpinOneHalf.h"

namespace Mpspp {

template<typename ModelType>
class ModelSelector {

	typedef typename ModelType::ParametersSolverType ParametersSolverType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename ModelType::SymmetryLocalType SymmetryLocalType;
	typedef typename ModelType::GeometryType GeometryType;
	typedef typename ModelType::ConcurrencyType ConcurrencyType;
	typedef typename ModelType::MpoLocalType MpoLocalType;

	typedef HubbardOneOrbital<ParametersSolverType,
							  InputValidatorType,
							  SymmetryLocalType,
							  GeometryType,
							  ConcurrencyType> HubbardOneOrbitalType;

	typedef HeisenbergSpinOneHalf<ParametersSolverType,
							  InputValidatorType,
							  SymmetryLocalType,
							  GeometryType,
							  ConcurrencyType> HeisenbergSpinOneHalfType;

public:

	typedef ModelBase<ParametersSolverType,
					  InputValidatorType,
					  SymmetryLocalType,
					  GeometryType,
					  ConcurrencyType> ModelBaseType;

	ModelSelector(const std::string& name)
	: name_(name),model_(0)
	{}

	~ModelSelector()
	{
		if (model_) delete model_;
	}

	const ModelBaseType& operator()(const ParametersSolverType& solverParams,
						  InputValidatorType& io,
						  const GeometryType& geometry,
						  ConcurrencyType& concurrency)
	{
		if (name_ == "HubbardOneBand") {
			model_ = new HubbardOneOrbitalType(solverParams,io,geometry,concurrency);
		} else if (name_ == "HeisenbergSpinOneHalf") {
			model_ = new HeisenbergSpinOneHalfType(solverParams,io,geometry,concurrency);
		} else {
			std::string str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "Unknown model " + name_ + "\n";
			throw std::runtime_error(str.c_str());
		}
		return *model_;
	}

private:

	std::string name_;
	ModelBaseType* model_;

}; // ModelSelector

} // namespace Mpspp

/*@}*/
#endif // MODEL_SELECTOR_H

