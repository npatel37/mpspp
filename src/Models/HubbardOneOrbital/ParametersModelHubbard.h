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

/*! \file ParametersModelHubbard.h
 *
 *  Contains the parameters for the Hubbard model and function to read them from a file
 *
 */
#ifndef PARAMETERSMODELHUBBARD_H
#define PARAMETERSMODELHUBBARD_H

namespace Mpspp {
//! Hubbard Model Parameters
template<typename Field>
struct ParametersModelHubbard {

	typedef typename PsimagLite::Vector<Field>::Type VectorType;

	template<typename IoInputType>
	ParametersModelHubbard(IoInputType& io)
	{
		io.read(hubbardU,"hubbardU");
		io.read(potentialV,"potentialV");

		try {
			io.read(potentialT,"PotentialT"); //level,beQuiet);
		} catch (std::exception& e) {}
		omega=0;
		try {
			io.readline(omega,"omega=");
		} catch (std::exception& e) {}
	}

	// Do not include here connection parameters
	// those are handled by the Geometry
	// Hubbard U values (one for each site)
	VectorType hubbardU;
	// Onsite potential values, one for each site
	VectorType potentialV;

	// for time-dependent H:
	VectorType potentialT;
	Field omega;

	// target number of electrons  in the system
	int nOfElectrons;
};

//! Function that prints model parameters to stream os
template<typename FieldType>
std::ostream& operator<<(std::ostream &os,
                         const ParametersModelHubbard<FieldType>& parameters)
{
	os<<"hubbardU\n";
	os<<parameters.hubbardU;
	os<<"potentialV\n";
	os<<parameters.potentialV;
	if (parameters.potentialT.size()==0) return os;

	// time-dependent stuff
	os<<"potentialT\n";
	os<<parameters.potentialT;
	os<<"omega="<<parameters.omega<<"\n";
	return os;
}
} // namespace Mpspp

/*@}*/
#endif
