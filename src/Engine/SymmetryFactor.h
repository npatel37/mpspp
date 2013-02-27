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

#ifndef SYMMETRY_FACTOR_H
#define SYMMETRY_FACTOR_H

#include "SymmetryComponent.h"

namespace Mpspp {

class SymmetryFactor {

public:

	typedef SymmetryComponent SymmetryComponentType;
	typedef typename SymmetryComponentType::IoInputType IoInputType;
	typedef typename SymmetryComponentType::PairType PairType;
	typedef typename SymmetryComponent::VectorType VectorIntegerType;

	enum {CORNER_LEFT = SymmetryComponentType::CORNER_LEFT,
		  CORNER_RIGHT = SymmetryComponentType::CORNER_RIGHT};

//	SymmetryFactor(IoInputType& io,size_t nk)
//		: left_(io,nk),right_(io,nk),super_(io,right_.size())
//	{
//		std::string str("SymmetryFactor: super=");
//		str += ttos(super_.size()) + " left=" + ttos(left_.size());
//		str += " right=" + ttos(right_.size()) + "\n";
//		std::cerr<<str;
//	}

	SymmetryFactor()
		: left_(SymmetryComponentType::COMPONENT_LEFT),
		  right_(SymmetryComponentType::COMPONENT_RIGHT),
		  super_(SymmetryComponentType::COMPONENT_SUPER)
	{}

	void growFirst(size_t site,
				   const std::vector<size_t>& quantumNumbers,
				   size_t nsites)
	{
		assert(site+1<nsites);
		std::vector<size_t> qn(1,0);
		SymmetryComponentType onesiteLeft(SymmetryComponentType::COMPONENT_LEFT,0,site,qn);
		left_ = onesiteLeft;

		size_t siteRight = nsites - 1 -site;
		SymmetryComponentType onesiteRight2(SymmetryComponentType::COMPONENT_RIGHT,0,siteRight,qn);
		SymmetryComponentType onesiteRight3(SymmetryComponentType::COMPONENT_LEFT,0,siteRight,quantumNumbers);
		right_.combine(onesiteRight3,onesiteRight2);

		super_.combine(left_,right_);
	}

	void grow(size_t site,
			  const std::vector<size_t>& quantumNumbers,
			  const SymmetryFactor& previous,
			  size_t nsites)
	{
		assert(site+1<nsites);
		SymmetryComponentType onesiteRight(SymmetryComponentType::COMPONENT_RIGHT,0,site,quantumNumbers);
		left_.combine(previous.left(),onesiteRight);

		size_t siteRight = nsites - 1 -site;
		SymmetryComponentType onesiteRight3(SymmetryComponentType::COMPONENT_LEFT,0,siteRight,quantumNumbers);
		if (site==0) {
			right_=previous.right();
		} else {
			right_.combine(onesiteRight3,previous.right());
		}

		super_.combine(left_,right_);
	}

	void moveRight(const SymmetryComponentType& oldLeft,
				  const SymmetryComponentType& onesite,
				  const SymmetryComponentType& oldRight)
	{
		left_.combine(oldLeft,onesite);
		right_ = oldRight;
		super_.combine(left_,right_);
	}

	void moveLeft(const SymmetryComponentType& oldLeft,
				  const SymmetryComponentType& onesite,
				  const SymmetryComponentType& oldRight)
	{
		left_ = oldLeft;
		right_.combine(onesite,oldRight);
		super_.combine(left_,right_);
	}

	template<typename SomeTruncationType>
	void truncate(size_t part,size_t cutoff,const SomeTruncationType& trunc)
	{
		if (part==ProgramGlobals::PART_LEFT)
			left_.truncate(cutoff,trunc);
		else
			right_.truncate(cutoff,trunc);
		super_.combine(left_,right_);
	}

	void set(size_t part,const SymmetryComponentType& c)
	{
		switch(part) {
		case SymmetryComponent::COMPONENT_LEFT:
			left_=c;
			break;
		case SymmetryComponent::COMPONENT_RIGHT:
			right_=c;
			break;
		case SymmetryComponent::COMPONENT_SUPER:
			super_=c;
			break;
		}
		super_.combine(left_,right_);
	}

	const SymmetryComponentType& super() const { return super_; }

	const SymmetryComponentType& left() const { return left_; }

	const SymmetryComponentType& right() const { return right_; }

	friend std::ostream& operator<<(std::ostream& os,const SymmetryFactor& symm);

private:

	SymmetryComponentType left_;
	SymmetryComponentType right_;
	SymmetryComponentType super_;
}; // SymmetryFactor


std::ostream& operator<<(std::ostream& os,const SymmetryFactor& symm)
{
	os<<symm.left_<<" ";
	os<<symm.right_<<" ";
	os<<symm.super_<<"\n";
	return os;
}

} // namespace Mpspp

/*@}*/
#endif // SYMMETRY_FACTOR_H

