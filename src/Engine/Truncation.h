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

#ifndef TRUNCATION_H
#define TRUNCATION_H

#include "Sort.h"
#include "ProgramGlobals.h"

namespace Mpspp {

template<typename ContractedLocalType>
class Truncation {
	
	typedef typename ContractedLocalType::MatrixProductOperatorType MpoLocalType;
	typedef typename MpoLocalType::MpsLocalType MpsLocalType;
	typedef typename MpsLocalType::VectorRealType VectorRealType;
	typedef typename VectorRealType::value_type RealType;
	typedef typename ProgramGlobals::Vector<size_t>::Type VectorIntegerType;
	typedef typename ContractedLocalType::SymmetryLocalType SymmetryLocalType;
	typedef typename SymmetryLocalType::SymmetryFactorType SymmetryFactorType;

public:

	Truncation(MpsLocalType& mps,ContractedLocalType& contracted)
		: mps_(mps),contracted_(contracted)
	{}

//	size_t size() const { return s_.size(); }

//	RealType& operator()(size_t i) { return s_[i]; }

	void operator()(SymmetryLocalType& symm,size_t site,size_t part,size_t cutoff)
	{
		if (part==ProgramGlobals::PART_LEFT) {
			if (symm(site).left().size()<=cutoff) return;
		} else {
			if (symm(site).right().size()<=cutoff) return;
		}
		mps_.truncate(site,part,cutoff);
		contracted_.truncate(site,part,cutoff);
		symm.truncate(site,part,cutoff);
	}

//	void order() const
//	{
//		Sort<VectorRealType> sort;
//		sort.sort(s_,perm_);
//	}

//	void print(std::ostream& os) const
//	{
//		os<<"s.size= "<<s_.size()<<"\n";
//		for (size_t i=0;i<s_.size();i++)
//			os<<s_[i]<<" ";
//		os<<"\n";
//	}

//	void set(const VectorRealType& s) { s_ = s; }

//	void truncate(size_t spaceSize)
//	{
//		std::cout<<"spaceSize= "<<spaceSize<<" truncate()= "<<s_.size()<<"   ";
//		for (size_t i=0;i<s_.size();i++) {
//			std::cout<<s_[i]<<" ";
//		}
//		std::cout<<"\n";
//	}

//	template<typename SomeMatrixType>
//	void recoverSvd(SomeMatrixType& mat,const SomeMatrixType& u,const SomeMatrixType& vt) const
//	{
//		size_t m = mat.n_row();
//		size_t n = vt.n_col();
//		size_t min = s_.size();
//		assert(u.n_col()>=min);
//		assert(vt.n_row()>=min);
//		for (size_t i=0;i<m;i++) {
//			for (size_t j=0;j<n;j++) {
//				mat(i,j) = 0.0;
//				for (size_t k=0;k<min;k++)
//					mat(i,j) += u(i,k) * s_[k] * vt(k,j);
//			}
//		}
//	}

private:

	MpsLocalType& mps_;
	ContractedLocalType& contracted_;

//	VectorRealType s_;
//	VectorIntegerType perm_;

}; // Truncation

} // namespace Mpspp

/*@}*/
#endif // TRUNCATION_H

