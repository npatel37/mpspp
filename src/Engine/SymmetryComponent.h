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

#ifndef SYMMETRY_COMPONENT_H
#define SYMMETRY_COMPONENT_H

#include "ProgramGlobals.h"
#include "IoSimple.h"
#include "Sort.h"

namespace Mpspp {

class SymmetryComponent {

public:

	typedef PsimagLite::IoSimple::In IoInputType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef PsimagLite::Vector<SizeType>::Type VectorIntegerType;

	enum {CORNER_LEFT,CORNER_RIGHT};
	enum {COMPONENT_LEFT,COMPONENT_RIGHT,COMPONENT_SUPER};

	SymmetryComponent(SizeType type)
	    : type_(type),leftSize_(0)
	{}

	SymmetryComponent(SizeType type,
	                  SizeType hilbert,
	                  SizeType site,
	                  const VectorIntegerType& quantumNumbers)
	    : type_(type),
	      leftSize_(hilbert),
	      block_(1,site),
	      quantumNumbers_(quantumNumbers)
	{
		if (quantumNumbers.size()==1) block_.resize(0);

		findPermutationAndPartition();
	}

	void setSite(SizeType site)
	{
		assert(block_.size()==1);
		block_[0] = site;
	}

	void setComponent(SizeType comp)
	{
		type_=comp;
	}

	void grow(SizeType site,
	          const VectorIntegerType& quantumNumbers,
	          SizeType nsites,
	          SizeType leftSize)
	{
		assert(type_==COMPONENT_RIGHT);
		SymmetryComponent sc1(type_,0,site,quantumNumbers);
		SymmetryComponent sc(type_);
		*this = sc1;
		SizeType site1 = site+1;
		for (SizeType i=site1;i<nsites;i++) {
			SymmetryComponent sc2(type_,0,i,quantumNumbers);
			sc.combine(*this,sc2);
			*this = sc;
		}
		this->leftSize_ = leftSize;
	}

	void combine(const SymmetryComponent& left,const SymmetryComponent& right)
	{
		//block_.clear();
		blockUnion(block_,left.block_,right.block_); //! B= pS.block Union X
		assert(isValidBlock());

		SizeType ns = left.size();
		SizeType ne = right.size();

		quantumNumbers_.clear();

		for (SizeType j=0;j<ne;j++) {
			for (SizeType i=0;i<ns;i++) {
				quantumNumbers_.push_back(left.quantumNumbers_[i]+right.quantumNumbers_[j]);
			}
		}

		if (ne==0) quantumNumbers_= left.quantumNumbers_;
		if (ns==0) quantumNumbers_= right.quantumNumbers_;

		// order quantum numbers of combined basis:
		findPermutationAndPartition();
		leftSize_ = left.size();
	}

	template<typename SomeTruncationType>
	void truncate(SizeType cutoff,const SomeTruncationType& trunc)
	{
		if (size()<=cutoff) return;
		trunc.vector(quantumNumbers_,cutoff);
		findPermutationAndPartition();
	}

	PsimagLite::String typeToString() const
	{
		switch (type_) {
		case COMPONENT_LEFT:
			return "left";
		case COMPONENT_RIGHT:
			return "right";
		case COMPONENT_SUPER:
			return "super";
		}
		return "unknown";
	}

	SizeType partitions() const { return partition_.size(); }

	SizeType partitionSize(SizeType i) const
	{
		assert(i+1<partition_.size());
		return partition_[i+1]-partition_[i];
	}

	SizeType partitionOffset(SizeType i) const
	{
		assert(i<partition_.size());
		return partition_[i];
	}

	PairType unpack(SizeType i) const
	{
		assert(i<permutation_.size());
		SizeType ip = permutation_[i];
		if (leftSize_==0 && type_==COMPONENT_LEFT)
			return PairType(0,ip);
		if (leftSize_==0 && type_==COMPONENT_RIGHT)
			return PairType(ip,0);
		div_t q = div(ip,leftSize_);
		return PairType(q.rem,q.quot);
	}

	SizeType pack(SizeType a1,SizeType a2) const
	{
		assert(a1+a2*leftSize_<permutationInverse_.size());
		if (leftSize_==0 && type_==COMPONENT_LEFT) {
			assert(a1==0);
			return permutationInverse_[a2];
		}
		if (leftSize_==0 && type_==COMPONENT_RIGHT) {
			assert(a2==0);
			return permutationInverse_[a1];
		}
		return permutationInverse_[a1+a2*leftSize_];
	}

	SizeType size() const
	{
		return quantumNumbers_.size();
	}

	SizeType qn(SizeType state) const
	{
		assert(state<quantumNumbers_.size());
		return quantumNumbers_[state];
	}

	SizeType split() const { return leftSize_; }

	const VectorIntegerType& block() const { return block_; }

	const bool operator==(const SymmetryComponent& other) const
	{
		bool b1 = (type_ == other.type_);
		bool b2 = (leftSize_ == other.leftSize_);
		bool b3 = (block_ == other.block_);
		bool b4 = (quantumNumbers_ == other.quantumNumbers_);
		bool b5 = (partition_ == other.partition_);
		bool b6 = (permutation_ == other.permutation_);
		bool b7 = (permutationInverse_ == other.permutationInverse_);
		return (b1 && b2 && b3 && b4 && b5 && b6 && b7);
	}

	friend std::ostream& operator<<(std::ostream& os,
	                                const SymmetryComponent& symm);

private:

	// Match with SaveInternal in DMRG++'s Basis.h
	template<typename IoInputter>
	void loadInternal(IoInputter& io)
	{
		//int x=0;
		//useSu2Symmetry_=false;
		//io.readline(x,"#useSu2Symmetry=");
		//if (x>0) useSu2Symmetry_=true;
		io.read(block_,"#BLOCK");
		io.read(quantumNumbers_,"#QN");
		//io.read(electrons_,"#ELECTRONS");
		//io.read(electronsOld_,"#0OLDELECTRONS");
		io.read(partition_,"#PARTITION");
		io.read(permutationInverse_,"#PERMUTATIONINVERSE");
		permutation_.resize(permutationInverse_.size());
		for (SizeType i=0;i<permutation_.size();i++)
			permutation_[permutationInverse_[i]]=i;
		/*dmrgTransformed_=false;
		if (useSu2Symmetry_)
			symmSu2_.load(io);
		else
			symmLocal_.load(io);*/
	}

	void findPermutationAndPartition()
	{

		permutation_.resize(size());

		PsimagLite::Sort<VectorIntegerType> sort;
		sort.sort(quantumNumbers_,permutation_);

		findPartition();

		permutationInverse_.resize(permutation_.size());
		for (SizeType i=0;i<permutationInverse_.size();i++)
			permutationInverse_[permutation_[i]]=i;
	}

	//! Finds a partition of the basis given the effecitve quantum numbers
	//! Find a partition of the basis given the effecitve quantum numbers
	void findPartition()
	{
		assert(quantumNumbers_.size()>0);

		SizeType qtmp = quantumNumbers_[0]+1;
		partition_.clear();
		for (SizeType i=0;i<size();i++) {
			if (quantumNumbers_[i]!=qtmp) {
				partition_.push_back(i);
				qtmp = quantumNumbers_[i];
			}
		}

		partition_.push_back(size());
	}

	//! A = B union C
	template<typename Block>
	void blockUnion(Block &A,Block const &B,Block const &C)
	{
		A=B;
		for (SizeType i=0;i<C.size();i++) A.push_back(C[i]);
	}

	bool isValidBlock() const
	{
		if (block_.size()<2) return true;
		for (SizeType i=0;i<block_.size()-1;i++) {
			if (block_[i]>=block_[i+1])
				return false;
		}
		return true;
	}

	SizeType type_;
	SizeType leftSize_;
	VectorIntegerType block_;
	VectorIntegerType quantumNumbers_;
	VectorIntegerType partition_;
	VectorIntegerType permutation_;
	VectorIntegerType permutationInverse_;
}; // SymmetryComponent

std::ostream& operator<<(std::ostream& os,const SymmetryComponent& symm)
{
	os<<"type="<<symm.typeToString()<<" leftSize_= ";
	os<<symm.leftSize_<<" size= "<<symm.quantumNumbers_.size();
	os<<" block= ";
	for (SizeType i=0;i<symm.block_.size();i++)
		os<<symm.block_[i]<<" ";
	os<<" ";
	return os;
}

} // namespace Mpspp

/*@}*/
#endif // SYMMETRY_COMPONENT_H

