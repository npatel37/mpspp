#include "Provenance.h"

std::ostream& operator<<(std::ostream& os,const Provenance&)
{
	os<<"DMRG++ version "<<MPSPP_VERSION<<"\n";
	os<<"PsimagLite version "<<PSIMAGLITE_VERSION<<"\n";
	return os;
}

