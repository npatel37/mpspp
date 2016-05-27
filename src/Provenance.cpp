#include "Provenance.h"

std::ostream& operator<<(std::ostream& os,const Provenance&)
{
	os<<"MPS++ version "<<MPSPP_VERSION<<"\n";
	os<<"PsimagLite version "<<PSIMAGLITE_VERSION<<"\n";
	return os;
}

