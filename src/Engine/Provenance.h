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

DISCLAIMER

UT-Battelle, LLC AND THE GOVERNMENT MAKE NO REPRESENTATIONS AND DISCLAIM ALL
WARRANTIES, BOTH EXPRESSED AND IMPLIED.  THERE ARE NO EXPRESS OR IMPLIED
WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE
USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT, COPYRIGHT, TRADEMARK, OR
OTHER PROPRIETARY RIGHTS, OR THAT THE SOFTWARE WILL ACCOMPLISH THE INTENDED
RESULTS OR THAT THE SOFTWARE OR ITS USE WILL NOT RESULT IN INJURY OR DAMAGE.
THE USER ASSUMES RESPONSIBILITY FOR ALL LIABILITIES, PENALTIES, FINES, CLAIMS,
CAUSES OF ACTION, AND COSTS AND EXPENSES, CAUSED BY, RESULTING FROM OR ARISING
OUT OF, IN WHOLE OR IN PART THE USE, STORAGE OR DISPOSAL OF THE SOFTWARE.

*********************************************************

*/
/** \ingroup MPSPP */
/*@{*/

#ifndef PROVENANCE_H
#define PROVENANCE_H
#include "Version.h" // do not commit Version.h, it is created dynamically

class Provenance {

public:

}; // Provenance

std::ostream& operator<<(std::ostream& os,const Provenance &prov)
{
	os<<"MPS++: revision: "<<mpsppRevision<<"\n";
	os<<"MPS++: diff: "<<mpsppDiff<<"\n";
	os<<"PsimagLite: revision: "<<psimagLiteRevision<<"\n";
	os<<"PsimagLite: diff: "<<psimagLiteDiff<<"\n";
	return os;
}

/*@}*/
#endif // PROVENANCE_H

