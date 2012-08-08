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

/*! \file ProgramGlobals.h
 *
 *
 *
 */
#ifndef PROGRAM_LIMITS_H
#define PROGRAM_LIMITS_H

namespace Mpspp {
struct ProgramGlobals {
	//		static size_t const MaxNumberOfSites = 300; // max number of sites that a model can use
	//		static size_t const MaxLanczosSteps = 1000000; // max number of internal Lanczos steps
	//		static size_t const LanczosSteps = 200; // max number of external Lanczos steps
	//		static double const LanczosTolerance; // tolerance of the Lanczos Algorithm
	//		enum {INFINITE=0,EXPAND_ENVIRON=1,EXPAND_SYSTEM=2};
	//		enum {SYSTEM_SYSTEM,SYSTEM_ENVIRON,ENVIRON_SYSTEM,ENVIRON_ENVIRON};
	//		enum {SYSTEM,ENVIRON};
	//		enum {FERMION,BOSON};
}; // ProgramGlobals

//	double const ProgramGlobals::LanczosTolerance = 1e-12;
}; // namespace Mpspp
/*@}*/
#endif

