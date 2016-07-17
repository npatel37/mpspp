# Preliminaries
## Disclaimer and Licensing

MPS++ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
MPS++ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with MPS++. If not, see <http://www.gnu.org/licenses/>.
The full software license for MPS++ version 1.0.0
can be found in
file LICENSE.

## Please cite this work

MPS++ is a free and open source Matrix Product States (MPS) code.
The full software license for MPS++ version 0.1
can be found in
file LICENSE.
You are welcomed to use it and publish data
obtained with MPS++. If you do, please cite this
work. Explain How To Cite This Work. FIXME. TBW.

## Formulas used

For general discussion, see U. Schollwock, 
Ann. Phys. 326, 96 (2011)

## Code Integrity
Hash of the latest commit is also posted at
https://web.ornl.gov/~gz1/hashes.html

Latest commit should always be signed.
Keys at https://web.ornl.gov/~gz1/keys.html

### Required Software

- Item GNU C++

- Item (required) The LAPACK library.

 The configure.pl script will ask for the LDFLAGS variable
 to pass to the compiler/linker. If the linux platform was
 chosen the default/suggested LDFLAGS will include -llapack.
 If the osx platform was chosen the default/suggested LDFLAGS will
 include  -framework Accelerate.
 For other platforms the appropriate linker flags must be given.
 More information on LAPACK is at http://netlib.org/lapack/

- Item (required) PsimagLite.

 This is here https://github.com/g1257/PsimagLite/.
 You can do <code>git clone https://github.com/g1257/PsimagLite.git</code>
 in a separate directory
 outside of the MPS++ distribution. `configure.pl` will ask you where you put it.

- Item (optional) make or gmake
(only needed to use the Makefile)

- Item (optional) perl
(only needed to run the configure.pl script)

### Building MPS++
<pre>
 (go to someDirectory)
 git clone https://github.com/g1257/PsimagLite.git
 git clone https://github.com/g1257/mpspp.git
 cd PsimagLite/lib
 perl configure.pl
 (you may now edit Config.make)
 make
 cd ../../
 cd mpspp/src
 perl configure.pl
 (you may now edit Config.make)
 make
</pre>

### Running MPS++
 <code>./mpspp -f input.inp</code>

 Sample input files can be found under <code>TestSuite/inputs/</code>.

<code>configure.pl</code> creates the <code>Makefile</code>
 according to preferences in the file <code>Config.make</code>. 
If <code>Config.make</code> does not exist, <code>configure.pl</code>
 copies <code>Config.make.sample</code>
into <code>Config.make</code>, but if <code>Config.make</code>
 exists, <code>configure.pl</code> will not
overwrite it.

<!---
The MPS++ computer program simulates models of strongly correlated quantum lattice systems. The algorithm used is the density-matrix renormalization group group formulated in the matrix product state language. These models and methods are used to obtain ground state, dynamic, finite temperature, and time evolution properties of the underlying materials being simulated. These materials are usually one- and quasi-one dimensional strongly correlated electronic systems, such as Mott insulators, high-temperature superconductors, manganites and other transition metal oxides. 
-->


