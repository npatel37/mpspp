# MPS++ (v0.1)
# Platform: Linux
# MPI: 0

LDFLAGS =      /usr/lib64/liblapack.so.3 /usr/lib64/libblas.so.3      -lpthread 
CPPFLAGS = -Werror -Wall  -IEngine -I../../PsimagLite/src -I../../PsimagLite/src/Geometry  \
 -IModels/HubbardOneOrbital
CXX = g++  -O3 -DNDEBUG -g3
#Comment out line below for debugging (COMMENT ALSO THE strip commands): 
#CXX = g++ -g3 #ALSO COMMENT OUT strip command below
EXENAME = mpspp 
all: $(EXENAME)

mpspp:  mpspp.o gitrev
	$(CXX) -o mpspp mpspp.o $(LDFLAGS)  

# dependencies brought about by Makefile.dep
%.o: %.cpp Makefile gitrev Engine/Version.h
	$(CXX) $(CPPFLAGS) -c $< 

Makefile.dep: Engine/Version.h mpspp.cpp
	$(CXX) $(CPPFLAGS) -MM mpspp.cpp  > Makefile.dep

Engine/Version.h: gitrev
	./gitrev > Engine/Version.h	

gitrev: gitrev.o
	$(CXX) -o gitrev gitrev.o $(LDFLAGS)

gitrev.o:
	$(CXX) $(CPPFLAGS) -c gitrev.cpp 

doc:
	cd ../doc; make

clean:
	rm -f core* $(EXENAME) *.o *.dep Engine/Version.h gitrev

include Makefile.dep
######## End of Makefile ########

