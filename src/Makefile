CXX = g++

CXXFLAGS = -O3 -std=c++11 -mfpmath=sse -msse -msse2 -msse3 -DEVIDENCE_SR -DEVIDENCE_PARS -I./ \
-DNDEBUG -W -pipe -Wundef -Winline --param large-function-growth=100000 -Wall -fPIC
ICCFLAGS = -O3 -xCORE-AVX-I -Wall

LINKFLAGS = -lm
GDLINKFLAGS = -lgd -lpng

CONTRAFOLD_SRCS = \
	Contrafold.cpp \
	ParameterManager.cpp \
	ComputationEngine.cpp \
	ComputationWrapper.cpp \
	OptimizationWrapper.cpp \
	OuterOptimizationWrapper.cpp \
	CGOptimizationWrapper.cpp \
	CGLinear.cpp \
	InnerOptimizationWrapper.cpp \
	InnerOptimizationWrapperEM.cpp \
	EM.cpp \
	InnerOptimizationWrapperLBFGS.cpp \
	LBFGS.cpp \
	InnerOptimizationWrapperSubgradientMethod.cpp \
	SubgradientMethod.cpp \
	InnerOptimizationWrapperStochasticGradient.cpp \
	LineSearch.cpp \
	InferenceEngine.cpp \
	FileDescription.cpp \
	Options.cpp \
	SStruct.cpp \
	GammaMLE.cpp \
	Utilities.cpp
	# SparseMatrix.cpp \

MAKECOORDS_SRCS = \
	MakeCoords.cpp \
	SStruct.cpp \
	LBFGS.cpp \
	LineSearch.cpp \
	Utilities.cpp

PLOTRNA_SRCS = \
	EncapsulatedPostScript.cpp \
	PlotRNA.cpp \
	SStruct.cpp \
	Utilities.cpp

SCOREPREDICTION_SRCS = \
	ScorePrediction.cpp \
	SStruct.cpp \
	Utilities.cpp

CONTRAFOLD_OBJS = $(CONTRAFOLD_SRCS:%.cpp=%.o)
MAKECOORDS_OBJS = $(MAKECOORDS_SRCS:%.cpp=%.o)
PLOTRNA_OBJS = $(PLOTRNA_SRCS:%.cpp=%.o)
SCOREPREDICTION_OBJS = $(SCOREPREDICTION_SRCS:%.cpp=%.o)

.PHONY: all viz api clean

all: contrafold score_prediction api_test
viz: make_coords plot_rna
api: api_test

contrafold: $(CONTRAFOLD_OBJS)
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) $(CONTRAFOLD_OBJS) $(LINKFLAGS) -o contrafold

Defaults.hpp: MakeDefaults.pl *.params.*
	perl MakeDefaults.pl contrafold.params.complementary contrafold.params.noncomplementary contrafold.params.profile

Contrafold.o: Contrafold.cpp Defaults.hpp
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) -c Contrafold.cpp

api_test: api_test.o contrafold_api.o FileDescription.o Options.o SStruct.o Utilities.o ParameterManager.o InferenceEngine.o ComputationEngine.o
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) FileDescription.o Options.o SStruct.o Utilities.o ParameterManager.o InferenceEngine.o ComputationEngine.o contrafold_api.o api_test.o $(LINKFLAGS) -o api_test

contrafold_api.o: contrafold_api.cpp
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) -c contrafold_api.cpp

api_test.o: api_test.cpp
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) -c api_test.cpp

make_coords: $(MAKECOORDS_OBJS)
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) $(MAKECOORDS_OBJS) $(LINKFLAGS) -o make_coords

plot_rna: $(PLOTRNA_OBJS)
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) $(PLOTRNA_OBJS) $(LINKFLAGS) $(GDLINKFLAGS) -o plot_rna

score_prediction: $(SCOREPREDICTION_OBJS)
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) $(SCOREPREDICTION_OBJS) $(LINKFLAGS) -o score_prediction

##########################################

gccathlon64:
	make all OTHERFLAGS="-march=athlon64 -fomit-frame-pointer -ffast-math -funroll-all-loops -funsafe-math-optimizations -fpeel-loops --param max-inline-insns-single=100000 --param inline-unit-growth=100000 --param large-function-growth=100000"

gccathlon64multi:
	make all CXX="mpiCC" OTHERFLAGS="-DMULTI -march=athlon64 -fomit-frame-pointer -ffast-math -funroll-all-loops -funsafe-math-optimizations -fpeel-loops --param max-inline-insns-single=100000 --param inline-unit-growth=100000 --param large-function-growth=100000"

gccathlon64multiprofile:
	make all CXX="mpiCC" OTHERFLAGS="-DMULTI -march=athlon64 -fomit-frame-pointer -ffast-math -funroll-all-loops -funsafe-math-optimizations -fpeel-loops --param max-inline-insns-single=100000 --param inline-unit-growth=100000 --param large-function-growth=100000 -pg -g"

intel:
	make all CXX="icpc" CXXFLAGS="$(ICCFLAGS)" OTHERFLAGS="-no-ipo -static"

intelmulti:
	make all LAMHCP="icpc" CXX="mpicc" CXXFLAGS="$(ICCFLAGS)" OTHERFLAGS="-DMULTI -no-ipo"

multi:
	make all CXX="mpiCC" OTHERFLAGS="-DMULTI"

gccp4:
	make all OTHERFLAGS="-march=pentium4 -fomit-frame-pointer -ffast-math -funroll-all-loops -funsafe-math-optimizations -fpeel-loops --param max-inline-insns-single=100000 --param large-function-growth=100000 --param inline-unit-growth=100000"

gccp4multi:
	make all CXX="mpiCC" OTHERFLAGS="-DMULTI -march=pentium4 -fomit-frame-pointer -ffast-math -funroll-all-loops -funsafe-math-optimizations -fpeel-loops --param max-inline-insns-single=100000 --param large-function-growth=100000 --param inline-unit-growth=100000"

gcc64:
	make all OTHERFLAGS="-m64 -funroll-all-loops -fomit-frame-pointer -ffast-math"

profile:
	make all OTHERFLAGS="-pg -g"

gccp4profile:
	make all OTHERFLAGS="-march=pentium4 -ffast-math -funroll-all-loops -funsafe-math-optimizations -fpeel-loops --param max-inline-insns-single=100000 --param large-function-growth=100000 --param inline-unit-growth=100000 -pg -g"

debug:
	make all CXXFLAGS="-g -fno-inline -W -Wall"

debugmulti:
	make all CXX="mpiCC" OTHERFLAGS="-DMULTI" CXXFLAGS="-g -fno-inline -W -Wall"

assembly:
	make all OTHERFLAGS="-Wa,-a,-ad"

%.o: %.cpp *.hpp *.ipp
	$(CXX) $(CXXFLAGS) $(OTHERFLAGS) -c $<

clean:
	rm -f contrafold make_coords plot_rna score_prediction *.o Defaults.hpp
