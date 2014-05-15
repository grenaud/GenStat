CXX      = g++ #-g -pg 
LIBGAB   = /home/gabriel_renaud/lib/
LIBTABIX = /home/gabriel_renaud/Software/tabix-0.2.6/

CXXFLAGS = -lm -O3 -Wall -Wunused-variable -lz -I${LIBGAB}  -I${LIBGAB}/VCFparser -I${LIBTABIX} -I${LIBGAB}/VCFparser/gzstream/ -c
LDFLAGS  = -lz


all: vcfcompute AlleleCounter.o ComputeDivergence.o testRandomCoordGenerator RandomGenomicCoord.o GenomicRange.o GenomicWindows.o ComputeFracAnc.o ComputeFracHetero.o testComputation MultiVCFParser.o Dstats.o DstatCounter.o generateCoords MSobject.o MSParser.o testMS divergenceMS DivergenceResult.o ComputeDivergence_core.o 


%.o: %.cpp
	${CXX} ${CXXFLAGS} $^ -o $@


testMS:	testMS.o ${LIBGAB}utils.o MSParser.o MSobject.o ${LIBGAB}VCFparser/gzstream/libgzstream.a
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS) 

divergenceMS:	divergenceMS.o ${LIBGAB}utils.o MSParser.o MSobject.o ${LIBGAB}VCFparser/gzstream/libgzstream.a
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS) 


generateCoords:	generateCoords.o ${LIBGAB}utils.o  GenomicWindows.o RandomGenomicCoord.o GenomicRange.o ${LIBGAB}VCFparser/gzstream/libgzstream.a
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

testComputation:	testComputation.o ${LIBGAB}utils.o ComputeFracAnc.o ComputeFracHetero.o ${LIBGAB}VCFparser/ReadTabix.o ${LIBGAB}VCFparser/VCFreader.o ${LIBGAB}VCFparser/CoreVCF.o ${LIBGAB}VCFparser/SimpleVCF.o ${LIBGAB}VCFparser/BAMTABLEreader.o ${LIBGAB}VCFparser/BAMTableObj.o ${LIBGAB}VCFparser/gzstream/libgzstream.a ${LIBTABIX}libtabix.a GenomicRange.o ${LIBGAB}VCFparser/FilterVCF.o  ${LIBGAB}VCFparser/SetVCFFilters.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

testRandomCoordGenerator:	testRandomCoordGenerator.o ${LIBGAB}utils.o   RandomGenomicCoord.o  GenomicRange.o GenomicWindows.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS) 

vcfcompute: vcfcompute.o ${LIBGAB}utils.o  ${LIBGAB}VCFparser/ReadTabix.o  ${LIBGAB}VCFparser/VCFreader.o ${LIBGAB}VCFparser/SimpleVCF.o  ${LIBGAB}VCFparser/CoreVCF.o ${LIBGAB}VCFparser/gzstream/libgzstream.a ${LIBGAB}VCFparser/BAMTABLEreader.o  ${LIBGAB}VCFparser/BAMTableObj.o ${LIBTABIX}libtabix.a  AlleleCounter.o ComputeDivergence.o RandomGenomicCoord.o GenomicRange.o ComputeFracAnc.o ComputeFracHetero.o  ${LIBGAB}VCFparser/FilterVCF.o  GenomicWindows.o MultiVCFParser.o Dstats.o DstatResult.o DstatCounter.o Dstat_core.o ${LIBGAB}VCFparser/SetVCFFilters.o DivergenceResult.o ComputeDivergence_core.o
	${CXX}  -o $@ $^ $(LDLIBS)  $(LDFLAGS)


clean :
	rm -f vcfcompute.o vcfcompute testRandomCoordGenerator testComputation generateCoords testMS divergenceMS  *.o


