CXX      = g++ #-g -pg 
LIBGAB   = libgab/
LIBTABIX = htslib/
LIBVCFPARSER = VCFparser/

#/home/gabriel_renaud/Software/tabix-0.2.6/

CXXFLAGS = -lm -O3 -Wall -Wunused-variable -lz -I${LIBGAB}  -I${LIBVCFPARSER}/ -I${LIBTABIX} -I${LIBGAB}/LIBVCFPARSER/gzstream/ -c
LDFLAGS  = -lz


all: vcfcompute AlleleCounter.o ComputeDivergence.o testRandomCoordGenerator RandomGenomicCoord.o GenomicRange.o GenomicWindows.o ComputeFracAnc.o ComputeFracHetero.o testComputation MultiVCFParser.o Dstats.o DstatCounter.o generateCoords MSobject.o MSParser.o testMS divergenceMS DivergenceResult.o ComputeDivergence_core.o 


%.o: %.cpp
	${CXX} ${CXXFLAGS} $^ -o $@


testMS:	testMS.o ${LIBGAB}utils.o MSParser.o MSobject.o ${LIBVCFPARSER}/gzstream/libgzstream.a
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS) 

divergenceMS:	divergenceMS.o ${LIBGAB}utils.o MSParser.o MSobject.o ${LIBVCFPARSER}/gzstream/libgzstream.a
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS) 


generateCoords:	generateCoords.o ${LIBGAB}utils.o  GenomicWindows.o RandomGenomicCoord.o GenomicRange.o ${LIBVCFPARSER}/gzstream/libgzstream.a
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS)

testComputation:	testComputation.o ${LIBGAB}utils.o ComputeFracAnc.o ComputeFracHetero.o ${LIBVCFPARSER}/ReadTabix.o ${LIBVCFPARSER}/VCFreader.o ${LIBVCFPARSER}/CoreVCF.o ${LIBVCFPARSER}/SimpleVCF.o ${LIBVCFPARSER}/BAMTABLEreader.o ${LIBVCFPARSER}/BAMTableObj.o ${LIBVCFPARSER}/gzstream/libgzstream.a ${LIBTABIX}libtabix.a GenomicRange.o ${LIBVCFPARSER}/FilterVCF.o  ${LIBVCFPARSER}/SetVCFFilters.o
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS) 

testRandomCoordGenerator:	testRandomCoordGenerator.o ${LIBGAB}utils.o   RandomGenomicCoord.o  GenomicRange.o GenomicWindows.o
	${CXX}  -o $@ $^ $(LDLIBS) $(LDFLAGS) 

vcfcompute: vcfcompute.o ${LIBGAB}utils.o  ${LIBVCFPARSER}/ReadTabix.o  ${LIBVCFPARSER}/VCFreader.o ${LIBVCFPARSER}/SimpleVCF.o  ${LIBVCFPARSER}/CoreVCF.o ${LIBVCFPARSER}/gzstream/libgzstream.a ${LIBVCFPARSER}/BAMTABLEreader.o  ${LIBVCFPARSER}/BAMTableObj.o ${LIBTABIX}libtabix.a  AlleleCounter.o ComputeDivergence.o RandomGenomicCoord.o GenomicRange.o ComputeFracAnc.o ComputeFracHetero.o  ${LIBVCFPARSER}/FilterVCF.o  GenomicWindows.o MultiVCFParser.o Dstats.o DstatResult.o DstatCounter.o Dstat_core.o ${LIBVCFPARSER}/SetVCFFilters.o DivergenceResult.o ComputeDivergence_core.o
	${CXX}  -o $@ $^ $(LDLIBS)  $(LDFLAGS)


clean :
	rm -f vcfcompute.o vcfcompute testRandomCoordGenerator testComputation generateCoords testMS divergenceMS  *.o


