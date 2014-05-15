
== Build instructions ==

First build the dependency:


cd libgab/
make utils.o FastQParser.o  FastQObj.o
cd ..
cd htslib/
make
cd ..


