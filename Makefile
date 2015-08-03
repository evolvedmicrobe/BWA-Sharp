



all: clean bwasharp
	mkdir build
	nuget restore Bio.BWA/Bio.BWA.sln
	mv bwa_src/libbwacsharp.so build/
	xbuild /p:Configuration=Release Bio.BWA/Bio.BWA.sln
	mv Bio.BWA/Bio.BWA/bin/Release/* build/

clean:
	rm -rf build
bwasharp:
	$(MAKE) -C bwa_src/