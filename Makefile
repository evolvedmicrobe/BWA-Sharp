



all: clean bwasharp
	mkdir build
	nuget restore Bio.BWA/Bio.BWA.sln
	cp bwa_src/libbwacsharp.so build/
	xbuild /p:Configuration=Release Bio.BWA/Bio.BWA.sln
	mv Bio.BWA/Bio.BWA/bin/Release/* build/

clean:
	rm -rf build
bwasharp:
	$(MAKE) -C bwa_src/

test: all
	nuget install -Prerelease NUnit.Console -Version 3.0.0-beta-3 -OutputDirectory testrunner
	cp bwa_src/libbwacsharp.so ./Bio.BWA/TestBWA/bin/Release/
	ls -lR ./Bio.BWA/TestBWA/bin/Release/
	cd ./Bio.BWA/TestBWA/bin/Release/ && mono ./../../../../testrunner/NUnit.Console.3.0.0-beta-3/tools/nunit-console.exe TestBWA.dll
	