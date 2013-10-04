# SVN: $Id: Makefile 175 2013-06-26 01:37:49Z mushkar $
all: dirs
	cd ./src; \
	make all; \
	cd ..

Fit: dirs
	cd ./src; \
	make Fit; \
	cd ..

dirs:
	mkdir -p ./bin; \
	mkdir -p ./obj; \
	mkdir -p ./lib; \

clean:
	cd ./src; \
	make clean; \
	cd ..

print:
	cd ./src; \
	make print; \
	cd ..

docs:
	cd ./src; \
	make docs; \
	cd ..
