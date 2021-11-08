VERSION=0_9

# In the following, define NAUTY to be the directory where the
# nauty include files are, and NAUTYLIB to be where libraries
# like nautyW1.a are. 

NAUTY=${NAUTY_HOME}
NAUTYLIB=${NAUTY}

SURGE=-g -I ${NAUTY} -DWORDSIZE=32 -DMAXN=WORDSIZE -DOUTPROC=surgeproc \
  -march=native -mtune=native -DPREPRUNE=surgepreprune \
  -DPRUNE=surgeprune -DGENG_MAIN=geng_main

all : surge surge_0 surge_1 surge_2 surge_3

surge : surge.c 
	gcc -o surge ${SURGE} -O3 -DZLIB \
	   surge.c geng.c planarity.c ${NAUTYLIB}/nautyW1.a -lz

surge_0 : surge.c plugin0.c
	gcc -o surge_0 ${SURGE} -O3 -DZLIB -DSURGEPLUGIN='"plugin0.c"' \
	   surge.c geng.c planarity.c ${NAUTYLIB}/nautyW1.a -lz

surge_1 : surge.c plugin1.c
	gcc -o surge_1 ${SURGE} -O3 -DZLIB -DSURGEPLUGIN='"plugin1.c"' \
	   surge.c geng.c planarity.c ${NAUTYLIB}/nautyW1.a -lz

surge_2 : surge.c plugin2.c
	gcc -o surge_2 ${SURGE} -O3 -DZLIB -DSURGEPLUGIN='"plugin2.c"' \
	   surge.c geng.c planarity.c ${NAUTYLIB}/nautyW1.a -lz

surge_3 : surge.c plugin3.c
	gcc -o surge_3 ${SURGE} -O3 -DZLIB -DSURGEPLUGIN='"plugin3.c"' \
	   surge.c geng.c planarity.c ${NAUTYLIB}/nautyW1.a -lz

# This one is unoptimized for easier debugging
surge_g : 
	gcc -o surge_g ${SURGE} -DZLIB \
	   surge.c geng.c planarity.c ${NAUTYLIB}/nautyW1.a -lz

canonsdf : canonsdf.c
	gcc -o canonsdf -I{NAUTY} canonsdf.c ${NAUTYLIB}/nautyL1.a

tarfile :
	rm -rf surge${VERSION}
	mkdir surge${VERSION}
	cp SDFformats.txt surge.c planarity.c geng.c makefile \
	  B[1-9].sdf canonsdf.c plugin[0-3].c surge.txt surge${VERSION}/
	touch surge${VERSION}/This_is_surge_${VERSION}.txt
	tar cvf surge${VERSION}.tar surge${VERSION}/
