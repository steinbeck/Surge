# Surge User Manual
## Table of Contents

## Introduction
Surge generates all chemical structures which share a given molecular formula, such as C<sub>10</sub>H<sub>16</sub>O. In the case of C<sub>10</sub>H<sub>16</sub>O, there are 452458 such structures.
Surge generates them as coordinateless chemical graphs, i.e. it determines which atoms are bonded to which and with which bond order. You can use open-source chemical toolkits such as Openbabel, rdkit or the CDK to read the output of Surge and convert the
Surge is an open source tool released under the Apache 2.0 License and relies on the Nauty package.  

## Installation
Surge binaries are available from the [https://github.com/steinbeck/Surge/releases](Surge release page) for the three major operating systems MS Windows, MacOS and Linux.

Alternatively, you can compile Surge yourself. That should work on most platforms for which a C compiler is available:
1. Download the latest Nauty release from http://users.cecs.anu.edu.au/~bdm/nauty/ and build it following the instrucations on the page.
2. Download surge.c from the [source releases on this GitHub page](https://github.com/steinbeck/Surge/releases) and put it into the nauty folder
3. Compile using the instructions at the beginning of surge.c. The following works on Linux, MacOS as well as with the https://MSYS2.org
Software Distribution and Building Platform for Windows. The latter was used to build the Windows release of Surge available on the [release page](https://github.com/steinbeck/Surge/releases) .
```
gcc -o surge -O3 -mpopcnt -march=native -DWORDSIZE=32 -DMAXN=WORDSIZE -DOUTPROC=surgeproc \
         -DPRUNE=surgeprune -DGENG_MAIN=geng_main surge.c geng.c nautyW1.a
```
You can build-in gzip output using the zlib library (https://zlib.net). Add -DZLIB to the compilation, and link with the zlib library eitherby adding -lz or libz.a . This will activate the -z command to gzip the output. The command to compile would then be:
```
gcc -o surge -O3 -mpopcnt -march=native -DWORDSIZE=32 -DMAXN=WORDSIZE -DOUTPROC=surgeproc -DZLIB -DPRUNE=surgeprune -DGENG_MAIN=geng_main surge.c geng.c nautyW1.a -lz
```



## Usage

### Surge I/O
Surge outputs its results either as SMILES, SD-Files (SDF) or in a Surge-specific concise line notation. Surge output can be piped directly into toolkits such as Openbabel to analyse it or generate images of chemical diagrams.
The command

```
surge -S C6H6 | obabel -i smi - -o png -xp 1000 > test.png
```

for example, generates a 1000x1000 pixel image of the 217 isomers of benzene.

### Badlists
