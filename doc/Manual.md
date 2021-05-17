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

Surge is a command line tool. Running `surge -u C10H16O` will generate the 452458 isomers of C<sub>10</sub>H<sub>16</sub>O in 0.1s on some vanilla flavor year-2021 PC. Running `surge -S C10H16O` outputs those structurs in SMILES format. Further formats supported are SD Files (SDF) and a concise Surge-specific format.  
For large sets of structures, the -z option for compressing the output in gzip format will come in handy.

`surge -help` will show all options:

```
Usage: surge [-oFILE] [-z] [-u|-A|-S] [-T] [-e#|-e#:#] [-m#/#] formula

Make chemical graphs from a formula.
  Known atoms are C,B,N,P,O,S,H,Cl,F,Br,I at their lowest valences.

  formula = a formula like C8H6N2

  -O#   Output stage: 1 after geng, 2 after vcolg, 3 after multig
        Default is to write SD format
  -S    Output in SMILES format
  -A    Output in alphabetical format
  -u    Just count, don't write
  -e# -e#:#  Restrict to given range of distinct bonds
  -t# -t#:#  Limit number of rings of length 3
  -f# -f#:#  Limit number of rings of length 4
  -p# -p#:#  Limit number of rings of length 5
  -b    Only rings of even length
  -T    Disallow triple bonds
  -B#,...,# Specify sets of substructures to avoid
     1 = no triple bonds in rings up to length 7
     2 = Bredt's rule for two rings ij with one bond in
           common (33, 34, 35, 36, 44, 45)
     3 = Bredt's rule for two rings ij with two bonds in
           common (i,j up to 56)
     4 = Bredt's rule for two rings of length 6 sharing three bonds
     5 = no substructures A=A=A (in ring or not)
     6 = no substructures A=A=A in rings up to length 8
     7 = no K_33 or K24 structure
     8 = none of cone of P4 or K4 with 3-ear
     9 = no atom on more than one ring of length 3 or 4
  -v    Write more information to stderr
  -m#/#  Do only a part. The two numbers are res/mod where 0<=res<mod.
  -oFILE Write the output to the given file rather than to stdout.
  -z     Write output in gzip format (only if compiled with zlib)
```

### Surge I/O
Surge outputs its results either as SMILES, SD-Files (SDF) or in a Surge-specific concise line notation. Surge output can be piped directly into toolkits such as Openbabel to analyse it or generate images of chemical diagrams.
The command

```
surge -S C6H6 | obabel -i smi - -o png -xp 1000 > test.png
```

for example, generates a 1000x1000 pixel image of the 217 isomers of benzene.

### Badlists
