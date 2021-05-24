# Surge User Manual
## Table of Contents
* [Surge User Manual](#surge-user-manual)
   * [Table of Contents](#table-of-contents)
   * [Introduction](#introduction)
   * [Installation](#installation)
   * [Usage](#usage)
      * [Surge I/O](#surge-io)
      * [Badlists](#badlists)

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
gcc  -o surge -g -O3 -DWORDSIZE=32 -DMAXN=WORDSIZE -DOUTPROC=surgeproc \
  -march=native -DPRUNE=surgeprune -DGENG_MAIN=geng_main \
  -DZLIB -DPREPRUNE=surgepreprune surge.c geng.c planarity.c nautyW1.a 
```
You can build-in gzip output using the zlib library (https://zlib.net). Add -DZLIB to the compilation, and link with the zlib library eitherby adding -lz or libz.a . This will activate the -z command to gzip the output. The command to compile would then be:
```
gcc  -o surge -g -O3 -DWORDSIZE=32 -DMAXN=WORDSIZE -DOUTPROC=surgeproc \
  -march=native -DPRUNE=surgeprune -DGENG_MAIN=geng_main \
  -DZLIB -DPREPRUNE=surgepreprune surge.c geng.c planarity.c nautyW1.a -lz
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
  -f# -f#:#  Limit number of cycles of length 4
  -p# -p#:#  Limit number of cycles of length 5
  -b    Only rings of even length (same as only cycles of even length)
  -T    Disallow triple bonds
  -P    Require planarity
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

### Description of available options

Options are indicated by a hyphen in the Unix style, such as `-u`.

They can be concatenated or written separately: `-u -t0` is the same as `-ut0`. The exception is the option `-oFILE` (specify
output file) where the file name is everything from the `-o` to the end of the option.

If options have values, they must be written against the option name and not separately: `-t2` is valid but `-t 2` is not.

#### Output options.

`-u`  Don't write any molecules; just generate them and report the number.

`-oFILE`  Specify a file name for output, to be used instead of the standard output. You can specify names with spaces or special characters by using quotes, for example `-o'Very Many Molecules.smi'`

`-S`  Output in basic SMILES format

`-A`  Output in "alphabetic format", which looks like this 

>4 3 C2O2H4 2-3 1-2 0=1

The first two numbers indicate the numbers of atoms and bond to follow, i.e. 4 atoms and 3 bonds in this example. The atoms are numbered in the order they appear in the name     

>0=C, 1=C, 2=O, 3=O
    
Then the bonds are listed with "-" for single, "=" for double and "#" for triple.

`-O#`  This is mostly for debugging purposes. # is a number 1, 2 or 3.  For -O1, simple graphs from the initial generation are written in graph6 format. For -O2, labelled simple graphs are written in a format matching the nauty utility vcolg. For example:
        
> 4 3 2 0 2 0  0 2 0 3 1 3
    
This means 4 atoms, 3 bonds, atom types 0,2,0 and bonds 0-2, 0-3, 1-3. The atom type are C=0, N=1, O=2, P=3, S=4, F=5, Cl=6, Br=7, I=8, B=9. Bond multiplicities have not been assigned at this stage. For `-O3`, fully constructed molecules are written in a format natching the nauty utility multig. It is the same as for `-O2` except that bond multiplicities are given (0=single, 1=double, 2=triple). For example:

> 4 3 2 0 2 0  0 2 0 0 3 0 1 3 1
   
Bond (0,2) is single, (0,3) is single, (1,3) is double.

If none of `-u`, `-S`, `-A` or `-O` are given, output is in SDF format.

`-z` Except for the cases of -u and -O1, the output is gzipped. This option is only available if surge has been built with the zlib library (see [Installation](#installation)).

`-m#/#` (where each # is a number) This option tells surge to do only a part of the generation. If the two numbers are res and mod, respectively, it must be that 0 <= res <= mod-1. For example, to split the generation into five parts which can be run independently, use `-m0/5`, `-m1/5`, `-m2/5`, `-m3/5` and `-m4/5` respectively on the parts.  The overhead in this splitting is very low. However, if the number of parts is large some of the parts may be a lot bigger than others.

`-v`  Write some additional statistics to standard error (STDERR).

#### Rings and cycles.

Our description of surge options takes care to distinguish two types of closed path, which we call "rings" and "cycles". A cycle is a closed path that does not repeat any atom except that it finishes on the same atom it started on.

A ring is a cycle with the additional restriction that there is no bond between any of its atoms apart from the bonds forming the cycle. (In graph theory language, a ring is a chordless cycle, also known as an induced cycle.)

Every ring is a cycle, but a cycle may or may not be a ring.

As an example, napthalene has two rings of length 6, but three cycles (the two rings, plus a cycle of length 10). The cycle of length 10 is not a ring because of the extra bond between two of its atoms. Cheminformaticians would hold a different view on this and distinguish between the Smallest Set of Smallest Rings and the Set of All Rings.  
![Image of Napthalene](naphtalene.png)

Four surge options restrict the cycles of the molecule.

`-t#`  or  `-t#:#`  (where each # is a number) Specify the allowed number of cycles of length 3. For example, `-t0`  or  `-t1:3`.

`-f#`  or  `-f#:#`  (where each # is a number) Specify the allowed number of cycles of length 4. For example, `-f0`  or  `-f1:3`.

`-p#`  or  `-p#:#`  (where each # is a number) Specify the allowed number of cycles of length 5. For example, `-p0`  or  `-p1:3`.

`-b`  Specify that no cycles of odd length are permitted. A simple exercise is that there is a cycle of odd length if and only if there is a ring of odd length, so an equivalent description is that no rings of odd length are permitted. (In graph theory language, the graph is bipartite.)

#### Other global structural restrictions

`-e#`  or  `-e#:#`  (where each # is a number) Specify the allowed number of bonds. Bonds are counted without regard to their multiplicity. Examples: `-e15`  or `-e12:13`. If the number of atoms is n and the number of bonds is e, then e - n + 1 is the minimum number of bonds that need to be broken to reach an acyclic structure.

`-P`  Require that the molecule be planar. That is, it can be drawn in the plane with no bonds crossing.

#### Forbidden substructures

`-T`  Forbid triple bonds

`-Blist` Forbid sets of substructures.  The argument of -B is a list of numbers separated by commas without spaces. For example, in `-B2,3,8`each number indicates a set of substructures that are forbidden (for the meaning of those numbers, read on). You can use `-B` more than once, for example `-B2,4,6` is the same as `-B4 -B6,2`. We will describe the meaning of each number in the following.  

`-B1`  Rings of length up to 7 have no triple bonds. This is equivalent to cycles of length up to 7 having no triple bonds.

`-B2`  Consider rings of length r and s which share one bond (i.e. fused rings). Let e be the common bond and let f be any bond belonging to one of the rings and sharing one atom with e. In the cases {r,s} = {3,3}, {3,4} and {3,5}, both e and f must be single bonds. In the cases {r,s} = {3,6}, {4,4} and {4,5}, f must be a single bond.
![Forbidden substructure examples for option B2](option-B2.png)

-B3  Consider rings of length r and s which share two bonds.
    Let e be one of the shared bonds and let f be a bond
    belonging to one of the rings and sharing one atom
    with e.
    In the cases {r,s} = {4,4}, {4,5}, {4,6}, {5,5} and {5,6},
    both e and f must be single bonds.
      <needs a picture: Molgen cases 25-34.>

-B4  Consider two rings of length 6 that share three bonds.
    Then any bond which lies in one of the rings and has
    an atom in the other ring must be a single bond.
      <needs a picture: Molgen case 35.>

-B5  No atom has two double bonds.
      <picture of  A=A=A>

-B6  No atom in a ring of length up to 8 has two double bonds.
      <needs a picture: Molgen cases 6-11.>

-B7  These are forbidden: two atoms with four common neighbours,
     and three atoms with three common neighbours.
      <needs a picture: Molgen cases 38-39.>

-B8  These are forbidden: a cycle of length 5 having an atom
     bonded to each of the other 4 atoms, a set of 4 atoms all
     bonded to each other sharing one bond with a cycle of length 4.
      <needs a picture: Molgen cases 36-37.>

-B9  Every atom lies on at most one ring of length 3 or 4.
     Equivalently, every atom lies on at most one cycle of length
     3 or 4.
      <picture of forbidden possibilities>

================================================================


### Surge I/O
Surge outputs its results either as SMILES, SD-Files (SDF) or in a Surge-specific concise line notation. Surge output can be piped directly into toolkits such as Openbabel to analyse it or generate images of chemical diagrams.
The command

```
surge -S C6H6 | obabel -i smi - -o png -xp 1000 > test.png
```

for example, generates a 1000x1000 pixel image of the 217 isomers of benzene.

### Badlists
