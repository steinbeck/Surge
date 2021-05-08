# Surge
## About
Surge is a chemical structure generator based on the Nauty package and thereby on the principles of canonical augmentation. 
More precisely, Surge generates all non-isomorphic constitutional isomers of a given molecular formula. 

## Usage
Surge is a command line tool. Running `surge -u C10H16O` will generate the 452458 isomers of C<sub>10</sub>H<sub>16</sub>O in 0.1s on some vanilla flavor year-2021 PC. Running `surge -S C10H16O` outputs those structurs in SMILES format. Further formats supported are SD Files (SDF) and a concise Surge-specific format.  
For large sets of structures, the -z option for compressing the output in gzip format will come in handy. 

`surge -help` will show all options:

```
Usage: surge [-oFILE] [-z|-u] [-a] [-T] [-e#|-e#:#] [-m#/#] formula

Make chemical graphs from a formula.
  Known atoms are C,B,N,P,O,S,H,Cl,F,Br,I at their lowest valences.

  formula = a formula like C8H6N2

  -O#   Output stage: 1 after geng, 2 after vcolg, 3 after multig
        Default is to write chemical formulae
  -S    Output in SMILES format
  -a    Output in alphabetical format
  -u    Just count, don't write
  -e#,  -e#:#  Restrict to given range of distinct bonds
  -t#,  -t#:#  Limit number of rings of length 3
  -f#,  -f#:#  Limit number of rings of length 4
  -p#,  -p#:#  Limit number of rings of length 5
  -T    Disallow triple bonds
  -v    Write more information to stderr
  -m#/#  Do only a part. The two numbers are res/mod where 0<=res<mod.
  -oFILE Write the output to the given file rather than to stdout.
  -z     Write output in gzip format (only if compiled with zlib)
  -G...  Anything to the end of the parameter is passed to geng
```
## Installation
### Option 1: Binary releases
Download one of the releases on this GitHub page and run it. 

### Option 2: Build from source code
1. Download the latest Nauty release from http://users.cecs.anu.edu.au/~bdm/nauty/ and build it following the instrucations on the page. 
2. Download surge.c from the [source releases on this GitHub page](https://github.com/steinbeck/Surge/releases) and put it into the nauty folder 
3. Compile using the instructions at the beginning of surge.c. The following works on Linux, MacOS as well as with the https://MSYS2.org
Software Distribution and Building Platform for Windows. The latter was used to build the Windows release of Surge available on the [release page](https://github.com/steinbeck/Surge/releases) . 
```
gcc -o surge -O3 -DWORDSIZE=32 -DMAXN=WORDSIZE -DOUTPROC=surgeproc \
         -DPRUNE=surgeprune -DGENG_MAIN=geng_main surge.c geng.c nautyW1.a
```
You can build-in gzip output using the zlib library (https://zlib.net). Add -DZLIB to the compilation, and link with the zlib library eitherby adding -lz or libz.a . This will activate the -z command to gzip the output. The command to compile would then be:
```
gcc -o surge -O3 -DWORDSIZE=32 -DMAXN=WORDSIZE -DOUTPROC=surgeproc \
         -DZLIB -DPRUNE=surgeprune -DGENG_MAIN=geng_main \
         surge.c geng.c nautyW1.a -lz
```

## Misc
Surge was developed by [Brendan McKay](http://users.cecs.anu.edu.au/~bdm) with the help of [Christoph Steinbeck](https://github.com/steinbeck) and [Mehmet Aziz Yirik](https://github.com/mehmetazizyirik). 

