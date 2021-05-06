# Surge
## About
Surge is a chemical structure generator based on the Nauty package and thereby on the principles of canonical augmentation. 
More precisely, Surge generates all non-isomorphic constitutional isomers of a given molecular formula. 

## Usage
Surge is a command line tool. Running `surge -u C10H16O` will generate the 452458 isomers of C<sub>10</sub>H<sub>16</sub>O in 0.1s on some vanilla flavoe year-2021 PC. Running `surge -S C10H16O` outputs those structurs in SMILES format. Further formats supported are SD Files (SDF) and a concise Surge-specific format.  
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
2. Download surge.c from this GitHub page and put it into the nauty folder 
3. Compile using the instructions at the beginning of surge.c. 

## Misc
Surge was developed by Brendan McKay with the help of Christoph Steinbeck and Mehmet Aziz Yirik. 

