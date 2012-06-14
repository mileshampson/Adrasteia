#!/bin/bash
echo Compiling.........................................
cd ../../Bin
swig -outdir ./ -python ../Src/Align_score.i 
swig -outdir ./ -python ../Src/agrepy.i
swig -outdir ./ -python ../Src/Bitwise.i
swig -outdir ./ -python ../Src/Statistics/gen_sequence.i
gcc -fPIC -c ../Src/*.c -I /usr/include/python2.4
gcc -fPIC -c ../Src/Statistics/*.c -I /usr/include/python2.4
ld -shared Align_score.o Align_score_wrap.o -o ./_Align_score.so 
ld -shared agrepy.o lagrepy.o sagrepy.o agrepy_wrap.o -o ./_agrepy.so 
ld -shared Bitwise.o Bitwise_wrap.o -o ./_Bitwise.so 
ld -shared gen_beta.o gen_dirch_mix.o gen_dirch.o gen_norm.o gen_sequence.o gen_sequence_wrap.o -o ./_GenSequence.so
cp ../Src/Adrasteia.py ./
echo Compiling Complete................................