#!/bin/sh

cp scripts/controlPlotter_SL.tex $1
cd $1
pdflatex controlPlotter_SL.tex
cp controlPlotter_SL.pdf ..
cd ..
pdflatex yields.tex
cd ../..
