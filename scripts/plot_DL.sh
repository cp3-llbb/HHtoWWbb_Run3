#!/bin/sh

cp scripts/controlPlotter_DL.tex $1
cd $1
pdflatex controlPlotter_DL.tex
cp controlPlotter_DL.pdf ..
cd ..
pdflatex yields.tex
cd ..
