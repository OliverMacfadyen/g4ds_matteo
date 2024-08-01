#!/bin/sh
# script to generate the pdf manual using g4dsman.tex as input

latex g4dsman.tex
latex g4dsman.tex
dvipdf g4dsman.dvi g4dsman.pdf
rm -f g4dsman.dvi g4dsman.aux g4dsman.toc g4dsman.log
