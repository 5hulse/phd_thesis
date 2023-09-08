ls text/generated/*.py | xargs -L1 python &
cp ~/Documents/DPhil/papers/cupid/main.pdf cupid-draft.pdf
pdftk cupid-draft.pdf cat 1 output cupid-draft-page1.pdf
pdftk cupid-draft.pdf cat 2-end output cupid-draft-page2-onwards.pdf
xelatex --shell-escape thesis &&
biber thesis &&
makeindex thesis.nlo -s nomencl.ist -o thesis.nls &&
xelatex --shell-escape thesis &&
xelatex --shell-escape thesis
