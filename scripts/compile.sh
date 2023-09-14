ls text/generated/*.py | xargs -L1 python &
cp ~/Documents/DPhil/papers/cupid/main_angew.pdf cupid-draft.pdf
xelatex --shell-escape thesis &&
biber thesis &&
makeindex thesis.nlo -s nomencl.ist -o thesis.nls &&
xelatex --shell-escape thesis &&
xelatex --shell-escape thesis
