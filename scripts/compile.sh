# ls text/generated/*.py | xargs -L1 python &
xelatex --shell-escape thesis &&
biber thesis &&
makeindex thesis.nlo -s nomencl.ist -o thesis.nls &&
xelatex --shell-escape thesis &&
xelatex --shell-escape thesis
