# ls text/generated/*.py | xargs -L1 python &
xelatex --shell-escape thesis_bound &&
biber thesis_bound &&
makeindex thesis_bound.nlo -s nomencl.ist -o thesis_bound.nls &&
xelatex --shell-escape thesis_bound &&
xelatex --shell-escape thesis_bound
