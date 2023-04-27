# personal_pcs="precision spectre"
# work_pcs="parsley.chem.ox.ac.uk belladonna.chem.ox.ac.uk"
# found=1
# for pc in $personal_pcs
# do
#     if [ $HOSTNAME = $pc ] ; then
#         NMRESPYPATH=/home/simon/Documents/DPhil/projects/NMR-EsPy
#         found=0
#     fi
# done

# if [ $found = 1 ]; then
#     for pc in $work_pcs
#     do
#         echo $pc
#         if [ $HOSTNAME = $pc ] ; then
#             NMRESPYPATH=/u/mf/jesu2901/DPhil/projects/spectral_estimation/NMR-EsPy
#             found=0
#         fi
#     done
# fi

# if [ $found = 1 ]; then
#     echo "Unknown PC: add to the script."
#     exit
# fi

# $NMRESPYPATH/docs/builddocs.sh
# cp $NMRESPYPATH/docs/_build/latex/nmr-espy.pdf .
xelatex --shell-escape thesis && biber thesis && makeindex thesis.nlo -s nomencl.ist -o thesis.nls && xelatex --shell-escape thesis && xelatex --shell-escape thesis
