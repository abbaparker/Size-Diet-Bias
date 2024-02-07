# Size-Diet-Bias

This repository contains the code used for the analyses and figures in the manuscript "An assessment of body size and dietary biases in fossil mammal assemblages of the Pleistocene of Eurasia." 

The R code is divided in 4 scripts, which work best when run sequentially.
1- Takes the input of a download of occurrences from the NOW database (or similar) and creates a site/occurrence matrix to use for faunal analyses. Appendix 1 figure.
2- Compares generic occurrence counts at fossil sites and the nearest modern Eco-ISEA3H hexagon community sample. Figs 1, 4, 5. Appendix 5 figure. 
3- Finds the body mass distributions at fossil sites, runs all models of bimodal body size distributions, . Figs. 2, 3. Appendix 2, 3 figures.
4- Counts the number of large/small, herbivore/non-herbivore genera in fossil and modern faunas. Runs GLLVM on these guild counts, comparing site types to modern. Figs. 6, 7. 

The Data folder contains all the data necessary for the fossil faunal analyses. Due to file size, the Eco-ISEA3H Phylacine samples (GenusMatrix for modern communities) are not available here; if you would like to see this data, email abigail.parker@helsinki.fi . Further information on Phylacine https://megapast2future.github.io/PHYLACINE_1.2/ , and Eco-ISEA3H https://etsin.fairdata.fi/dataset/552c6ac2-4677-4a7b-952f-1632b2b9c335
