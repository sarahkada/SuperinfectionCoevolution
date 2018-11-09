# SuperinfectionCoevolution
## Simulation code for the manuscript « Superinfection and the coevolution
## of parasite virulence and host recovery rate ».
## Kada & Lion, 2015.
## Dryad Digital Repository. https://doi.org/10.5061/dryad.8d6fg.
##

### Fig 1 ####
- use source("supercoev-isoclines.R ») 
- re-run the simulations (about 40 minutes).
- plot in tikz format (better aspect ratio)
Note: R is far from being the best tool to calculate isoclines, I would advise to use Mathematica for this purpose (code available in Mathematica folder).

### Fig 3, 4, 5, 6####
- open R
- For drawing figure 3, 4, 5, 6  and the figure 1A in the Annexe D, use the file supercoev > source("supercoev.R") as a source file. Warning: use R instead of Rstudio to run fig 4.
- either open the data file and plot the figure, 
- or run the simulations and draw the figure as eps or tikz. The simulations run for 10 to 90 minutes.

### Fig 7 ####
- use source("supercoev-superfunction.R") (only necessary for the eps output function)
- draw the figure by opening the csv file
OR
- re-run the simulations (less than 10 minutes).
- figure 7a was drown using the function written in the source file "supercoev-superfunction.R », with parameter used found in figure 7 legend.

### Fig 8 ####
- use source("supercoev-facilitation.R")
- open cdv table and plot figure 8
OR
- re-run simulations (less than 20 minutes)

Note that all the data that are used to draw the bistability zone are made under Mathematica, by counting the number of solution in parameter space. (see mathematica folder).
