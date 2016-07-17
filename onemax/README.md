# onemax
*Onemax* is a basic problem of evolutionary algorithms. The aim of the problem is
to find a individuals where its genes are all '1'. This is the first algorithm
that I learned, when I was introduced to the world of the Genetics Algorithms.
This algorithm is divided in different parts:
* **Initialization of variables**: the needed variables and the initial population will be initialized.
* **Fitness calculation**: the fitness will be calculated for the given problem (in this case the number of '1').
* **Tournament selection**: a combination of individuals will "fight" each other. The ones with the best fitness will be selected for the next step. In total two arrays of winners will be generated
* **Crossover**: using the winners (also called parents) the new population will be created. The individuals of the new population will be created as a product of the combination of the parents.
* **Mutation**: if decided, a gene of the new population's individuals will be modified.
* **Elitism**: to preserve that the best individual of the last population won't be deleted, will be saved on the new population.
