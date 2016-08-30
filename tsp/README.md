# Travel Salesman Problem
The Travel Salesman Problem (TSP) is famous problem that can be solved using evolutionary
algorithms. The goal of this algorithm is to visit each city and comeback to the first one
using the most efficient way (the shortest distance travelled). The main idea to solve this
problem is:

* As an input data, a couple of cities with them location (latitude and longitude) will be
  read. In this case cities of Germany.
* A couple of combinations of the cities will be generated, what it's called population.
  For each combination, the complete distance will be calculated.
* Using Tournament Selection the best combinations will be identified. These combinations will
  be modified using Crossover and Mutation to create the new population.
* The created new population will be used on the next iteration.
* After some iterations, the best combination of cities will be calculated.

On this algorithm the fitness is the distance that a combination of cities have.
Of course, the best fitness will be the minimum fitness of a population.

## Crossover method
The crossover method used is the combination of both parents using a crossing
point.

## Mutation method
On the mutation part, it's possible currently to use two methods:

* **Swapping**: using the function *swap_vect*, two elements of the vector (gene) will be
  swapped.
* **Mutation block**: a whole block of the vector will be swapped. The block size is
  configurable. This method works theoretically better, since will not destroy a
  good succession of elements.
