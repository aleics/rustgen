extern crate rand;

use rand::Rng;

/// Onemax is a basic problem of evolutionary algorithms. The aim of the problem is to find a individuals
/// where its genes are all '1'. This is the first algorithm that I learned, when I was introduced to the
/// world of the Genetics Algorithms.
/// This algorithm is divided in different parts:
/// * Initialization of variables: the needed variables and the initial population will be initializated.
/// * Fitness calculation: the fitness will be calculated for the given problem (in this case the number of '1').
/// * Tournament selection: a combination of individuals will "fight" each other. The ones with the best fitness
///   will be selected for the next step. In total two arrays of winners will be generated
/// * Crossover: using the winners (also called parents) the new population will be created. The individuals of the
///   new population will be created as a product of the combination of the parents.
/// * Mutation: if decided, a gene of the new population's individuals will be modified.
/// * Elitism: to preserve that the best individual of the last population won't be deleted, will be saved on the
///   new population.

fn main() {
    // initialization of variables
    let pop_size: u32 = 20;
    let n_genes: u32 = 10;
    let max_gens: u32 = 500;

    // declaration of the arrays for the fitness: the best and median fitness will be saved for every generation,
    // to be able to see the progress
    let mut max_fitness_arr: Vec<u32> = Vec::with_capacity(max_gens as usize);
    let mut median_fitness_arr: Vec<f32> = Vec::with_capacity(max_gens as usize);

    // population initialization
    let mut pop: Vec<Vec<u32>> = create_random_matrix(pop_size, n_genes, 0, 1);

    // iterate for the generations
    for gen in 0..max_gens {

        // get the fitness: number of '1'
        let mut fitness: Vec<u32> = Vec::with_capacity(pop_size as usize);
        for gene in pop.iter() {
            let mut sum: u32 = 0;
            for v in gene.iter() {
                sum = sum + v;
            }

            fitness.push(sum);
        }

        // get the best indiv and its fitness
        let (max_fitness, best_indiv_idx) = get_max(&fitness);
        if best_indiv_idx == -1 {
            break;
        }
        max_fitness_arr.push(max_fitness);

        // save the best individual
        let mut best_indiv: Vec<u32> = Vec::new();
        best_indiv.clone_from(&pop[best_indiv_idx as usize]);

        // calculate the median fitness of the current generation
        let median_fitness = get_median(&fitness);
        median_fitness_arr.push(median_fitness);

        // Tournament selection
        // get the combination of individuals
        let matchup_a = create_random_matrix(pop_size, 2, 0, pop_size - 1);
        let matchup_b = create_random_matrix(pop_size, 2, 0, pop_size - 1);

        // get the winners
        let parent_a = get_parent(matchup_a, &fitness);
        let parent_b = get_parent(matchup_b, &fitness);

        // Crossover
        // determine randomly if the crossover will be made
        let do_xover = create_random_vec(pop_size, 0, 1);

        let mut new_pop: Vec<Vec<u32>> = Vec::with_capacity(pop_size as usize);
        for i in 0..pop.len() {
            let mut idx: usize;
            let mut new_gene: Vec<u32> = Vec::with_capacity(n_genes as usize);
            if do_xover[i] == 1 { // if crossover
                let crosspoint: u32 = rand::thread_rng().gen_range(0, n_genes + 1); // create a random crosspoint
                for j in 0..n_genes { // create the new gene product of both parents
                    if j <= crosspoint {
                        idx = parent_a[i] as usize;
                        new_gene.push(pop[idx][j as usize]);
                    } else {
                        idx = parent_b[i] as usize;
                        new_gene.push(pop[idx][j as usize]);
                    }
                }
            } else { // if not crossover was decided, get the parentA
                idx = parent_a[i] as usize;
                new_gene.clone_from(&pop[idx]);
            }
            new_pop.push(new_gene);
        }

        // Mutation
        // create a matrix to determine if a gene will be modified
        let do_mutation: Vec<Vec<u32>> = create_random_matrix(pop_size, n_genes, 0, 1);
        for i in 0..pop_size {
            for j in 0..n_genes {
                if do_mutation[i as usize][j as usize] == 1 { // if the gene has to be modified
                    if new_pop[i as usize][j as usize] == 1 {
                        new_pop[i as usize][j as usize] = 0;
                    } else {
                        new_pop[i as usize][j as usize] = 1;
                    }
                }
            }
        }

        pop = new_pop;

        // Elitism : save the best individual for the next generation
        pop[0] = best_indiv;

        if max_fitness == n_genes { // if ideal fitness is found
            println!("value found in {} generations", gen);
            println!("best fitness = {}", max_fitness);
            println!("median fitness = {}", median_fitness);
            break;
        }
    }
}

/// Create a random matrix between an specific range
fn create_random_matrix(size_row: u32, size_col: u32, min_range: u32, max_range: u32) -> Vec<Vec<u32>> {
    let mut mat: Vec<Vec<u32>> = Vec::with_capacity(size_row as usize);

    for _ in 0..size_row {
        let mut row: Vec<u32> = Vec::with_capacity(size_col as usize);
        for _ in 0..size_col {
            let num: u32 = rand::thread_rng().gen_range(min_range, max_range + 1);
            row.push(num);
        }
        mat.push(row);
    }

    mat
}

/// Create a random vector between an specific range
fn create_random_vec(size: u32, min_range: u32, max_range: u32) -> Vec<u32> {
    let mut v: Vec<u32> = Vec::with_capacity(size as usize);
    for _ in 0..size {
    	let num: u32 = rand::thread_rng().gen_range(min_range, max_range + 1);
        v.push(num);
    }
    v
}

/// Get the winner of a matchup given a fitness
/// Note: every matchup is an array, where every element is indexing to a population individual
fn get_parent(matchup: Vec<Vec<u32>>, fitness: &Vec<u32>) -> Vec<u32>{
    let mut parent: Vec<u32> = Vec::new();
    for i in matchup.iter() {
        let idx_first: usize = i[0] as usize;
        let idx_second: usize = i[1] as usize;

        if fitness[idx_first] >= fitness[idx_second] {
            parent.push(i[0]);
        } else {
            parent.push(i[1]);
        }
    }
    parent
}

/// Get the maximum value of a vector, and it's position
fn get_max(vec: &Vec<u32>) -> (u32, i32) {
    let mut max: u32 = 0;
    let mut index_max: i32 = -1;

    for i in 0..1 {
        for j in 0..vec.len() {
            if vec[i] > vec[j] {
                if vec[i] > max {
                    max = vec[i];
                    index_max = i as i32;
                }
            } else {
                if vec[j] > max {
                    max = vec[j];
                    index_max = j as i32;
                }
            }

        }
    }
    (max, index_max)
}

/// Get the median value of a vector
fn get_median(vec: &Vec<u32>) -> (f32) {
    let mut sum: u32 = 0;
    for i in vec.iter() {
        sum = sum + i;
    }
    (sum as f32) / (vec.len() as f32)
}
