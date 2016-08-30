/// The Travel Salesman Problem (TSP) is famous problem that can be solved using evolutionary
/// algorithms. The goal of this algorithm is to visit each city and comeback to the first one
/// using the most efficient way (the shortest distance travelled). The main idea to solve this
/// problem is:
/// * As an input data, a couple of cities with them location (latitude and longitude) will be
///   read. In this case cities of Germany.
/// * A couple of combinations of the cities will be generated, what it's called population.
///   For each combination, the complete distance will be calculated.
/// * Using Tournament Selection the best combinations will be identified. These combinations will
///   be modified using Crossover and Mutation to create the new population.
/// * After some generations, the best combination of cities will be calculated.


extern crate csv;
extern crate rustc_serialize;
extern crate rand;

use rand::Rng;
use std::f32;

#[derive(RustcDecodable)]
struct City {
    name: String,
    pop: u32,
    long: f32,
    lat: f32
}

fn main() {
    let arr_cities: Vec<City> = read_datafile("data/cities.csv");

    let pop_size: usize = 10;
    let n_genes: usize = arr_cities.len();
    let max_gens: usize = 1500;

    let poss = create_succ_vec(n_genes);
    let mut best_fitness_arr: Vec<f32> = Vec::with_capacity(max_gens as usize);
    let mut median_fitness_arr: Vec<f32> = Vec::with_capacity(max_gens as usize);

    let mut mat_dist: Vec<Vec<f32>> = Vec::with_capacity(n_genes);
    for i in 0..n_genes {
        let mut pr_arr_dist: Vec<f32> = Vec::with_capacity(n_genes);
        for j in 0..n_genes {
            pr_arr_dist.push(calcul_eucl_distance(arr_cities[i].long,
                                                  arr_cities[i].lat,
                                                  arr_cities[j].long,
                                                  arr_cities[j].lat));
        }

        mat_dist.push(pr_arr_dist);
    }

    let mut pop: Vec<Vec<u32>> = Vec::with_capacity(n_genes);
    for _ in 0..pop_size {
        pop.push(randperm(&poss));
    }

    for gen in 0..max_gens {

        let mut fitness: Vec<f32> = vec![0.0; pop_size];
        for i in 0..pop_size {
            for j in 0..(n_genes-1) {
                fitness[i] = fitness[i] + mat_dist[pop[i][j] as usize][pop[i][j + 1] as usize];
            }
            fitness[i] = fitness[i] + mat_dist[pop[i][n_genes-1] as usize][pop[i][1] as usize];
        }

        let (min_fitness, best_indiv_idx) = get_min_f32(&fitness);
        if best_indiv_idx == -1 {
            break;
        }

        best_fitness_arr.push(min_fitness);

        // save the best individual
        let mut best_indiv: Vec<u32> = Vec::new();
        best_indiv.clone_from(&pop[best_indiv_idx as usize]);

        // calculate the median fitness of the current generation
        let median_fitness = get_median(&fitness);
        median_fitness_arr.push(median_fitness);

        // Tournament selection
        // get the combination of individuals
        let matchup_a = create_random_matrix(pop_size, 2, 0, (pop_size - 1) as u32);
        let matchup_b = create_random_matrix(pop_size, 2, 0, (pop_size - 1) as u32);

        // get the winners
        let parent_a = get_parent_sm(matchup_a, &fitness);
        let parent_b = get_parent_sm(matchup_b, &fitness);

        // Crossover
        let rate_cross: f32 = 0.8;

        let mut new_pop: Vec<Vec<u32>> = Vec::with_capacity(pop_size);
        for i in 0..pop.len() {
            let mut idx: usize;
            let mut new_gene: Vec<u32> = Vec::with_capacity(n_genes);

            let do_xover: f32 = rand::thread_rng().gen_range(0.0, 1.0); // determine randomly if the crossover will be made

            if do_xover > rate_cross { // if crossover
                let crosspoint: u32 = rand::thread_rng().gen_range(0, (n_genes as u32) + 1); // create a random crosspoint

                for j in 0..n_genes { // create the new gene product of both parents
                    if (j as u32) <= crosspoint {
                        idx = parent_a[i] as usize;

                        if contains(&new_gene, pop[idx][j]) == false { // if the new gene is not present
                            new_gene.push(pop[idx][j]);
                        }
                    } else {
                        idx = parent_b[i] as usize;
                        if contains(&new_gene, pop[idx][j]) == false { // if the new gene is not present
                            new_gene.push(pop[idx][j]);
                        }
                    }
                }
            } else { // if not crossover was decided, get the parentA
                idx = parent_a[i] as usize;
                new_gene.clone_from(&pop[idx]);
            }

            new_gene = complete_gene(&new_gene, n_genes as u32); // complete the gene if cities are missing

            new_pop.push(new_gene);
        }

        for i in 0..new_pop.len() {
            // to use a swapping method uncomment this blog and comment the mutation_block function
            /*let point_a: usize = rand::thread_rng().gen_range(0, new_pop[i].len());
            let point_b: usize = rand::thread_rng().gen_range(0, new_pop[i].len());

            new_pop[i] = swap_vect(&new_pop[i], point_a, point_b);*/
            new_pop[i] = mutation_block(&new_pop[i], 3);
        }

        pop = new_pop; // copy the population for the next iteration

        // Elitism : save the best individual for the next generation
        pop[0] = best_indiv;

        // final plot
        if gen == max_gens-1 {
            println!("best fitness progression:");
            for i in 0..max_gens {
                println!("{} ", best_fitness_arr[i]);
            }

            println!("best indiv: ");
            for i in 0..pop[0].len() {
                print!("{} - ", arr_cities[pop[0][i] as usize].name);
            }
            print!("{}", arr_cities[pop[0][0] as usize].name);
            println!("");
        }
    }
}

// Read the CSV data file
fn read_datafile(path: &'static str) -> Vec<City> {
    let mut arr_cities: Vec<City> = Vec::new();
    let mut reader_cities = csv::Reader::from_file(path).unwrap();
    for city in reader_cities.decode() {
        let city: City = city.unwrap();
        arr_cities.push(city);
    }
    arr_cities
}

// Create a matrix of zeros
fn create_zeros_matrix(size_row: usize, size_col: usize) -> Vec<Vec<f32>> {
    let mut mat: Vec<Vec<f32>> = Vec::with_capacity(size_row as usize);

    for _ in 0..size_row {
        let mut row: Vec<f32> = Vec::with_capacity(size_col as usize);
        for _ in 0..size_col {
            row.push(0.0);
        }
        mat.push(row);
    }

    mat
}

/// Create a random matrix between an specific range
fn create_random_matrix(size_row: usize, size_col: usize, min_range: u32, max_range: u32) -> Vec<Vec<u32>> {
    let mut mat: Vec<Vec<u32>> = Vec::with_capacity(size_row as usize);

    for _ in 0..size_row {
        let mut row: Vec<u32> = Vec::with_capacity(size_col);
        for _ in 0..size_col {
            let num: u32 = rand::thread_rng().gen_range(min_range, max_range + 1);
            row.push(num);
        }
        mat.push(row);
    }

    mat
}

/// Create a random vector between an specific range
fn create_random_vec(size: usize, min_range: u32, max_range: u32) -> Vec<u32> {
    let mut v: Vec<u32> = Vec::with_capacity(size);
    for _ in 0..size {
    	let num: u32 = rand::thread_rng().gen_range(min_range, max_range + 1);
        v.push(num);
    }
    v
}

// Create a vector of successive numbers: (0, 1, 2, ...)
fn create_succ_vec(max: usize) -> Vec<u32> {
    let mut vect: Vec<u32> = Vec::with_capacity(max);
    for i in 0..max {
        vect.push(i as u32);
    }
    vect
}

// Calculate the Eucledian distance between two points
fn calcul_eucl_distance(long1: f32, lat1: f32, long2: f32, lat2:f32) -> f32 {
    ((long1 - long2).powf(2.0) + (lat1 - lat2).powf(2.0)).sqrt()
}

// Return a vector that represents a random permutation of the vector given
fn randperm(vect: &Vec<u32>) -> Vec<u32> {
    let mut res = vect.clone();
    rand::thread_rng().shuffle(&mut res);

    res
}

// Check if a vector contains a value
fn contains(vec: &Vec<u32>, val: u32) -> bool {
    for v in vec.iter() {
        if *v == val {
            return true;
        }
    }
    false
}

/// Get the minimum value of a vector, and it's position
fn get_min_f32(vec: &Vec<f32>) -> (f32, i32) {
    let mut min: f32 = f32::INFINITY;
    let mut index_min: i32 = -1;

    for i in 0..1 {
        for j in 0..vec.len() {
            if vec[i] < vec[j] {
                if vec[i] < min {
                    min = vec[i];
                    index_min = i as i32;
                }
            } else {
                if vec[j] < min {
                    min = vec[j];
                    index_min = j as i32;
                }
            }

        }
    }
    (min, index_min)
}

/// Get the median value of a vector
fn get_median(vec: &Vec<f32>) -> (f32) {
    let mut sum: f32 = 0.0;
    for i in vec.iter() {
        sum = sum + i;
    }
    (sum) / (vec.len() as f32)
}

/// Get the winner of a matchup given a fitness
/// Note: every matchup is an array, where every element is indexing to a population individual
fn get_parent_sm(matchup: Vec<Vec<u32>>, fitness: &Vec<f32>) -> Vec<u32>{
    let mut parent: Vec<u32> = Vec::new();
    for i in matchup.iter() {
        let idx_first: usize = i[0] as usize;
        let idx_second: usize = i[1] as usize;

        if fitness[idx_first] <= fitness[idx_second] {
            parent.push(i[0]);
        } else {
            parent.push(i[1]);
        }
    }
    parent
}

// Complete the genoma of an individual with the missing genes
fn complete_gene(input: &Vec<u32>, max: u32) -> Vec<u32> {
    let mut output = input.clone();
    for i in 0..max {
        if contains(&output, i) == false {
            output.push(i as u32);
        }
    }

    output
}

// Swap two members in a vector
fn swap_vect(gene: &Vec<u32>, point_a: usize, point_b: usize) -> Vec<u32> {
    let mut copy_gene = gene.clone();
    copy_gene.swap(point_a, point_b);

    copy_gene.clone()
}

// Realize a mutation block of the individual given its genes and a block size
fn mutation_block(gene: &Vec<u32>, block_size: usize) -> Vec<u32> {
    let n_genes: usize = gene.len();
    let mut out: Vec<u32> = gene.clone();

    let num: f32 = rand::thread_rng().gen_range(0.0, 1.0);
    let rate_mut: f32 = 1.0/(n_genes as f32);

    if num > rate_mut {
        let point_a: usize = rand::thread_rng().gen_range(0, n_genes - block_size);
        let point_b: usize = rand::thread_rng().gen_range(0, n_genes - block_size);

        out = swap_block_vect(&out, point_a, point_b, block_size);
    }

    out
}

// Realize a swapping of a full block in a vector given two indexes and a block size
fn swap_block_vect(gene: &Vec<u32>, point_a: usize, point_b: usize, block_size: usize) -> Vec<u32> {
    if ((point_b + block_size) >= gene.len()) || ((point_a + block_size) >= gene.len()) {
        panic!("block size and point bigger than gene size. info: point_a: {}, point_b: {}, block_size: {}", point_a, point_b, block_size);
    }

    let mut copy_gene = gene.clone();

    for i in 0..block_size {
        copy_gene.swap(point_a + i, point_b + i);
    }

    copy_gene.clone()
}
