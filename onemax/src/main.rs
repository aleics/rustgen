//onemax implementation as an exercise to learn Rust
extern crate rand;

use rand::Rng;

fn main() {
    //initialization of variables
    let pop_size: u32 = 10;
    let n_genes: u32 = 8;
    let rate_mut = 1/n_genes; //rate_mut: f32
    let max_gens: u32 = 500;

    //declaration of the arrays for the fitness
    let mut max_fitness_arr: Vec<u32> = Vec::with_capacity(max_gens as usize);
    let mut median_fitness_arr: Vec<f32> = Vec::with_capacity(max_gens as usize);

    //population initialization
    let mut pop: Vec<Vec<u32>> = create_random_matrix(pop_size, n_genes, 0, 1);  

    //iterate for the generations
    for _ in 0..max_gens {
        let mut fitness: Vec<u32> = Vec::with_capacity(pop_size as usize);

        for gene in pop.iter() { //get the fitness: number of '1'
            let mut sum: u32 = 0;
            for v in gene.iter() {
                sum = sum + v;
            }
            //println!("{}", sum);
            fitness.push(sum);
        }

        let (max_fitness, best_indiv) = get_max(&fitness); //get the best indiv and its fitness
        let median_fitness = get_median(&fitness); //get the median fitness

        max_fitness_arr.push(max_fitness);
        median_fitness_arr.push(median_fitness);

        println!("maxfitness = {}", max_fitness);
        println!("medianfitness = {}", median_fitness);

        //Tournament selection between the individuals
        let matchup_a = create_random_matrix(pop_size, 2, 0, pop_size - 1);
        let matchup_b = create_random_matrix(pop_size, 2, 0, pop_size - 1);

        //get the winners
        let parent_a = get_parent(matchup_a, &fitness);
        let parent_b = get_parent(matchup_b, &fitness);

        //crossover
        let do_xover = create_random_vec(pop_size, 0, 1); //determine randomly if the crossover will be made

        let mut new_pop: Vec<Vec<u32>> = Vec::with_capacity(pop_size as usize);
        for i in 0..pop.len() {
            let mut idx: usize;
            let mut new_gene: Vec<u32> = Vec::with_capacity(n_genes as usize);
            if do_xover[i] == 1 { //if crossover
                let crosspoint: u32 = rand::thread_rng().gen_range(0, n_genes + 1); //create a random crosspoint
                for j in 0..n_genes { //create the new gene product of both parents
                    if j <= crosspoint {
                        idx = parent_a[i] as usize;
                        new_gene.push(pop[idx][j as usize]);
                    } else {
                        idx = parent_b[i] as usize;
                        new_gene.push(pop[idx][j as usize]);
                    }
                }
            } else { //if not crossover get the parentA
                idx = parent_a[i] as usize;
                new_gene.clone_from(&pop[idx]);
            }
            new_pop.push(new_gene);
        }

        pop = new_pop;        
    }
}

//create a random matrix between an specific range
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

    return mat;
}

//create a random vector between an specific range
fn create_random_vec(size: u32, min_range: u32, max_range: u32) -> Vec<u32> {
    let mut v: Vec<u32> = Vec::with_capacity(size as usize);
    for _ in 0..size {
    	let num: u32 = rand::thread_rng().gen_range(min_range, max_range + 1);
        v.push(num);        
    }
    return v;
}

//get the winner of a matchup given a fitness
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
    return parent;
}

//get the maximum value of a vector, and it's position
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
    return (max, index_max); 
}

//get the median value of a vector
fn get_median(vec: &Vec<u32>) -> (f32) {
    let mut sum: u32 = 0;
    for i in vec.iter() {
        sum = sum + i;
    }
    return (sum as f32) / (vec.len() as f32);
}
