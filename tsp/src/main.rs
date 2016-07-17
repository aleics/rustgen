extern crate csv;
extern crate rustc_serialize;

#[derive(RustcDecodable)]
struct City {
    name: String,
    pop: u32,
    long: f32,
    lat: f32
}

fn main() {
    let mut arr_cities: Vec<City> = read_datafile("data/cities.csv");

    let pop_size: usize = 10;
    let n_genes: usize = arr_cities.len();
    let max_gens: usize = 500;
}

fn read_datafile(path: &'static str) -> Vec<City> {
    let mut arr_cities: Vec<City> = Vec::new();
    let mut reader_cities = csv::Reader::from_file(path).unwrap();
    for city in reader_cities.decode() {
        let city: City = city.unwrap();
        arr_cities.push(city);
    }
    arr_cities
}
