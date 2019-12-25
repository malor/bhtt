use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use ordered_float::NotNan;
use statrs::statistics::OrderStatistics;

const QUANTILES: [f64; 14] = [
    0.0, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.75, 0.9, 0.95, 0.99, 1.0,
];

pub struct Dataset {
    values: Vec<f64>,
    quantiles: BTreeMap<NotNan<f64>, f64>,
}

impl Dataset {
    pub fn from_file(filename: &str) -> std::io::Result<Dataset> {
        let file = File::open(filename)?;

        let values: Vec<f64> = BufReader::new(file)
            .lines()
            .map(|line| line.unwrap())
            .filter(|line| !line.starts_with("#"))
            .map(|line| line.parse::<f64>().unwrap())
            .collect();

        // quantile() sorts the vector internally, so we create a copy to
        // preserve the original order of values
        let mut values_copy = values.clone();
        let quantiles = QUANTILES
            .iter()
            .map(|q| (NotNan::new(*q).unwrap(), *(&mut values_copy.quantile(*q))))
            .collect();

        Ok(Dataset {
            values: values,
            quantiles: quantiles,
        })
    }

    pub fn values(&self) -> &[f64] {
        &self.values
    }

    pub fn quantiles(&self) -> &BTreeMap<NotNan<f64>, f64> {
        &self.quantiles
    }
}
