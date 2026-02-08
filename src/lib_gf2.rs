#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GF2(u8);

use std::fmt;

impl GF2 {
    /// Creates a new GF2 element from a u8 value
    pub fn new(val: u8) -> Self {
        GF2(val % 2)
    }

    /// Returns the zero element of GF2
    pub fn zero() -> Self {
        GF2(0)
    }

    /// Returns the one element of GF2
    pub fn one() -> Self {
        GF2(1)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl std::ops::Add for GF2 {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        GF2(self.0 ^ other.0)
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl std::ops::Mul for GF2 {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        GF2(self.0 & other.0)
    }
}

impl fmt::Display for GF2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[derive(Debug, Clone)]
pub struct MatrixGF2 {
    data: Vec<Vec<GF2>>,
}

impl MatrixGF2 {
    /// Creates a new empty matrix
    pub fn new() -> Self {
        MatrixGF2 { data: Vec::new() }
    }

    /// Converts a Vec<Vec<u64>> to a MatrixGF2
    pub fn from_u64(data: Vec<Vec<u64>>) -> Self {
        let mut matrix = MatrixGF2 { data: Vec::new() };
        for row in data {
            let new_row = row.into_iter().map(|val| GF2::new(val as u8)).collect();
            matrix.data.push(new_row);
        }
        matrix
    }

    /// Returns the number of rows
    pub fn num_rows(&self) -> usize {
        self.data.len()
    }

    /// Returns the number of columns
    pub fn num_cols(&self) -> usize {
        self.data.first().map_or(0, |row| row.len())
    }

    /// Gets the element at position (i, j)
    pub fn get(&self, i: usize, j: usize) -> Option<GF2> {
        self.data.get(i).and_then(|row| row.get(j)).copied()
    }

    /// Computes the transpose of the matrix
    pub fn transpose(&self) -> Self {
        if self.data.is_empty() {
            return MatrixGF2::new();
        }

        let rows = self.num_rows();
        let cols = self.num_cols();
        let mut transposed = MatrixGF2 {
            data: vec![vec![GF2::zero(); rows]; cols],
        };

        for i in 0..rows {
            for j in 0..cols {
                transposed.data[j][i] = self.data[i][j];
            }
        }
        transposed
    }

    /// Swaps two rows of the matrix
    pub fn swap_rows(&mut self, i: usize, j: usize) {
        self.data.swap(i, j);
    }

    /// Adds one row to another (mod 2)
    pub fn add_row_to(&mut self, dest: usize, src: usize) {
        if dest < self.num_rows() && src < self.num_rows() {
            for k in 0..self.data[dest].len() {
                self.data[dest][k] = self.data[dest][k] + self.data[src][k];
            }
        }
    }
}

impl Default for MatrixGF2 {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Display for MatrixGF2 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for row in &self.data {
            for val in row {
                write!(f, "{} ", val)?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

/// Computes the Reduced Row Echelon Form (RREF) of a matrix
/// Returns the RREF matrix and a vector of pivot column indices
pub fn reduced_row_echelon_form(mut matrix: MatrixGF2) -> (MatrixGF2, Vec<usize>) {
    let num_rows = matrix.num_rows();
    let num_cols = matrix.num_cols();

    let mut pivot_out = 0;
    let mut pivots = Vec::new();

    for i in 0..num_cols {
        let pivot_in = (pivot_out..num_rows).find(|&j| matrix.get(j, i) != Some(GF2::zero()));

        if let Some(pivot_row) = pivot_in {
            matrix.swap_rows(pivot_out, pivot_row);

            for k in 0..num_rows {
                if k != pivot_out && matrix.get(k, i) == Some(GF2::one()) {
                    matrix.add_row_to(k, pivot_out);
                }
            }

            pivots.push(i);
            pivot_out += 1;
        }
    }

    (matrix, pivots)
}

/// Finds a basis for the kernel of a matrix over GF2 using Gauss-Jordan elimination
pub fn find_kernel_mod2_gauss_jordan(matrix: Vec<Vec<u64>>) -> Vec<Vec<u64>> {
    println!("Finding kernel basis using Gauss-Jordan elimination");

    let matrix_gf2 = MatrixGF2::from_u64(matrix).transpose();
    let total_elements = matrix_gf2.num_rows() * matrix_gf2.num_cols();

    if total_elements < 1000 {
        println!("Matrix:\n{}", matrix_gf2);
    }

    let (rref_matrix, pivots) = reduced_row_echelon_form(matrix_gf2);

    if total_elements < 1000 {
        println!("RREF:\n{}", rref_matrix);
    }
    println!("Pivots: {:?}", pivots);

    let num_vars = rref_matrix.num_cols();
    let free_vars: Vec<usize> = (0..num_vars).filter(|i| !pivots.contains(i)).collect();

    if total_elements < 1000 {
        println!("Free variables: {:?}", free_vars);
    }

    free_vars
        .into_iter()
        .map(|free_var| {
            let mut kernel_vector = vec![0; num_vars];
            kernel_vector[free_var] = 1;

            for (row, &pivot_col) in pivots.iter().enumerate() {
                if rref_matrix.get(row, free_var) == Some(GF2::one()) {
                    kernel_vector[pivot_col] = 1;
                }
            }
            kernel_vector
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gf2_operations() {
        assert_eq!(GF2::new(0) + GF2::new(0), GF2::zero());
        assert_eq!(GF2::new(0) + GF2::new(1), GF2::one());
        assert_eq!(GF2::new(1) + GF2::new(0), GF2::one());
        assert_eq!(GF2::new(1) + GF2::new(1), GF2::zero());

        assert_eq!(GF2::new(0) * GF2::new(0), GF2::zero());
        assert_eq!(GF2::new(0) * GF2::new(1), GF2::zero());
        assert_eq!(GF2::new(1) * GF2::new(0), GF2::zero());
        assert_eq!(GF2::new(1) * GF2::new(1), GF2::one());

        assert_eq!(GF2::new(2) + GF2::new(3), GF2::one());
        assert_eq!(GF2::new(2) * GF2::new(3), GF2::zero());
    }

    #[test]
    fn test_matrix_operations() {
        let data = vec![vec![1, 0, 1], vec![0, 1, 1]];
        let matrix = MatrixGF2::from_u64(data);

        assert_eq!(matrix.num_rows(), 2);
        assert_eq!(matrix.num_cols(), 3);
        assert_eq!(matrix.get(0, 0), Some(GF2::one()));
        assert_eq!(matrix.get(1, 2), Some(GF2::one()));

        let transposed = matrix.transpose();
        assert_eq!(transposed.num_rows(), 3);
        assert_eq!(transposed.num_cols(), 2);
        assert_eq!(transposed.get(0, 1), Some(GF2::zero()));
        assert_eq!(transposed.get(1, 0), Some(GF2::zero()));
    }

    #[test]
    fn test_rref() {
        // Test with a simple matrix
        let data = vec![vec![5, 6, 7], vec![10, 11, 13]];
        let matrix = MatrixGF2::from_u64(data);
        let (rref, pivots) = reduced_row_echelon_form(matrix);

        assert_eq!(pivots, vec![0, 1]);
        assert_eq!(rref.get(0, 0), Some(GF2::one()));
        assert_eq!(rref.get(1, 1), Some(GF2::one()));
    }

    #[test]
    fn test_find_kernel() {
        // Test with a matrix that has a non-trivial kernel
        let data = vec![vec![1, 0, 1, 1], vec![0, 1, 1, 1], vec![1, 1, 0, 0]];

        let kernel = find_kernel_mod2_gauss_jordan(data);
        assert_eq!(kernel.len(), 1); // Should have 1 free variable
        assert_eq!(kernel[0], vec![1, 1, 1]); // Expected kernel vector
    }

    #[test]
    fn test_empty_matrix() {
        let empty_matrix = MatrixGF2::new();
        assert_eq!(empty_matrix.num_rows(), 0);
        assert_eq!(empty_matrix.num_cols(), 0);

        let transposed = empty_matrix.transpose();
        assert_eq!(transposed.num_rows(), 0);
        assert_eq!(transposed.num_cols(), 0);

        let (rref, pivots) = reduced_row_echelon_form(empty_matrix);
        assert_eq!(rref.num_rows(), 0);
        assert_eq!(pivots.len(), 0);
    }
}
