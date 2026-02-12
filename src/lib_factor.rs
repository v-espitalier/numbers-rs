use rand::{RngExt, SeedableRng};

use crate::lib_gf2::find_kernel_mod2_gauss_jordan;
use crate::lib_utils::{bit_length, digit_length, gcd, generate_primes_sieve, is_perfect_square};
use std::f64::consts::E;

/// Computes Q(x) = (x + t)^2 - n for the quadratic sieve algorithm.
#[inline]
pub fn quadratic_residue(x: u64, t: u64, n: u64) -> u64 {
    (x + t)
        .checked_mul(x + t)
        .and_then(|square| square.checked_sub(n))
        .unwrap_or(0)
}

/// Factorizes y using the factor base qs_b.
/// Returns Some(exponents) if y is B-smooth, None otherwise.
pub fn factor_with_base(y: u64, factor_base: &[u64]) -> Option<Vec<u64>> {
    // Case where y is directly in the factor base
    if let Some(pos) = factor_base.iter().position(|&b| b == y) {
        let mut exponents = vec![0; factor_base.len()];
        exponents[pos] = 1;
        return Some(exponents);
    }

    let mut exponents = vec![0; factor_base.len()];
    let mut remaining = y;

    for (i, &prime) in factor_base.iter().enumerate() {
        let mut exponent = 0;
        while remaining > 0 && remaining.is_multiple_of(prime) {
            remaining /= prime;
            exponent += 1;
        }
        exponents[i] = exponent;
    }

    if remaining == 1 {
        Some(exponents)
    } else {
        // remaining != 1 => Some factors of y are not in factor_base
        None
    }
}

/// Attempts to find a non-trivial factor of n using valid combinations.
/// Returns Some(factor) if a non-trivial factor is found, None otherwise.
pub fn find_factor(
    relations: Vec<(u64, u64)>,
    valid_combinations: Vec<Vec<u64>>,
    t: u64,
    n: u64,
) -> Option<u64> {
    for combination in valid_combinations {
        let mut x_product: u64 = 1;
        let mut y_squared_product: u64 = 1;

        for (i, &weight) in combination.iter().enumerate() {
            if weight == 1 {
                let (x_i, q_xi) = relations[i];
                x_product = x_product
                    .checked_mul(x_i + t)
                    .and_then(|p| p.checked_rem(n))?;
                y_squared_product = y_squared_product
                    .checked_mul(q_xi)
                    .and_then(|p| p.checked_rem(n))?;
            }
        }

        if let Some(y) = is_perfect_square(y_squared_product) {
            let y_mod = y % n;

            // Check both (x - y) and (x + y) for non-trivial factors
            for &sign in &[-1, 1] {
                let diff = if sign == -1 {
                    (x_product + n - y_mod) % n // Equivalent to (x - y) mod n
                } else {
                    (x_product + y_mod) % n // Equivalent to (x + y) mod n
                };

                let d = gcd(diff, n);
                if 1 < d && d < n {
                    return Some(d);
                }
            }
        }
    }

    None
}

/// Computes the optimal p_max parameter for the quadratic sieve.
/// alpha: smoothness parameter (typically 1.0)
pub fn compute_p_max(n: u64, alpha: f64) -> u64 {
    let log_n = (n as f64).ln();
    let log_log_n = log_n.ln();
    let l_n = E.powf(alpha * (log_n * log_log_n).sqrt());
    l_n.floor() as u64
}

/// Computes the optimal x_max parameter for the quadratic sieve.
pub fn compute_x_max(n: u64) -> u64 {
    let log_n = (n as f64).ln();
    let log_log_n = log_n.ln();
    let l_n = E.powf((log_n * log_log_n).sqrt());
    ((2 * n) as f64).sqrt() as u64 * l_n.powf(1.0 / 2.0_f64.sqrt()) as u64
}

/// Main quadratic sieve factorization function.
/// Returns Some(factor) if a non-trivial factor is found, None otherwise.
pub fn quadratic_sieve_factorization(
    n: u64,
    opt_p_max: Option<u64>,
    opt_x_max: Option<u64>,
) -> Option<u64> {
    println!("Attempting to factor n = {}", n);
    let length_bits = bit_length(n);
    let length_digits = digit_length(n);
    println!("n has {} digits and {} bits", length_digits, length_bits);

    // Skip trivial cases
    if n.is_multiple_of(2) {
        return Some(2);
    }

    let p_max = opt_p_max.unwrap_or_else(|| compute_p_max(n, 1.0));
    let x_max = opt_x_max.unwrap_or_else(|| compute_x_max(n));

    println!("Using p_max = {}, x_max = {}", p_max, x_max);

    let factor_base = generate_primes_sieve(p_max);
    println!("Factor base: {} primes up to {}", factor_base.len(), p_max);

    let mut relations = Vec::new();
    let mut exponent_matrix = Vec::new();

    let t = (n as f64).sqrt().floor() as u64;
    println!("Searching for relations with t = {}", t);

    // Collect relations
    for x in 1..=x_max {
        let y = quadratic_residue(x, t, n);
        if let Some(exponents) = factor_with_base(y, &factor_base) {
            relations.push((x, y));
            exponent_matrix.push(exponents);
        }
    }

    let relation_count = exponent_matrix.len();
    let matrix_size = relation_count * factor_base.len();

    if matrix_size < 1000 {
        println!("Found {} relations:", relation_count);
        println!("Relations: {:?}", relations);
        println!("Exponent matrix: {:?}", exponent_matrix);
    } else {
        println!(
            "Found {} relations (matrix size: {})",
            relation_count, matrix_size
        );
    }

    // Find valid combinations in the kernel
    let valid_combinations = find_kernel_mod2_gauss_jordan(exponent_matrix);

    if matrix_size < 1000 {
        println!("Valid combinations: {:?}", valid_combinations);
    } else {
        println!("Found {} valid combinations", valid_combinations.len());
    }

    find_factor(relations, valid_combinations, t, n)
}

/// Factorizes an integer `n` using Pollard's Rho algorithm.
pub fn pollard_rho_factorization(n: u128, seed: u64) -> u128 {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed); // Initialize RNG with seed

    loop {
        // Generate a random constant `c` in the range [2, n-1]
        let c: u128 = rng.random_range(2..n);

        // Initialize x, y, and d for the algorithm
        let mut x: u128 = 2;
        let mut y: u128 = 2;
        let mut d: u128 = 1;

        // Loop until a non-trivial factor is found
        while d == 1 {
            // Iterate the function f(x) = (x^2 + c) mod n
            x = (x.wrapping_mul(x).wrapping_add(c)) % n;
            y = (y.wrapping_mul(y).wrapping_add(c)) % n;
            y = (y.wrapping_mul(y).wrapping_add(c)) % n;

            // Compute GCD of |x-y| and n
            d = gcd(x.abs_diff(y), n);
        }

        // Return the non-trivial factor if found
        if d != n {
            return d;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quadratic_residue() {
        assert_eq!(quadratic_residue(1, 3, 10), 16 - 10); // (1+3)^2 - 10 = 16-10 = 6
        assert_eq!(quadratic_residue(2, 3, 10), 25 - 10); // (2+3)^2 - 10 = 25-10 = 15
    }

    #[test]
    fn test_factor_with_base() {
        let factor_base = vec![2, 3, 5];
        assert_eq!(factor_with_base(2, &factor_base), Some(vec![1, 0, 0]));
        assert_eq!(factor_with_base(6, &factor_base), Some(vec![1, 1, 0]));
        assert_eq!(factor_with_base(30, &factor_base), Some(vec![1, 1, 1]));
        assert_eq!(factor_with_base(7, &factor_base), None);
    }

    #[test]
    fn test_small_numbers() {
        assert_eq!(
            quadratic_sieve_factorization(15, Some(5), Some(10)),
            Some(3)
        );
        assert_eq!(
            quadratic_sieve_factorization(35, Some(7), Some(10)),
            Some(5)
        );
        assert_eq!(
            quadratic_sieve_factorization(143, Some(11), Some(20)),
            Some(11)
        );

        assert_eq!(pollard_rho_factorization(15, 42), 3_u128);
        assert_eq!(pollard_rho_factorization(15, 47), 5_u128);
        assert_eq!(pollard_rho_factorization(35, 42), 7_u128);
        assert_eq!(pollard_rho_factorization(143, 42), 13_u128);
        assert_eq!(pollard_rho_factorization(1000730021, 42), 10007_u128);
    }

    #[test]
    fn test_compute_parameters() {
        assert!(compute_p_max(100, 1.0) > 0);
        assert!(compute_x_max(100) > 0);
    }
}
