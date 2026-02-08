pub mod lib_factor;
pub mod lib_gf2;
pub mod lib_utils;

use lib_factor::quadratic_sieve_factorization;
use std::time::Instant;

/// Helper function to find a trivial factor (for display purposes)
fn find_trivial_factor(n: u64) -> Option<u64> {
    (2..)
        .take_while(|&i| i * i <= n)
        .find(|&i| n.is_multiple_of(i))
}

fn main() {
    // Test numbers for factorization
    // Small numbers (easy to factor):
    // let n: u64 = 71 * 59;      // 4189
    // let n: u64 = 31 * 37;      // 1147
    // let n: u64 = 37 * 101;     // 3737
    // let n: u64 = 97 * 101;     // 9797
    // let n: u64 = 3 * 37 * 101; // 11211

    // Medium numbers (require more resources):
    // let n: u64 = 37 * 257;     // 9509
    let n: u64 = 1009 * 10007; // 10106063 (takes 2-3 seconds)

    // Large numbers (require larger parameters):
    // let n: u64 = 101 * 1009;   // 101909
    // let n: u64 = 10007 * 100003; // 1000730021 (very long)

    // Optional parameters for quadratic sieve
    // For small numbers, default values are sufficient:
    let opt_p_max: Option<u64> = None; // None uses automatically calculated value
    // let opt_p_max: Option<u64> = Some(200); // To force a specific value

    // For larger numbers, increasing x_max can help
    // let opt_x_max: Option<u64> = Some(1000); // For medium numbers
    // let opt_x_max: Option<u64> = Some(3310); // For large numbers
    let opt_x_max: Option<u64> = None; // None uses automatically calculated value

    println!(
        "Attempting to factor n = {} ({} × {})",
        n,
        find_trivial_factor(n).unwrap_or(n),
        n / find_trivial_factor(n).unwrap_or(1)
    );

    let start = Instant::now();

    match quadratic_sieve_factorization(n, opt_p_max, opt_x_max) {
        Some(f) => {
            let duration = start.elapsed();
            println!("Found factor f = {} for n = {}", f, n);
            println!("Time elapsed: {:.2?}", duration);

            if n.is_multiple_of(f) {
                let cofactor = n / f;
                println!("Factor is valid: {} × {} = {}", f, cofactor, f * cofactor);
            } else {
                println!("ERROR: Factor is not valid: {} does not divide {}.", f, n);
            }
        }
        None => {
            println!("No factor found.");
        }
    }
}
