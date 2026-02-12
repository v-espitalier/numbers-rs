use std::cmp::PartialEq;
use std::ops::Rem;

/// Calculates the number of bits required to represent a number in binary.
/// Returns 0 for n = 0.
pub fn bit_length(n: u64) -> u32 {
    if n == 0 { 0 } else { 64 - n.leading_zeros() }
}

/// Calculates the number of decimal digits required to represent a number.
/// Returns 0 for n = 0.
pub fn digit_length(n: u64) -> u32 {
    n.checked_ilog10().map_or(0, |x| x + 1)
}

/// Generates all prime numbers up to p_max using the Sieve of Eratosthenes.
/// Returns an empty vector if p_max < 2.
pub fn generate_primes_sieve(p_max: u64) -> Vec<u64> {
    if p_max < 2 {
        return Vec::new();
    }
    if p_max == 2 {
        return vec![2];
    }

    let mut prime_list = vec![2];
    for i in (3..=p_max).step_by(2) {
        let sqrt_i = (i as f64).sqrt() as u64;

        if !prime_list.iter().any(|&j| j <= sqrt_i && i % j == 0) {
            prime_list.push(i);
        }
    }

    prime_list
}

/// Checks if a number is a perfect square.
/// Returns Some(y) if n = y^2, None otherwise.
pub fn is_perfect_square(n: u64) -> Option<u64> {
    if n == 0 {
        return Some(0);
    }

    // Binary search
    let mut left = 1;
    let mut right = n;

    while left <= right {
        let mid = left + (right - left) / 2;
        let square = mid.checked_mul(mid)?;

        match square.cmp(&n) {
            std::cmp::Ordering::Equal => return Some(mid),
            std::cmp::Ordering::Less => left = mid + 1,
            std::cmp::Ordering::Greater => right = mid - 1,
        }
    }

    None
}

/// Computes the greatest common divisor of two numbers using Euclid's algorithm.
/// Returns 0 if both inputs are 0.
#[inline]
pub fn gcd<T>(mut a: T, mut b: T) -> T
where
    T: PartialEq + Rem<Output = T> + From<u8> + Copy,
{
    let zero = T::from(0);
    while b != zero {
        let temp = b;
        b = a % b;
        a = temp;
    }
    a
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bit_length() {
        assert_eq!(bit_length(0), 0);
        assert_eq!(bit_length(1), 1);
        assert_eq!(bit_length(2), 2);
        assert_eq!(bit_length(255), 8);
        assert_eq!(bit_length(u64::MAX), 64);
    }

    #[test]
    fn test_digit_length() {
        assert_eq!(digit_length(0), 0);
        assert_eq!(digit_length(1), 1);
        assert_eq!(digit_length(9), 1);
        assert_eq!(digit_length(10), 2);
        assert_eq!(digit_length(999), 3);
        assert_eq!(digit_length(1000), 4);
        assert_eq!(digit_length(u64::MAX), 20);
    }

    #[test]
    fn test_generate_primes_sieve() {
        assert_eq!(generate_primes_sieve(1), vec![]);
        assert_eq!(generate_primes_sieve(2), vec![2]);
        assert_eq!(generate_primes_sieve(10), vec![2, 3, 5, 7]);
        assert_eq!(
            generate_primes_sieve(30),
            vec![2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        );
    }

    #[test]
    fn test_is_perfect_square() {
        assert_eq!(is_perfect_square(0), Some(0));
        assert_eq!(is_perfect_square(1), Some(1));
        assert_eq!(is_perfect_square(4), Some(2));
        assert_eq!(is_perfect_square(144), Some(12));
        assert_eq!(is_perfect_square(10_000), Some(100));
        assert_eq!(is_perfect_square(4_294_836_225), Some(65_535));
        assert_eq!(is_perfect_square(u64::MAX - 1), None);

        assert_eq!(is_perfect_square(2), None);
        assert_eq!(is_perfect_square(3), None);
        assert_eq!(is_perfect_square(10), None);
        assert_eq!(is_perfect_square(1000), None);
    }

    #[test]
    fn test_gcd() {
        assert_eq!(gcd(48, 18), 6);
        assert_eq!(gcd(18, 48), 6); // Commutativity
        assert_eq!(gcd(0, 5), 5);
        assert_eq!(gcd(5, 0), 5);
        assert_eq!(gcd(0, 0), 0);
        assert_eq!(gcd(1, 1), 1);
        assert_eq!(gcd(17, 5), 1);
        assert_eq!(gcd(100, 75), 25);
        assert_eq!(gcd(987654321, 123456789), 9);
        assert_eq!(gcd(u64::MAX, 1), 1);
        assert_eq!(gcd(u64::MAX, u64::MAX), u64::MAX);
    }
}
