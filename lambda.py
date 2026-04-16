from functools import cache

def main():
    print(diff_monomial((4,1,1,1)))
    print(diff_monomial((4,)))
    print(make_admissible((1,5,20)))
    print(max(diff_monomial((1000,3,1,1))))
    print(len(diff_monomial((1000,3,1,1))))

def binom_mod2(n, k):
    """Returns True if n choose k is 1 mod 2, False if 0."""
    if k < 0 or k > n:
        return False
    # Lucas' theorem bitwise trick
    return (n - k) & k == 0

@cache
def reduce_pair(i, j):
    """
    Takes an inadmissible pair (i, j) where 2i < j.
    Returns a set of admissible tuples (the polynomial).
    """
    if (2 * i >= j):
        raise ValueError("pair is already admissible")

    result = set()
    
    n = j - 2*i - 1

    for k in range(0,n // 2 + 1):
        if binom_mod2(n - k - 1, k) and 2*i+1+k >= 0:
            result.add((i+n-k,2*i+1+k))
    
    return result

def diff(k):
    return reduce_pair(-1, k)

@cache
def make_admissible(monomial):
    """
    Takes a single tuple and strictly reduces it to an admissible polynomial.
    Returns a set of tuples.
    """
    # Our worklist starts with the raw monomial
    worklist = {monomial}
    admissible_poly = set()
    
    while worklist:
        # Pop an arbitrary monomial to process
        m = worklist.pop()
        
        # Scan for the first inadmissible pair
        inadmissible_index = -1
        for i in range(len(m) - 1):
            if 2 * m[i] < m[i+1]:
                inadmissible_index = i
                break
                
        if inadmissible_index == -1:
            # No inadmissible pairs found; the monomial is fully admissible.
            # XOR it into the result to handle any late-stage cancellations.
            admissible_poly ^= {m}
        else:
            # We found an inadmissible pair at index 'i'.
            i = inadmissible_index
            left = m[:i]
            right = m[i+2:]
            
            # Get the replacement polynomial for the bad pair
            new_pairs = reduce_pair(m[i], m[i+1])
            
            # Construct the new monomials and handle Modulo 2 arithmetic
            for pair in new_pairs:
                new_m = left + pair + right
                
                # Modulo 2 Cancellation Logic:
                if new_m in admissible_poly:
                    # It already survived previously, so 1 + 1 = 0
                    admissible_poly.remove(new_m)
                elif new_m in worklist:
                    # We generated it twice in the intermediate steps, 1 + 1 = 0
                    worklist.remove(new_m)
                else:
                    # It is completely new, add it to the worklist to be scanned
                    worklist.add(new_m)
                    
    return admissible_poly

@cache
def diff_monomial_raw(monomial):
    """
    Applies the Leibniz rule to a full monomial tuple.
    Returns a set of tuples (a polynomial).
    """
    result = set()
    
    for i in range(len(monomial)):
        prefix = monomial[:i]
        suffix = monomial[i+1:]
        
        # Calculate d() for the specific generator at index i
        d_gen = diff(monomial[i])
        
        # Distribute the prefix and suffix over the resulting polynomial
        for term in d_gen:
            # Concatenate the tuples together
            new_monomial = prefix + term + suffix
            
            # XOR into the result set to handle any modulo 2 cancellations
            result ^= {new_monomial}
            
    return result

def diff_monomial(monomial):
    """
    Computes the true differential of a monomial, returning an admissible polynomial.
    """
    # 1. Get the raw Leibniz rule output (likely inadmissible)
    raw_polynomial = diff_monomial_raw(monomial)
    
    final_admissible_polynomial = set()
    
    # 2. Reduce every term to admissibility
    for raw_m in raw_polynomial:
        final_admissible_polynomial ^= make_admissible(raw_m)
        
    return final_admissible_polynomial

def generate_admissible_basis(max_dim):
    """
    Generates all admissible monomials up to max_dim.
    Returns a list of strictly sorted tuples.
    """
    basis = []

    def dfs(current_monomial, current_sum):
        # If the monomial is not empty, add it to our basis
        if current_monomial:
            basis.append(tuple(current_monomial))
        
        # Determine the upper bound for the next generator
        if not current_monomial:
            # No previous element, bounded only by the remaining dimension
            max_next_val = max_dim
        else:
            # Bounded by both admissibility and the remaining dimension
            last_val = current_monomial[-1]
            max_next_val = min(2 * last_val, max_dim - current_sum)
            
        # Branch out to all valid next generators
        for k in range(1, max_next_val + 1):
            current_monomial.append(k)
            dfs(current_monomial, current_sum + k)
            current_monomial.pop() # Backtrack

    # Start the recursive search with an empty list and sum of 0
    dfs([], 0)
    
    # The Curtis algorithm requires a strict lexicographical ordering.
    # Python's default sorting for lists of tuples does exactly this!
    basis.sort()
    
    return basis

if __name__ == "__main__":
    main()