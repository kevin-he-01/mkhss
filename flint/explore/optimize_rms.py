from typing import Callable
from sage.all import *
import itertools

# Boolean function: function mapping vector of bits to an integer
# Example: f(x1, x2, x3) = x1 AND (x2 OR x3) + 3
def interpolate_boolean_function(n: int, func: Callable[[list[int]], int]):
    # n variables
    R = PolynomialRing(QQ, names=['x' + str(i) for i in range(n)])
    vars = R.gens()
    truth_table = []
    for input in itertools.product([0, 1], repeat=n):
        # Evaluate the function
        output = func(list(input))
        # Append the input and output to the truth table
        truth_table.append((input, output))
    print("Truth table:")
    for input, output in truth_table:
        print(f"Input: {input}, Output: {output}")
    accumulator_poly = R(0)
    for input, output in truth_table:
        # Create a polynomial for each row of the truth table
        poly = R(1)
        for i in range(n):
            if input[i] == 1:
                poly *= vars[i]
            else:
                poly *= (1 - vars[i])
        # Add the polynomial to the polynomial ring
        accumulator_poly += output * poly
    # Print the polynomial
    print("Polynomial:")
    print(accumulator_poly)
    print(R)

# interpolate_boolean_function(2, lambda x: int(x[0] and x[1]))
# interpolate_boolean_function(2, lambda x: int(x[0] or x[1]))
# interpolate_boolean_function(2, lambda x: int(x[0] ^ x[1]))

# interpolate_boolean_function(4, lambda x: int(x[0] and x[1]) + int(x[2] and x[3]))
# interpolate_boolean_function(4, lambda x: int(x[0] ^ x[1]) + int(x[2] ^ x[3]))
