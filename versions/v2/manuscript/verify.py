from fractions import Fraction
from itertools import combinations, permutations
import json


def matrix(rows):
    return [[Fraction(entry) for entry in row] for row in rows]


def multiply(left, right):
    return [[sum(left[row][inner] * right[inner][column]
                 for inner in range(len(right)))
             for column in range(len(right[0]))] for row in range(len(left))]


def trace(value):
    return sum(value[index][index] for index in range(len(value)))


def kronecker(left, right):
    return [[left[row // len(right)][column // len(right[0])]
             * right[row % len(right)][column % len(right[0])]
             for column in range(len(left[0]) * len(right[0]))]
            for row in range(len(left) * len(right))]


def determinant(value):
    result = Fraction(0)
    for ordering in permutations(range(len(value))):
        inversions = sum(ordering[left] > ordering[right]
                         for left in range(len(value))
                         for right in range(left + 1, len(value)))
        product = Fraction((-1) ** inversions)
        for row, column in enumerate(ordering):
            product *= value[row][column]
        result += product
    return result


def identity(size):
    return matrix([[int(row == column) for column in range(size)]
                   for row in range(size)])


def flag(index):
    return matrix([[int(row == column == index) for column in range(4)]
                   for row in range(4)])


def run():
    projectors = [matrix([[1, 0], [0, 0]]),
                  matrix([['1/5', '2/5'], ['2/5', '4/5']]),
                  matrix([['1/5', '-2/5'], ['-2/5', '4/5']])]
    for projector in projectors:
        assert multiply(projector, projector) == projector
        assert trace(projector) == 1
        assert projector[0][1] == projector[1][0]
    noise = Fraction(1, 1000)
    pure_weight = (1 - noise) / 3
    products = [kronecker(projector, flag(index))
                for index, projector in enumerate(projectors)]
    state = [[sum(pure_weight * product[row][column] for product in products)
              + noise / 8 * int(row == column)
              for column in range(8)] for row in range(8)]
    assert trace(state) == 1
    alice = [identity(2), matrix([[1, 0], [0, -1]]), matrix([[0, 1], [1, 0]])]
    coefficients = [[trace(multiply(state, kronecker(observable, flag(index))))
                     for index in range(4)] for observable in alice]
    coefficients.append([Fraction(0)] * 4)
    expected = matrix([
        ['1333/4000', '1333/4000', '1333/4000', '1/4000'],
        ['333/1000', '-999/5000', '-999/5000', 0],
        [0, '333/1250', '-333/1250', 0], [0, 0, 0, 0]])
    assert coefficients == expected
    for bob_row in range(4):
        for bob_column in range(4):
            if bob_row != bob_column:
                unit = matrix([[int(row == bob_row and column == bob_column)
                                for column in range(4)] for row in range(4)])
                for observable in alice:
                    assert trace(multiply(state, kronecker(observable, unit))) == 0
    for bob_index in range(4):
        block = [[state[row * 4 + bob_index][column * 4 + bob_index]
                  for column in range(2)] for row in range(2)]
        assert block[0][1] == block[1][0]
    full_coefficients = [row + [Fraction(0)] * 12 for row in coefficients]
    minors = []
    nonzero_rows = set()
    for rows in combinations(range(4), 3):
        for columns in combinations(range(16), 3):
            minor = determinant([[full_coefficients[row][column] for column in columns]
                                 for row in rows])
            if minor:
                minors.append(minor)
                nonzero_rows.add(rows)
    assert nonzero_rows == {(0, 1, 2)}
    expected_minors = [Fraction(147815037, 1562500000),
                       Fraction(110889, 5000000000),
                       -Fraction(110889, 5000000000),
                       Fraction(332667, 12500000000)]
    assert minors == expected_minors
    norm_squared = sum(minor * minor / 8 for minor in minors)
    assert norm_squared == Fraction(559341805929555807, 500000000000000000000)
    bound_squared = Fraction(9, 8192)
    assert norm_squared > bound_squared
    assert minors[0] > Fraction(3, 32)
    endpoint = Fraction(1, 200)
    margin = 512 * (4 - endpoint) * (1 - endpoint) ** 2 - 2025
    assert margin == Fraction(574, 15625) > 0
    tetrahedral_signs = matrix([[1, 1, 1], [1, -1, -1],
                               [-1, 1, -1], [-1, -1, 1]])
    for axis in range(3):
        assert sum(signs[axis] for signs in tetrahedral_signs) == 0
        for other in range(3):
            assert sum(signs[axis] * signs[other] for signs in tetrahedral_signs) == 4 * int(axis == other)
    fourth_order_ratios = {}
    for dimension in range(5, 9):
        for parameter in [Fraction(0), Fraction(1, 1000)]:
            identity_row = [(1 - parameter) / 4 + parameter / dimension] * 4
            identity_row += [parameter / dimension] * (dimension - 4)
            pauli_rows = [[(1 - parameter) * signs[axis] / 4 for signs in tetrahedral_signs]
                          + [Fraction(0)] * (dimension - 4) for axis in range(3)]
            identity_energy = sum(entry * entry for entry in identity_row) / 2
            assert identity_energy == (1 - parameter) ** 2 / 8 + (2 * parameter - parameter ** 2) / (2 * dimension)
            pauli_energies = [sum(entry * entry for entry in row) / 6 for row in pauli_rows]
            assert pauli_energies == [(1 - parameter) ** 2 / 24] * 3
            for row in pauli_rows:
                assert sum(left * right for left, right in zip(identity_row, row)) == 0
            for axis in range(3):
                for other in range(axis + 1, 3):
                    assert sum(left * right for left, right in zip(pauli_rows[axis], pauli_rows[other])) == 0
            fourth_order_squared = identity_energy
            for energy in pauli_energies:
                fourth_order_squared *= energy
            threshold_squared = Fraction((dimension - 1) ** 3, 11664 * dimension ** 4)
            assert fourth_order_squared > threshold_squared
            fourth_order_ratios[f'D={dimension},t={parameter}'] = format(float(fourth_order_squared / threshold_squared), '.10f')
    uniform_ratio = Fraction(124, 125 * 110592) / Fraction(64, 11664 * 625)
    assert uniform_ratio == Fraction(4185, 4096) > 1
    sharp_bound_squared = Fraction(1, 110592)
    sharp_dimensions = []
    for dimension in range(4, 13):
        gram_diagonal = [Fraction(1, 8)] + [Fraction(1, 24)] * 3
        weighted_gram = [gram_diagonal[0]] + [3 * value for value in gram_diagonal[1:]]
        assert weighted_gram == [Fraction(1, 8)] * 4
        gram_determinant = Fraction(1)
        weighted_determinant = Fraction(1)
        for value, weighted in zip(gram_diagonal, weighted_gram):
            gram_determinant *= value
            weighted_determinant *= weighted
        assert gram_determinant == sharp_bound_squared
        assert weighted_determinant == 27 * sharp_bound_squared == Fraction(1, 4096)
        for index in range(dimension):
            bob_marginal = Fraction(1, 4) if index < 4 else Fraction(0)
            assert (bob_marginal == Fraction(1, dimension)) == (dimension == 4)
        sharp_dimensions.append(dimension)
    printed_tail_bound = Fraction(3, 8)
    actual_tail_squared = Fraction(9, 24)
    assert actual_tail_squared > printed_tail_bound ** 2
    assert actual_tail_squared == Fraction(1, 2) * Fraction(3, 4)
    coherent_qubit_cardinalities = []
    for cardinality in range(2, 13):
        overlap = Fraction(cardinality - 2, 2 * (cardinality - 1))
        gram = [[Fraction(1) if row == column else overlap
                 for column in range(cardinality)] for row in range(cardinality)]
        eigenvalue_ones = 1 + (cardinality - 1) * overlap
        eigenvalue_orthogonal = 1 - overlap
        assert eigenvalue_ones > 0 and eigenvalue_orthogonal > 0
        assert sum(gram[0]) == eigenvalue_ones
        if cardinality <= 4:
            coherent_qubit_cardinalities.append(cardinality)
    assert max(coherent_qubit_cardinalities) == 4
    print(json.dumps({
        'method': 'exact rational matrix arithmetic',
        'state': 'rho_(1/1000), dimensions 2 by 4',
        'projectors': 'three symmetric trace-one idempotents',
        'density_trace': str(trace(state)),
        'correlation_block': 'all entries match Equation (3)',
        'compound_minors_checked': 2240,
        'compound_nonzero_entries': len(minors),
        'compound_nonzero_rows': len(nonzero_rows),
        'single_minor_exceeds_bound': True,
        'squared_norm_to_squared_bound_ratio': format(float(norm_squared / bound_squared), '.12f'),
        'noise_interval_endpoint': str(endpoint),
        'endpoint_polynomial_margin': str(margin),
        'tetrahedral_fourth_order_ratios': fourth_order_ratios,
        'tetrahedral_uniform_squared_ratio_lower_bound': str(uniform_ratio),
        'sharp_fourth_order_bound_squared': str(sharp_bound_squared),
        'sharp_equality_dimensions_checked': sharp_dimensions,
        'weighted_gram_diagonal_at_equality': ['1/8'] * 4,
        'appendix_D31_printed_tail_bound': str(printed_tail_bound),
        'appendix_D31_actual_tail_squared': str(actual_tail_squared),
        'appendix_D31_corrected_tail_bound_squared': '3/8',
        'coherent_qubit_design_cardinality_bound': 4,
        'coverage': 'finite checks for the fixed state and endpoint; the interval follows from the written monotonicity argument'
    }, indent=2))


if __name__ == '__main__':
    run()
