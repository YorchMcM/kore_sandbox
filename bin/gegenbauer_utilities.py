import numpy as np
import numpy.linalg as la
import scipy.sparse as ss
import scipy.linalg as las
from scipy import fft
from scipy.special import factorial
from typing import Union, Callable, Any
from warnings import warn

from numpy import pi as PI
TWOPI = 2*PI

'''

########################################################################################################################
########################################################################################################################
##########                                                                                                    ##########
##########                                        CHEBYSHEV COEFFICIENTS                                      ##########
##########                                                                                                    ##########
########################################################################################################################
########################################################################################################################

'''


def node(k: int, n: int, N: int) \
        -> float:
    return np.cos((2 * n + 1) * k * PI / 2 / N)


def compute_coefficients(function: Union[Callable | np.ndarray],
                         truncation_order: int,
                         cutoff_threshold: float = 1e-9,
                         function_arguments: Union[Any | None] = None) \
        -> np.ndarray:

    if type(function) is not np.ndarray:
        if function_arguments is None:
            evals = np.array([function(node(1, n, truncation_order + 1)) for n in np.array(range(truncation_order + 1))])
        else:
            evals = np.array([function(node(1, n, truncation_order + 1), function_arguments) for n in
                              np.array(range(truncation_order + 1))])
    else:
        evals = function.copy()

    coeffs = fft.dct(evals) / (truncation_order + 1)
    coeffs[0] = coeffs[0] / 2.0

    coeffs[abs(coeffs) < cutoff_threshold] = 0.0

    return coeffs.flatten()


def coefficients_of_powers_of_x(power: int,
                                truncation_order: Union[int|None] = None) \
        -> np.ndarray:

    if truncation_order is None:
        truncation_order = power
    elif truncation_order < 0:
        raise ValueError('(coefficients_of_powers_of_x): Negative truncation order. Invalid')

    if truncation_order < power:
        warn('(coefficients_of_powers_of_x): Truncation order below significant number of coefficients. Information '
             'will be lost.')
        to_return = np.zeros(power+1)
    else:
        to_return = np.zeros(truncation_order + 1)

    if power == 0:
        to_return[0] = 1.0
    elif power == 1:
        to_return[:2] = np.array([0.0, 1.0])
    elif power == 2:
        to_return[:3] = np.array([0.5, 0.0, 0.5])
    elif power == 3:
        to_return[:4] = np.array([0.0, 0.75, 0.0, 0.25])
    elif power == 4:
        to_return[:5] = np.array([0.375, 0.0, 0.5, 0.0, 0.125])
    elif power == 5:
        to_return[:6] = np.array([0.0, 0.875, 0.0, 0.3125, 0.0, 0.0625])
    elif power == 6:
        to_return[:7] = np.array([0.4375, 0.0, 0.59375, 0.0, 0.1875, 0.0, 0.03125])
    else:
        warn('(chebyshev_coefficients_of_powers_of_x): Power not available. Returning zeros. '
             'Use compute_chebyshev_coefficients.')

    if truncation_order < power:
        to_return = to_return[:truncation_order+1]
    elif truncation_order > len(to_return)-1:
        to_return = np.pad(to_return, (0, truncation_order - len(to_return) + 1))

    return to_return


def expand_coefficients(coefficients: np.ndarray,
                        x: Union[np.ndarray | None] = None,
                        cutoff_threshold: float = 1e-9) \
        -> np.ndarray:

    if x is None:
        x = np.round(np.arange(-1.0, 1.0 + 0.01, 0.01), 2)
    elif min(x) < -1.0 or max(x) > 1.0:
        # If the provided domain overflows the [-1,1] interval, it is re-scaled and shifted.
        x = 2 * (x - min(x)) / (max(x) - min(x)) - 1.0

    coefficients = coefficients.copy()
    coefficients[abs(coefficients) < cutoff_threshold] = 0.0

    nodes = np.arccos(x)
    polynomial_orders = np.nonzero(coefficients)[0]

    expansion = np.zeros_like(x)
    for order in polynomial_orders:
        expansion = expansion + coefficients[order] * np.cos(order * nodes)

    return expansion


'''

########################################################################################################################
########################################################################################################################
##########                                                                                                    ##########
##########                                          GEGENBAUER OPERATORS                                      ##########
##########                                                                                                    ##########
########################################################################################################################
########################################################################################################################

'''


def chebyshev_multiplication_operator(factor_coefficients: np.ndarray,
                                      truncation_order: Union[int | None] = None,
                                      return_intermediate_matrices: bool = False) \
        -> Union[np.ndarray | tuple]:
    '''
    This function computes the matrix :math:`\\boldsymbol{\\mathcal{M}_0[a]}` that multiplies two functions in
    Chebyshev space. It is the direct implementation of the equation right after Eq.(2.7) of Olver and Townsend (2013).

    :param factor_coefficients: Chebyshev coefficients of factor :math:`a`.
    :param truncation_order: Maximum order present in the expansion of the result. Sets the size of the operator.
    :param return_intermediate_matrices: If True, returns the Toeplitz and almost-Hankel matrices separately,
    as well as whole operator itself. Defaults to False.
    :return:

    '''

    if truncation_order is None:
        truncation_order = len(factor_coefficients) - 1

    # Auxiliary padding. Always good to have enough zeros in case they are needed.
    factor_coefficients = np.pad(factor_coefficients, (0, int(2 * truncation_order + 1) - len(factor_coefficients)))

    # Toeplitz matrix
    rowcol = factor_coefficients[:truncation_order + 1]
    rowcol[0] = 2.0 * rowcol[0]
    toeplitz = las.toeplitz(rowcol)

    # Hankel matrix
    col = factor_coefficients[:truncation_order + 1]
    row = factor_coefficients[truncation_order:int(2 * truncation_order + 1)]
    hankel = las.hankel(col, row)
    hankel[0, :] = np.zeros(truncation_order + 1)

    # Select output
    if return_intermediate_matrices:
        to_return = (toeplitz, hankel, 0.5 * (toeplitz + hankel))
    else:
        to_return = 0.5 * (toeplitz + hankel)

    return to_return


def derivative_operator(derivative_order: int,
                        truncation_order: int) \
        -> ss.csr_matrix:
    '''
    Returns the :math:`\\mathbf{\\mathcal D}_\\lambda` matrix of derivatives in Gegenbauer space. Direct
    implementation of its definition in Olver and Townsend (2013).

    :param derivative_order: The order :math:`\\lambda` of the derivative.
    :param truncation_order: The maximum order of the expansion of the result. Sets the size of the operator.

    '''

    if derivative_order < 0:
        raise ValueError('(derivative_operator): Negative derivative order. Invalid.')

    elif derivative_order > truncation_order:
        warn('(derivative_operator): Derivative order higher than truncation order. Result is the zero array.')
        to_return = np.zeros([truncation_order+1, truncation_order+1])

    elif derivative_order == 0:
        to_return = np.eye(truncation_order+1)

    else:
        diagonal = np.array(range(derivative_order, truncation_order + 1))
        to_return = ss.diags(diagonal, offsets=derivative_order).A
        to_return = 2.0 ** (derivative_order - 1) * factorial(derivative_order - 1) * to_return

    return ss.csr_matrix(to_return)


def rotation_operator(basis_index: int,
                      truncation_order: int) \
        -> ss.csr_matrix:
    '''
    Direct implementation of :math:`\\mathbf{\\mathcal{S}}_\\lambda` from Olver and Townsend(2013). Pre-multiplies a
    vector of coefficients of a function in the :math:`C^{(\\lambda)}` and returns the vector of coefficients in the
    :math:`C^{(\\lambda+1)}` basis.

    :param basis_index: The Gegenbauer family index :math;`\\lambda`.
    :param truncation_order: The maximum order of the expansion of the result. Sets the size of the operator.

    '''

    if basis_index < 0:
        raise ValueError('(rotation_operator): Negative basis order. Invalid.')

    if basis_index == 0:

        diagonal0 = 0.5 * np.ones(truncation_order + 1)
        diagonal0[0] = 1.0
        diagonal2 = -0.5 * np.ones(truncation_order - 1)

    else:

        diagonal0 = basis_index / (basis_index + np.array(range(truncation_order + 1)))
        diagonal2 = -basis_index / (basis_index + np.array(range(2, truncation_order + 1)))

    return ss.diags([diagonal0, diagonal2], offsets=[0, 2], format = 'csr')


def multiplication_operator(basis_index: int,
                            coefficients: np.ndarray,
                            vector_parity: int,
                            truncation_order: Union[int|None] = None) \
        -> Union[np.ndarray|ss.csr_array]:
    '''
    This function computes the :math:`\\mathcal{\\mathbf{M}}_\\lambda` operator, as described by Olver and Townsend (2013).

    - If :math:`\\lambda = 0`, the operator reduces to a Toeplitz + almost Hankel operators. In this implementation, this is performed by calling a function specially written for that.
    - If :math:`\\lambda = 1`, the operator reduces to a Toeplitz + Hankel operators. However, it can be also be written as the sum of Toeplitz operators, each one adding to the previous outside of the first row and column. This is what has been done in this implementation.
    - Otherwise, each term in the operator is given by the expression between Eq.(3.7) and Eq.(3.8) of Olver and Townsend (2013).

    The present implementation also considers the parity of all elements. In the application of interest,
    this operator will be pre-multiplying the :math:`\\lambda`-th derivative of an eigenvector, so the parity of (a)
    the eigenvector itself, (b) the derivative and (c) the function represented in this multiplication operator are
    all taken into account to set certain rows and columns to zero.

    In the cases for :math:`\\lambda = 0,1`, this has been done by creating the full operators with the SciPy
    Toeplitz and Hankel functions, and then setting the appropriate rows and columns to 0. In the general case,
    where the operator is to be filled element by element, the useful elements have been identified at the beginning
    and only those have been computed, skipping all the rest.

    :param basis_index: The index of the Gegenbauer family, :math:`\lambda`.
    :param coefficients: The :math:`\lambda`-th Gegenabuer coefficients of the factor represented in :math:`\mathbf{
        \mathcal{M}[a]}`, :math:`a_k`.
    :param vector_parity: The parity of the eigenvector.
    :param truncation_order: The truncation order, setting the size of the operator. If ``None``, it infers the order
        from the size of ``coefficients``. Defaults to None.

    '''

    if basis_index < 0:
        raise ValueError('(multiplication_operator): Negative basis order. Invalid.')

    if vector_parity not in [-1, 0, 1]:
        raise ValueError('(multiplication_operator): Invalid vector parity. Only allowed is 0, 1 or -1. Got '
                         + str(vector_parity) + '.')

    if truncation_order is None:
        truncation_order = len(coefficients)-1

    if sum(abs(coefficients)) == 0.0:

        to_return = np.zeros([truncation_order+1, truncation_order+1])

    else:

        # The parities of all individual elements. The parity of the factor represented in this multiplication
        # operator is assessed using (the parity of) its last non-zero coefficient.
        last_nonzero = np.nonzero(coefficients)[0][-1]
        derivative_parity = int(1 - 2 * (basis_index % 2))
        function_parity = int(1 - 2 * (last_nonzero % 2))
        operator_parity = vector_parity * derivative_parity * function_parity

        if basis_index == 0: # Toeplitz + almost Hankel

            to_return = chebyshev_multiplication_operator(coefficients, truncation_order)

        elif basis_index == 1: # Toeplitz + Hankel

            # Auxiliary padding. Always good to have enough zeros in case they are needed.
            coefficients = np.pad(coefficients, (0, int(2 * truncation_order + 1) - len(coefficients)))
            to_return = las.toeplitz(coefficients[:truncation_order+1])
            for idx in range(1,truncation_order+1):
                to_return[idx:,idx:] = to_return[idx:,idx:] + las.toeplitz(coefficients[2*idx:truncation_order+idx+1])

        else:

            # In the general case, parity is taken care of at the very beginning.
            if vector_parity * derivative_parity == 1:
                column_range = range(0, truncation_order + 1, 2)
            elif vector_parity * derivative_parity == -1:
                column_range = range(1, truncation_order + 1, 2)
            else:
                column_range = range(truncation_order + 1)

            if operator_parity == 1:
                row_range = range(0, truncation_order + 1, 2)
            elif operator_parity == -1:
                row_range = range(1, truncation_order + 1, 2)
            else:
                row_range = range(truncation_order + 1)

            coefficients = coefficients[:last_nonzero+1]

            # Now we loop over all useful entries of the operator.
            to_return = np.zeros([truncation_order+1,truncation_order+1])
            for row in row_range:
                for col in column_range:

                    # Range of s as defined by Olver and Townsend (2013), between Eq.(3.7) and Eq.(3.8).
                    range_of_s = np.array(range(max(0,col-row),col+1))
                    range_of_p = 2 * range_of_s + row - col
                    # The range of those indices are truncated according to the maximum "significant" order of the
                    # provided factor coefficients.
                    range_of_s = range_of_s[range_of_p < len(coefficients)]
                    range_of_p = range_of_p[range_of_p < len(coefficients)]

                    if len(range_of_s) == 0:
                        to_return[row, col] = 0.0

                    else:
                        array_of_multiplication_coefficients = multiplication_coefficients(basis_index,
                                                                                           range_of_s,
                                                                                           col,
                                                                                           range_of_p)

                        to_return[row,col] = np.sum(coefficients[range_of_p] * array_of_multiplication_coefficients)

        if (basis_index == 0 or basis_index == 1) and (vector_parity != 0):

            # If the Toeplitz + Hankel simplifications were used, the parity is taken care of at the end by removing
            # all unnecessary rows/columns.

            if vector_parity * derivative_parity == 1:
                columns_to_remove = np.array(range(1, truncation_order + 1, 2)).astype(int)
            else:
                columns_to_remove = np.array(range(0, truncation_order + 1, 2)).astype(int)

            if operator_parity == 1:
                rows_to_remove = np.array(range(1, truncation_order + 1, 2)).astype(int)
            else:
                rows_to_remove = np.array(range(0, truncation_order + 1, 2)).astype(int)

            to_return[rows_to_remove, :] = np.zeros([len(rows_to_remove), truncation_order + 1])
            to_return[:, columns_to_remove] = np.zeros([truncation_order + 1, len(columns_to_remove)])

    return ss.csr_matrix(to_return)


def multiplication_coefficients(basis_index: int,
                                range_of_subindex: np.ndarray,
                                idx1: int,
                                range_of_idx2: np.ndarray) \
        -> np.ndarray:
    '''
    Computes all the coefficients :math:`c_s^\\lambda(j,k)` for the matrix operator entry. This is the direct
    implementation of the recursion right below Eq.(3.9) from Olver and Townsend (2013). Note that what is called
    :math:`j` in this function according to the definition of :math:`c_s^\\lambda(j,k)` corresponds to the **column** of
    the matrix, confusingly denoted :math:`k` outside of this function

    :param basis_index: The basis index, :math:`\\lambda`.
    :param range_of_subindex: All :math:`s` involved in the sum.
    :param idx1: The :math:`j` involved in the sum
    :param range_of_idx2: All :math:`k` involved in the sum.
    '''

    to_return = np.zeros(len(range_of_subindex))
    to_return[0] = starting_multiplication_coefficient(basis_index,
                                                       int(range_of_subindex[0]),
                                                       idx1,
                                                       int(range_of_idx2[0]))

    for t in range(1, len(range_of_subindex)):
        s = range_of_subindex[t-1]
        idx2 = range_of_idx2[t-1]

        numerator = ((idx1 + idx2 + basis_index - s) *
                     (basis_index + s) *
                     (2 * basis_index + idx1 + idx2 - s) *
                     (idx2 - s + basis_index))

        denominator = ((idx1 + idx2 + basis_index - s + 1) *
                       (s + 1) * (basis_index + idx1 - s - 1) *
                       (basis_index + idx1 + idx2 - s) *
                       (idx2 - s + 1))

        to_return[t] = to_return[t-1] * numerator / denominator

    return to_return


def starting_multiplication_coefficient(basis_index: int,
                                        subindex: int,
                                        idx1: int,
                                        idx2: int) \
        -> float:
    '''
    Computes the starting :math:`c_s^\\lambda(j,k)` coefficient for the matrix operator entry. This is the direct
    implementation of Eq.(3.9) from Olver and Townsend (2013). Note that what is called :math:`j` in this function
    according to the definition of :math:`c_s^\\lambda(j,k)` corresponds to the **column** of the matrix,
    confusingly denoted :math:`k` outside of this function

    :param basis_index: The basis index, :math:`\\lambda`.
    :param subindex: The :math:`s` subindex.
    :param idx1: The :math:`j` index.
    :param idx2: The :math:`k` index.
    '''

    to_return = (idx1 + idx2 + basis_index - 2*subindex) / (idx1 + idx2 + basis_index - subindex)

    for t in range(subindex):
        to_return = to_return * (basis_index+t) / (1+t)
        to_return = to_return * (2*basis_index+idx1+idx2-2*subindex+t) / (basis_index+idx1+idx2-2*subindex+t)

    for t in range(idx1-subindex):
        to_return = to_return * (basis_index+t) / (1+t)
        to_return = to_return * (idx2-subindex+1+t) / (idx2-subindex+basis_index+t)

    return to_return
