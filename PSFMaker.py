from astropy.io import fits
import numpy as np

def objective_function(p, x, y, degree):
    num_coefficients = (degree + 1) * (degree + 2) // 2
    value = 0
    counter = 0
    for a in range(1, degree + 2):
        for b in range(1, degree + 2):
            if (a - 1) + (b - 1) <= degree:
                value += p[counter] * x**(a - 1) * y**(b - 1)
                counter += 1
    return value

def read_shopt(shoptFile):
    f = fits.open(shoptFile)
    polyMatrix = f[0].data
    degree = f[1].data['POLYNOMIAL_DEGREE'][0]
    s_matrix = f[6].data['s_MATRIX']
    g1_matrix = f[6].data['g1_MATRIX']
    g2_matrix = f[6].data['g2_MATRIX']
    return polyMatrix, degree, s_matrix, g1_matrix, g2_matrix

def p(u,v, polMatrix, degree):
    degree = int(degree)
    psf = np.zeros((polMatrix.shape[0], polMatrix.shape[1]))
    for i in range(polMatrix.shape[0]):
        for j in range(polMatrix.shape[1]):
            psf[i,j] = objective_function(polMatrix[i,j,:], u, v, degree)
    return psf/np.sum(psf)

def analytic_profile(u, v, s_matrix, g1_matrix, g2_matrix, radial_function):
    s = (s_matrix[0] * u**3 + s_matrix[1] * v**3 + s_matrix[2] * u**2 * v + s_matrix[3] * v**2 * u + s_matrix[4] * u**2 + s_matrix[5] * v**2 + s_matrix[6] * u * v + s_matrix[7] * u + s_matrix[8] * v + s_matrix[9])
         
    g1 = (g1_matrix[0] * u**3 + g1_matrix[1] * v**3 + g1_matrix[2] * u**2 * v + g1_matrix[3] * v**2 * u + g1_matrix[4] * u**2 + g1_matrix[5] * v**2 + g1_matrix[6] * u * v + g1_matrix[7] * u + g1_matrix[8] * v + g1_matrix[9])
          
    g2 = (g2_matrix[0] * u**3 + g2_matrix[1] * v**3 + g2_matrix[2] * u**2 * v + g2_matrix[3] * v**2 * u + g2_matrix[4] * u**2 + g2_matrix[5] * v**2 + g2_matrix[6] * u * v + g2_matrix[7] * u + g2_matrix[8] * v + g2_matrix[9])
          
    return radial_function(s, g1, g2)

