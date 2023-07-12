from astropy.io import fits

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
    degree = f[1].data['polynomial degree'][0]
    return polyMatrix, degree

def p(u,v, polMatrix, degree):
    psf = zeros((polMatrix.shape[0], polMatrix.shape[1]))
    for i in range(len(polMatrix.shape[0])):
        for j in range(len(polMatrix.shape[1])):
            psf[i,j] = objective_function(polMatrix[i,j,:], u, v, degree)
    return psf
            

