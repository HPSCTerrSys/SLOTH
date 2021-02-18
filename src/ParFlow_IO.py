import numpy as np
from struct import pack, unpack

################################################################################
################################ PFB ###########################################
################################################################################
def write_packed(f, fmt, val):
    f.write(pack(fmt, val))

def create_pfb(filename, var, delta=(1, 1, 1), subgrids=(1, 1, 1)):
    print(var.shape)
    nz, ny, nx = var.shape
    dz, dy, dx = delta
    sz, sy, sx = subgrids

    filepfb = open(filename, 'wb')

    # Write start indices of global domain in x, y, z direction
    write_packed(filepfb, '>d', 0)
    write_packed(filepfb, '>d', 0)
    write_packed(filepfb, '>d', 0)

    # Write number of global gridpoints in x, y, z direction
    write_packed(filepfb, '>i', nx)
    write_packed(filepfb, '>i', ny)
    write_packed(filepfb, '>i', nz)

    # Write delta x, delta y and delta z
    write_packed(filepfb, '>d', dx)
    write_packed(filepfb, '>d', dy)
    write_packed(filepfb, '>d', dz)

    nSubGrid = np.prod(subgrids)

    nnx = int(nx / sx)
    nny = int(ny / sy)
    nnz = int(nz / sz)

    # Write the subgrid grid ID
    write_packed(filepfb, '>i', nSubGrid)

    for iz in np.arange(sz)*nnz:
        for iy in np.arange(sy)*nny:
            for ix in np.arange(sx)*nnx:
                #print(ix,iy,iz, nnx,nny,nnz)
                # Write start indices in x, y, z direction
                write_packed(filepfb, '>i', int(ix))
                write_packed(filepfb, '>i', int(iy))
                write_packed(filepfb, '>i', int(iz))

                # Write number of grid points in x, y and z direction for this subgrid
                write_packed(filepfb, '>i', nnx)
                write_packed(filepfb, '>i', nny)
                write_packed(filepfb, '>i', nnz)

                # Write the relative(to global) grid refinement in this subgrid
                # 0=same resolution as global
                write_packed(filepfb, '>i', 0)
                write_packed(filepfb, '>i', 0)
                write_packed(filepfb, '>i', 0)

                # Assuming the data is stored in 3D array called varArray of global size nz*ny*nx
                fmt = ">%dd" % (nnz*nny*nnx)
                filepfb.write(pack(fmt, *var[iz:iz+nnz,
                                               iy:iy+nny,
                                               ix:ix+nnx].flatten()))

    filepfb.close()

def read_pfb(filename):
	with open(filename, "rb") as f:
		f = open(filename, "rb")

		# read meta informations of datafile
		meta_inf = np.fromfile(f, dtype='>f8', count = 3)
		x1 = meta_inf[0]
		y1 = meta_inf[1]
		z1 = meta_inf[2]


		meta_inf = np.fromfile(f, dtype='>i4', count = 3)
		nx = meta_inf[0]
		ny = meta_inf[1]
		nz = meta_inf[2]
		nn = nx * ny * nz


		meta_inf = np.fromfile(f, dtype='>f8', count = 3)
		dx = meta_inf[0]
		dy = meta_inf[1]
		dz = meta_inf[2]


		meta_inf = np.fromfile(f, dtype='>i4', count = 1)
		nsubgrid = meta_inf[0]

		data =  np.ndarray(shape=(nz,ny,nx), dtype='>f8')

		for s in range(nsubgrid):

			meta_inf = np.fromfile(f, dtype='>i4', count = 9)
			ix = meta_inf[0]
			iy = meta_inf[1]
			iz = meta_inf[2]
			# print("---{0} Start Index (X,Y,Z):".format(s+1), ix, iy, iz)

			nx = meta_inf[3]
			ny = meta_inf[4]
			nz = meta_inf[5]
			nn = nx*ny*nz
			# print("---{0} Dimensions (X,Y,Z):".format(s+1), nx, ny, nz)

			rx = meta_inf[6]
			ry = meta_inf[7]
			rz = meta_inf[8]
			# print("---{0} Offsets (X,Y,Z):".format(s+1), rx, ry, rz)

			tmp_data = np.fromfile(f, dtype='>f8', count=nn).reshape((nz,ny,nx))

			data[iz:iz+nz, iy:iy+ny, ix:ix+nx] = tmp_data

	return data

def vanGenuchten(refP, sSat, sRes, nVanG, aVanG):
    """  Calculates the degree of saturation as a function of the pressure head.


    The degree of saturation is calculated as a function of the pressure head
    according to M. Th. van Genuchten.
    Name: A Closed‐form Equation for Predicting the Hydraulic Conductivity of
            Unsaturated Soils
    DOI: https://doi.org/10.2136/sssaj1980.03615995004400050002x

    Parameters:
    refP:       Pressure head [L]
    sSat:       Relative saturated water content [-]
    sRes:       Relative residual saturation [-]
    nVanG:      Non linearity coefficient of the soil [-]
    aVanG:      Air entry values of the soil [L^-1]

    Returns:
    vanG:

    """
    mVanG = 1 - (1 / nVanG)

    vanG = ( (sSat - sRes) / ( 1 + (aVanG * np.absolute(refP))**(nVanG) )**mVanG ) + sRes

    vanG = np.where(refP>=0., 1., vanG) # avoid unsaturated values where water is ponding

    return vanG

