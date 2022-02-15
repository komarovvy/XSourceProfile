import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple

from x_constants import PHASE_TO_MU
from geometry import VecYZ




def cipg_set_int(cell_length=0.01, Imax=1.8e6, sigma=0.055, sigma_cutoff=3., subbining=2):
    '''
    Calculate intensities for Central Isometric Parallel Gaussian X-ray source
     on a source grid.
     
    Parameters
    ----------
    cell_length : float, optional
        Units: mm. Length of a unit cell of the grid. The default is 0.01.
    Imax : float , optional
        Units: count. Amplitude of the I distribution. The default is 1.8e6.
    sigma : float, optional
        Units: mm. Spread of the I distribution. The default is 0.055.
    sigma_cutoff : float, optional
        Units: pcs, to define a size of non-zero I. The default is 3.
    subbining : int, optional
        Units: pcs, subdivision of the grid cells for increasig of I precision. The default is 2.

    Returns
    -------
    tuple of grid_I, grid_dimention
        Integral intensity on the grid.
    '''
    #print(f'The source grid cell length {source.cell_length} mm')
    assert type(subbining) == int
    grid_length = 2 * sigma * sigma_cutoff
    grid_dimention = int(grid_length / cell_length)
    # the grid must have odd number of cells in col/row
    grid_dimention += 1 + grid_dimention % 2
    #print(f'Grid dimentions {grid_dimention}x{grid_dimention}, the length is {grid_length:.3f} mm ({grid_length/sigma:.1f} sigma)')
    subgrid_dimention = grid_dimention * subbining
    #print(f'Subgrid dimentions {subgrid_dimention}x{subgrid_dimention}')
    
    def rsq_ij(i, j, i_centre, cell_len_sq):
        result = (i - i_centre)**2
        result += (j - i_centre)**2
        result *= cell_len_sq
        return result
   
    subcell_area = (cell_length / subbining)**2
    _rsq_ij = lambda i, j: rsq_ij(i, j, i_centre=(subgrid_dimention - 1)/2, 
                                  cell_len_sq=subcell_area)
    
    radii_sq = np.array([_rsq_ij(i,j) for i in range(subgrid_dimention) 
                                      for j in range(subgrid_dimention)])
    #print("Radii array:\n", np.reshape(radii_sq, (subgrid_dimention, subgrid_dimention)))
    
    sigma_sq = sigma**2
    subgrid_norma = 1 / np.pi / sigma_sq
    subgrid_I = subgrid_norma * np.exp(-radii_sq / sigma_sq) * subcell_area
    subgrid_I = np.reshape(subgrid_I, (subgrid_dimention, subgrid_dimention))
    #print(f'Generated subgrid of {subgrid_I.shape} with sum intensity {np.sum(subgrid_I):.3f}')
    
    grid_I = np.array([np.sum(subgrid_I[i : i+subbining, j : j+subbining])
                       for i in range(0, subgrid_dimention, subbining) 
                       for j in range(0, subgrid_dimention, subbining)])
    grid_I *= Imax
    grid_I = np.reshape(grid_I, (grid_dimention, grid_dimention))
    #print(f'Generated grid of {source.grid_I.shape} with sum intensity {np.sum(source.grid_I)}')
    #print('Grid of partial intensities:\n', source.grid_I)

    return grid_I, VecYZ(grid_dimention, grid_dimention)
    
    
    
def table_set_int(source, Iij_2D_arr=(0.,0.,0., 0.,1.8e6,0., 0.,0.,0.)):
    print(f'{SOURCE_TYPE["Tab"]["description"]} source type is not implenented...')


SOURCE_TYPE = {
    'CIPG': {'description': 'Centered isotropic parallel Gaussian', 
             'init_func': cipg_set_int, 
             'param_list': {'Imax', 'sigma', 'sigma_cutoff'},
             },
    'Tab': {'description': 'Individual table values',
            'init_func': table_set_int, 
            'param_list': {'Iij_2D_arr'},
            },
}


class GridXRaySourceProfile():
    def __init__(self, wavelength='MoKaw', cell_len=0.01, center_shift=(0., 0.),
                 source_type='CIPG', 
                 source_param={'Imax': 1.8e6, 'sigma':0.055, 'sigma_cutoff':3.}):
        if source_type not in SOURCE_TYPE:
            raise ValueError(f'Source type "{source_type}" not implemented.')
        self.wavelength = wavelength   
        self.cell_length = cell_len
        if set(source_param.keys()).issubset(SOURCE_TYPE[source_type]['param_list']):
            #!!! should set .grid_I, .grid_len, and .grid_dim properties
            self.grid_I, self.grid_dim =\
                SOURCE_TYPE[source_type]['init_func'](self.cell_length, **source_param)
            self.type = source_type
        else:
            raise ValueError(f'Wrong parameters for "{source_type}" source type.' +
                             f'Should be {SOURCE_TYPE[source_type]["param_list"]}.')
        # shifts are correct to +- 1/2 of a cell of the grid!!!
        if len(center_shift) != 2 or not all([isinstance(i, (int,float)) for i in center_shift]):
            raise ValueError(f'Center shift should be a pair of floats (yz vector in mm)')
        shift_indeces = VecYZ(*center_shift) / self.cell_length
        self.origin_index = round((self.grid_dim - VecYZ(1, 1)) / 2. - shift_indeces)
        #self.origin_index = VecYZ((self.grid_dim.y - 1) // 2, (self.grid_dim.z - 1) // 2)
        
    def coord_to_index(self, r: VecYZ) -> VecYZ:
        i = r / self.cell_length + self.origin_index
        return round(i)
    
    def index_to_c_coord(self, i: VecYZ) -> VecYZ:
        # center of the cell with index i
        r = (i - self.origin_index) * self.cell_length
        return r

    def index_to_l_coord(self, i: VecYZ) -> VecYZ:
        # the lower corner of the cell with index i
        r = (i - self.origin_index - VecYZ(0.5, 0.5)) * self.cell_length
        return r

    def index_to_u_coord(self, i: VecYZ) -> VecYZ:
        # the upper corner of the cell with index i
        r = (i - self.origin_index + VecYZ(0.5, 0.5)) * self.cell_length
        return r
    
    def show_grid(self):
        fig, ax = plt.subplots()
        cs = ax.pcolormesh(self.grid_I)
        cbar = fig.colorbar(cs)
        
        plt.show()
    
    def show_I(self, unit='mm'):
        if unit == 'mm':
            unit_transform_factor = 1.
        elif unit == 'um':
            unit_transform_factor = 1000.
        else:
            raise ValueError(f'Unknown units for length: {unit}')
        fig, ax = plt.subplots()
        ax.invert_xaxis()
        ax.set_aspect('equal', adjustable='box')
        #TODO change when make it non-uniform and shifted!
        r_min = self.index_to_l_coord(VecYZ(0, 0))
        r_max = self.index_to_u_coord(self.grid_dim - VecYZ(1, 1))
        y_range = np.linspace(r_min.y, r_max.y, self.grid_dim.y + 1) * unit_transform_factor
        z_range = np.linspace(r_min.z, r_max.z, self.grid_dim.z + 1) * unit_transform_factor
        cs = ax.pcolormesh(y_range, z_range, self.grid_I)
        plt.scatter(0., 0., color='red', marker= '+')
        plt.scatter(0., 100., color='red', marker= '^')
        plt.scatter(100., 0., color='red', marker= '<')
        cbar = fig.colorbar(cs)
        plt.show()

    
if __name__ == '__main__':
    test_source = GridXRaySourceProfile(cell_len=0.01, center_shift=(-0.10, 0.08), 
                                        source_param={'Imax': 100., 'sigma':0.06, 'sigma_cutoff':3.})
    test_i = VecYZ(18, 18)
    print(f'Indexes {test_i} corresponds to coordinates of center {test_source.index_to_c_coord(test_i)}')
    print(f' with corners at:')
    print(f' {test_source.index_to_l_coord(test_i)} and {test_source.index_to_u_coord(test_i)}')
    
    test_source.show_I('um')
        