import numpy as np
import matplotlib.pyplot as plt

from x_constants import PHASE_TO_MU



def cipg_set_int(source, Imax=1.8e6, sigma=0.055, sigma_cutoff=3., subbining=2):
    #print(f'The source grid cell length {source.cell_length} mm')
    assert type(subbining) == int
    grid_length = 2 * sigma * sigma_cutoff
    grid_dimention = int(grid_length / source.cell_length)
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
   
    subcell_area = (source.cell_length / subbining)**2
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
    
    source.grid_I = np.array([np.sum(subgrid_I[i : i+subbining, j : j+subbining])
                              for i in range(0, subgrid_dimention, subbining) 
                              for j in range(0, subgrid_dimention, subbining)])
    source.grid_I *= Imax
    source.grid_I = np.reshape(source.grid_I, (grid_dimention, grid_dimention))
    source.grid_len = grid_length
    source.grid_dim = grid_dimention
    #print(f'Generated grid of {source.grid_I.shape} with sum intensity {np.sum(source.grid_I)}')
    #print('Grid of partial intensities:\n', source.grid_I)

    return np.sum(subgrid_I)

    
def cimd_set_int(cell_lenght=0.001, Imax=1.8e6, diam=0.3, subbining=2):
    #grid_dimention=int(1 + diam/cell_lenght % 2)
    grid_dimention = 301
    subgrid_dimention = grid_dimention * subbining
    def rsq_ij(i, j, i_centre, cell_len_sq):
        result = (i - i_centre)**2
        result += (j - i_centre)**2
        result *= cell_len_sq
        return result
    sub_cell_area= (cell_lenght/subbining)**2
    _rsq_ij = lambda i, j: rsq_ij(i, j, i_centre=(subgrid_dimention - 1) / 2,
                                  cell_len_sq=sub_cell_area)
    radii_sq = np.array([_rsq_ij(i, j) for i in range(subgrid_dimention)
                         for j in range(subgrid_dimention)])
    radii_sq = np.reshape(radii_sq, (subgrid_dimention, subgrid_dimention))
    print(radii_sq.shape)
    for x in range(0,len(radii_sq)):
        for y in range (0,len(radii_sq)):
            if radii_sq[x,y] > sub_cell_area * (subgrid_dimention - (subgrid_dimention - 1) / 2)**2:
                radii_sq[x,y]=0
            else:
                radii_sq[x,y] = 1

    radii_sq = np.array([np.sum(radii_sq[i: i + subbining, j: j + subbining])
                              for i in range(0, subgrid_dimention, subbining)
                              for j in range(0, subgrid_dimention, subbining)])
    cell_Int = Imax / sum(radii_sq)
    radii_sq = np.reshape(radii_sq, (grid_dimention, grid_dimention))
    radii_sq = cell_Int * radii_sq
    return radii_sq
# fig, ax = plt.subplots()
# cs = ax.pcolormesh(cimd_set_int())
# cbar = fig.colorbar(cs)
# plt.show()

def cidg_set_int(cell_length=1, Imax=1.8e6, sigmai=50, sigmaj=50, sigma_cutoff=3., subbining=2, cor_cof=0.0):
    # print(f'The source grid cell length {source.cell_length} mm')
    assert type(subbining) == int
    if sigmai >= sigmaj:
        grid_length = 2 * sigmai * sigma_cutoff
    else:
        grid_length = 2 * sigmaj * sigma_cutoff

    grid_dimention = int(grid_length / cell_length)
    grid_dimention += 1 + grid_dimention % 2
    subgrid_dimention = grid_dimention * subbining
    subgrid_I = np.array([(1 / (2 * np.pi * (sigmai * subbining) * (sigmaj * subbining) * (1 - cor_cof ** 2) ** 0.5)) * np.exp((-1 / (2 - 2 * cor_cof ** 2)) * ((i - subgrid_dimention / 2) ** 2 / (sigmai * subbining) ** 2 - cor_cof * abs(i - subgrid_dimention / 2) * abs(j - subgrid_dimention) / (sigmai * subbining) / (sigmaj * subbining) + (j - subgrid_dimention / 2) ** 2 / (sigmaj * subbining) ** 2))
                       for i in range(subgrid_dimention)
                       for j in range(subgrid_dimention)])
    subgrid_I = subgrid_I.reshape(subgrid_dimention, subgrid_dimention)
    grid_I = np.array([np.sum(subgrid_I[i: i + subbining-1, j: j + subbining-1])
                       for i in range(0, subgrid_dimention, subbining)
                       for j in range(0, subgrid_dimention, subbining)])
    i_cell = Imax / sum(grid_I)
    grid_I = grid_I.reshape(grid_dimention, grid_dimention)
    grid_I = grid_I * i_cell
    return grid_I
# fig, ax = plt.subplots()
# cs = ax.pcolormesh(cidg_set_int())
# cbar = fig.colorbar(cs)
# plt.show()

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
    def __init__(self, wavelength='MoKaw', 
                 cell_len=0.01, grid_len=None, grid_dim=None,
                 source_type='CIPG', 
                 source_param={'Imax': 1.8e6, 'sigma':0.055, 'sigma_cutoff':3.}):
        if source_type not in SOURCE_TYPE:
            raise ValueError(f'Source type "{source_type}" not implemented.')
        self.wavelength = wavelength
        
        self.cell_length = cell_len
        if grid_dim:
            # the grid size is defined via number of the grid cells
            # ???
            pass
        elif grid_len:
            # the grid dimentions are defined via size of the whole grid
            # ???
            pass
        else:
            # the grid dimentions are defined by source type
            if set(source_param.keys()).issubset(SOURCE_TYPE[source_type]['param_list']):
                # should set .grid_I, .grid_len, and .grid_dim properties
                SOURCE_TYPE[source_type]['init_func'](self, **source_param)
                self.type = source_type
            else:
                raise ValueError(f'Wrong parameters for "{source_type}" source type.' +
                                 f'Should be {SOURCE_TYPE[source_type]["param_list"]}.')
        self.center_index = (self.grid_dim - 1) // 2
        
    def coord_to_index(self, x):
        i = x / self.cell_length + self.center_index
        return round(i)
    
    def index_to_c_coord(self, i):
        # center of the cell with index i
        x = (i - self.center_index) * self.cell_length
        return x

    def index_to_l_coord(self, i):
        # the lower corner of the cell with index i
        x = (i - self.center_index - 0.5) * self.cell_length
        return x

    def index_to_u_coord(self, i):
        # the upper corner of the cell with index i
        x = (i - self.center_index + 0.5) * self.cell_length
        return x
        
    
    
if __name__ == '__main__':
    test_source = GridXRaySourceProfile()
    test_i = 17, 17 
    lcu = ((test_source.index_to_l_coord(test_i[0]), test_source.index_to_l_coord(test_i[1]),),
           (test_source.index_to_c_coord(test_i[0]), test_source.index_to_c_coord(test_i[1]),),
           (test_source.index_to_u_coord(test_i[0]), test_source.index_to_u_coord(test_i[1]),),
           )
    print(f'Indexes {test_i} corresponds to coordinates of center ({lcu[1][0]:.3f},{lcu[1][0]:.3f})')
    print(f' with the lower corner at ({lcu[0][0]:.3f},{lcu[0][0]:.3f}) and the upper one at ({lcu[2][0]:.3f},{lcu[2][0]:.3f})')
    print(f'Their back-calculated indices are:')
    for x, y in lcu:
        print(f'({test_source.coord_to_index(x)},{test_source.coord_to_index(y)})', end='; ')
    print('\nFinish!')
        