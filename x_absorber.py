# description of the absorbing sample

import numpy as np
import matplotlib.pyplot as plt

from x_constants import PHASE_TO_MU, MATERIALS, WAVELENGTH
from geometry import Vec3D, RectYZ, EulerAng


#TODO put to appropriate place!
def calc_mu(material, wavelength):
    #'Fe, MoKaw': (36.173*2 + 36.781*1)/3 * 7.86 * 0.1,
    phase = material + ', ' + wavelength
    try:
        result = PHASE_TO_MU[phase]
    except:
        raise ValueError(f'Unknown phase: {phase}')
    return result


#TODO check the disagreement between Py and Excel transmissions!!!
def generate_T_mask(dim, step, mu, R):
    #TODO add binning to get more reliable T values
    print("dim, step, mu, R", dim, step, mu, R)
    
    Rsq = R ** 2
    step_sq = step ** 2

    half = dim / 2
    #TODO make rsq_ij(0,0) == rsq_ij(dim-1,dim-1) ! 
    rsq_ij = lambda i, j: ((i - half)**2 + (j - half)**2) * step_sq
    
    # mask_rij = np.array([rsq_ij(i,j)**0.5 for i in range(dim) 
    #                                       for j in range(dim)])
    # print("Chord mask slice:\n", mask_rij.reshape((dim, dim))[::30,::30])
    
    chord = lambda rsq: 2. * (Rsq - rsq)**0.5 if rsq < Rsq else 0.
    
    mask_ch = np.array([chord(rsq_ij(i,j)) for i in range(dim) 
                                           for j in range(dim)])
    
    return np.exp(-mask_ch * mu).reshape((dim, dim))
    

class BallAbsorber():
    # all distances in mm, angles in radians!
    def __init__(self, shift=(0.e-3, 0.e-3, 0.e-3), diameter=300.e-3, 
                       mu=None, material='Fe', wavelength='MoKaw', 
                       mask_step=10.e-3):
        if all(isinstance(val, (int, float)) for val in shift) and\
           len(shift) == 3:
            self.shift = Vec3D(shift)
            self.rot_shift = None
        else:
            raise ValueError('Sample shift must be a tuple of floats of the length 3')
        #TODO add checkings
        if type(diameter) is int or type(diameter) is float:
            if diameter > 0.:
                self.radius = diameter / 2
            else:
                raise ValueError(f'Ball diameter should be above 0! D={diameter} is given')
        else:
            raise ValueError('Ball diameter should be float!')

        self.material = material
        self.wavelength = wavelength

        if mu:
            self.mu = mu
        elif material in MATERIALS and wavelength in WAVELENGTH:
            self.mu = calc_mu(material, wavelength)
        else:
            self.mu = None
        
        self.mask_step = mask_step
        mask_dim = round(2 * self.radius / self.mask_step)
        self.mask_dim = mask_dim + mask_dim % 2
        
        if self.mu:
            self.T_mask = generate_T_mask(dim=self.mask_dim, step=self.mask_step, 
                                          mu=self.mu, R=self.radius)
        
    def __str__(self):
        info = f"Material: {self.material}, wavelength: {self.wavelength}, mu[mm-1] = {self.mu:.2f}"
        info += f"\nmask array dimentions: {self.mask_dim} x {self.mask_dim}, mask step[mm]: {self.mask_step}"
        info += f"\nradius[mm]: {self.radius}, shift[mm]: {self.shift}"
        return info

    def calc_T_integral(self, rect=RectYZ(yl=0., zl=0., yu=0.1, zu=0.1)):
        pass      
    
    #TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    def rotate_euler(self, eu_ang: EulerAng):
        #calc self.rot_shift from eu_ang and self.shift
        pass
    
    #TODO make i<->r transforms using self.rot_shift value !!!
    # def coord_to_index(self, r: VecYZ) -> VecYZ:
    #     i = r / self.cell_length + self.origin_index
    #     return round(i)
    
    # def index_to_c_coord(self, i: VecYZ) -> VecYZ:
    #     # center of the cell with index i
    #     r = (i - self.origin_index) * self.cell_length
    #     return r

    # def index_to_l_coord(self, i: VecYZ) -> VecYZ:
    #     # the lower corner of the cell with index i
    #     r = (i - self.origin_index - VecYZ(0.5, 0.5)) * self.cell_length
    #     return r

    # def index_to_u_coord(self, i: VecYZ) -> VecYZ:
    #     # the upper corner of the cell with index i
    #     r = (i - self.origin_index + VecYZ(0.5, 0.5)) * self.cell_length
    #     return r
    
    # def show_T(self, unit='mm'):
    #     if unit == 'mm':
    #         unit_transform_factor = 1.
    #     elif unit == 'um':
    #         unit_transform_factor = 1000.
    #     else:
    #         raise ValueError(f'Unknown units for length: {unit}')
    #     fig, ax = plt.subplots()
    #     ax.invert_xaxis()
    #     ax.set_aspect('equal', adjustable='box')
    #     #TODO change when make it non-uniform and shifted!
    #     r_min = self.index_to_l_coord(VecYZ(0, 0))
    #     r_max = self.index_to_u_coord(self.grid_dim - VecYZ(1, 1))
    #     y_range = np.linspace(r_min.y, r_max.y, self.grid_dim.y + 1) * unit_transform_factor
    #     z_range = np.linspace(r_min.z, r_max.z, self.grid_dim.z + 1) * unit_transform_factor
    #     cs = ax.pcolormesh(y_range, z_range, self.grid_I)
    #     plt.scatter(0., 0., color='red', marker= '+')
    #     plt.scatter(0., 100., color='red', marker= '^')
    #     plt.scatter(100., 0., color='red', marker= '<')
    #     cbar = fig.colorbar(cs)
    #     plt.show()
    


if __name__ == '__main__':
    test_ball = BallAbsorber(diameter=380.e-3)
    print(test_ball)
    # print("Transmission mask slice:\n", test_ball.T_mask[::30,::30])
    # print("Transmission log-mask slice:\n", np.log(test_ball.T_mask[::30,::30]))
    fig, ax = plt.subplots()
    cs = ax.pcolormesh(np.log(test_ball.T_mask))
    cbar = fig.colorbar(cs)
    plt.show()
    