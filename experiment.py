import numpy as np
import matplotlib.pyplot as plt
from x_source import cidg_set_int
from geometry import Vec3D
from sample import BallAbsorber

x_beam = cidg_set_int()
abs_ball = BallAbsorber
shift_beam = np.array([0.,0.])#mm
vec = Vec3D(-0.3, 0., 0.)#mm
print(vec)
rot_vec = vec.rot_euler('XYZ', [-54.74, 0., 45.])
rot_vec = rot_vec.rot_euler('XYZ', [0., 0., 180.])
rot_vec = rot_vec.rot_euler('XYZ', [0., 0., 0.])
final_shift = -shift_beam + rot_vec.proj_yz
final_shift = final_shift * 1000 #to microns
final_shift = [int(round(final_shift[0])), int(round(final_shift[1]))]
print(rot_vec)
print(final_shift)

#final_shift = np.array([0,-200])

def set_experiment (x_beam, transm_mask, shift):
    if shift[0] == 0 and shift[1] == 0:
        b = x_beam * transm_mask
    elif shift[0] >= x_beam.shape[0] or shift[1] >= x_beam.shape[0]:
        b = x_beam
    elif shift[0] > 0 and shift[1] == 0:
        b = x_beam * np.hstack([np.ones((x_beam.shape[0], shift[0])),transm_mask[:, :-shift[0]]])
    elif shift[0] < 0 and shift[1] == 0:
        b = x_beam * np.hstack([transm_mask[:, -shift[0]:], np.ones((x_beam.shape[0], -shift[0]))])
    elif shift[0] == 0 and shift[1] > 0:
        b = x_beam * np.vstack([transm_mask[shift[1]:, ], np.ones((shift[1], x_beam.shape[0]))])
    elif shift[0] == 0 and shift[1] < 0:
        b = x_beam * np.vstack([np.ones((-shift[1], x_beam.shape[0])), transm_mask[:shift[1], ]])
    elif shift[0] > 0 and shift[1] > 0:
        b = x_beam * np.vstack([np.hstack([np.ones((x_beam.shape[0]-shift[1],shift[0])), transm_mask[shift[1]:, :-shift[0]]]), np.ones((shift[1],x_beam.shape[0]))])
    elif shift[0] < 0 and shift[1] > 0:
        b = x_beam * np.vstack([np.hstack([transm_mask[shift[1]:, -shift[0]:], np.ones((x_beam.shape[0] - shift[1], -shift[0]))]), np.ones((shift[1], x_beam.shape[0]))])
    elif shift[0] < 0 and shift[1] < 0:
        b = x_beam * np.vstack([np.ones((-shift[1], x_beam.shape[0])), np.hstack([transm_mask[:shift[1], -shift[0]:], np.ones((x_beam.shape[0] + shift[1], -shift[0]))])])
    elif shift[0] > 0 and shift[1] < 0:
        b = x_beam * np.vstack([np.ones((-shift[1], x_beam.shape[0])), np.hstack([np.ones((x_beam.shape[0] + shift[1], shift[0])), transm_mask[:shift[1], :-shift[0]]])])
    else:
        raise ValueError('shift must be [y,z]')
    return b



for phi in range(0,360,10):
    vec = Vec3D(-0.3, 0., 0.)  # mm

    rot_vec = vec.rot_euler('XYZ', [-54.74, 0., phi])
    rot_vec = rot_vec.rot_euler('XYZ', [0., 0., 180.])
    final_shift = final_shift * 1000  # to microns
    final_shift = [int(round(final_shift[0])), int(round(final_shift[1]))]
    set_experiment(x_beam,abs_ball().T_mask, final_shift)
    print(final_shift)
    fig, ax = plt.subplots()
    cs = ax.pcolormesh(np.log(set_experiment(x_beam, abs_ball().T_mask, final_shift)))
    cbar = fig.colorbar(cs)
    plt.show()
