import numpy as np
import matplotlib.pyplot as plt
from x_source import cidg_set_int
from geometry import Vec3D
from sample import BallAbsorber

x_beam = cidg_set_int()
abs_ball = BallAbsorber
shift_beam = np.array([0.,0.])#mm
vec = Vec3D(0.3, 0., 0.)#mm
print(vec)
rot_vec = vec.rot_euler('xyz', [0., 54.74, 0.])
rot_vec = rot_vec.rot_euler('xyz', [0., 0., 0.])
rot_vec = rot_vec.rot_euler('xyz', [0., 0., 0.])
final_shift = -shift_beam + rot_vec.proj_yz
final_shift = final_shift * 1000 #to microns
final_shift = [int(round(final_shift[0])), int(round(final_shift[1]))]
print(final_shift)

#final_shift = np.array([-100,-100])


if final_shift[0] == 0 and final_shift[1] == 0:
    b = x_beam * abs_ball().T_mask
elif final_shift[0] >= x_beam.shape[0] or final_shift[1] >= x_beam.shape[0]:
    b = x_beam
elif final_shift[0] > 0 and final_shift[1] == 0:
    b = x_beam * np.hstack([np.ones((x_beam.shape[0], final_shift[0])),abs_ball().T_mask[:, :-final_shift[0]]])
elif final_shift[0] < 0 and final_shift[1] == 0:
    b = x_beam * np.hstack([abs_ball().T_mask[:, x_beam.shape[0]+final_shift[0]:], np.ones((x_beam.shape[0], x_beam.shape[0]+final_shift[0]))])
elif final_shift[0] == 0 and final_shift[1] > 0:
    b = x_beam * np.vstack([abs_ball().T_mask[final_shift[1]:, ], np.ones((final_shift[1], x_beam.shape[0]))])
elif final_shift[0] == 0 and final_shift[1] < 0:
    b = x_beam * np.vstack([np.ones((final_shift[1] + x_beam.shape[0], x_beam.shape[0])), abs_ball().T_mask[:-final_shift[1], ]])
elif final_shift[0] > 0 and final_shift[1] > 0:
    b = x_beam * np.vstack([np.hstack([np.ones((x_beam.shape[0]-final_shift[1],final_shift[0])), abs_ball().T_mask[final_shift[1]:, :-final_shift[0]]]), np.ones((final_shift[1],x_beam.shape[0]))])
elif final_shift[0] < 0 and final_shift[1] > 0:
    b = x_beam * np.vstack([np.hstack([abs_ball().T_mask[final_shift[1]:, -final_shift[0]:], np.ones((x_beam.shape[0] - final_shift[1], -final_shift[0]))]), np.ones((final_shift[1], x_beam.shape[0]))])
elif final_shift[0] < 0 and final_shift[1] < 0:
    b = x_beam * np.vstack([np.ones((-final_shift[1], x_beam.shape[0])), np.hstack([abs_ball().T_mask[:final_shift[1], -final_shift[0]:], np.ones((x_beam.shape[0] + final_shift[1], -final_shift[0]))])])
elif final_shift[0] > 0 and final_shift[1] < 0:
    b = x_beam * np.vstack([np.ones((-final_shift[1], x_beam.shape[0])), np.hstack([np.ones((x_beam.shape[0] + final_shift[1], final_shift[0])), abs_ball().T_mask[:final_shift[1], :-final_shift[0]]])])
else:
    raise ValueError('shift must be [y,z]')
    pass
#b = np.ones((x_beam.shape[0] + final_shift[1], final_shift[0]))
#b = abs_ball().T_mask[:final_shift[1], :-final_shift[0]]
fig, ax = plt.subplots()
cs = ax.pcolormesh(np.log(b))
cbar = fig.colorbar(cs)
plt.show()
