import numpy as np
import matplotlib.pyplot as plt
from x_source import cidg_set_int
from geometry import Vec3D
from sample import BallAbsorber
from  fabio_test_bruker100 import read_run


x_beam = cidg_set_int()
abs_ball = BallAbsorber
shift_beam = np.array([0.,0.])#mm
vec = Vec3D(-0.3, 0., 0.)#mm
print(vec)
rot_vec = vec.rot_euler('XYZ', [-54.74, 0., float(read_run(1,50).phi)])
rot_vec = rot_vec.rot_euler('XYZ', [0., 0., 0.])
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

def compare_int (rns_rng, frm_rng, initial_shift, beam_shift): #shifts are mm
    compare_table = open('compare_book.txt', "w")
    compare_table.write('OM\tPHI\tDIF_SUM_I\n')
    tran_mask = abs_ball().T_mask
    for run in rns_rng:
        for frm in frm_rng:
            rot_vec = initial_shift.rot_euler('XYZ', [-54.74, 0., float(read_run(run, frm).phi)])
            rot_vec = rot_vec.rot_euler('XYZ', [0., 0., float(read_run(run, frm).om)])
            ball_yz_shift = rot_vec.proj_yz
            shift = np.array(ball_yz_shift).tolist() #- beam_shift
            shift[0] = shift[0] * 1000  # to microns
            shift[1] = shift[1] * 1000
            shift = [int(round(shift[0])), int(round(shift[1]))]
            diff_int = read_run(run,frm).sum_I - np.sum(set_experiment(x_beam,tran_mask,shift), dtype=np.float32)
            compare_table.write(f"{read_run(run,frm).om}\t{read_run(run,frm).phi}\t{diff_int}\n")
            print(diff_int, read_run(run, frm).phi, shift)
    print(f"run {run} is done")
    raise compare_table.close()
compare_int(range(1,6,1),range(1,360,1), Vec3D(-0.3 ,0. ,0.), [0,0])



print(read_run(1,91).sum_I)
print(np.sum(set_experiment(x_beam,abs_ball().T_mask,final_shift), dtype=np.float32))
print(read_run(1,91).sum_I - np.sum(set_experiment(x_beam,abs_ball().T_mask,final_shift), dtype=np.float32))

# fig, ax = plt.subplots()
# cs = ax.pcolormesh(np.log(set_experiment(x_beam, abs_ball().T_mask, [0,0])))
# cbar = fig.colorbar(cs)
# plt.show()