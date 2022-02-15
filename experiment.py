import numpy as np
import matplotlib.pyplot as plt
from x_source import cidg_set_int
from geometry import Vec3D
from sample import BallAbsorber
from  fabio_test_bruker100 import read_run
import time

x_beam = np.around(cidg_set_int())
abs_ball = BallAbsorber().T_mask
shift_beam = np.array([0.10,-0.25])#mm
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
    elif abs(shift[0]) >= x_beam.shape[0] or abs(shift[1]) >= x_beam.shape[0]:
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
        raise x_beam.shape
    return b
def set_experiment_test (x_beam,transm_mask,shift):
    if shift[0] == 0 and shift[1] == 0:
        x_beam = x_beam * transm_mask
    elif abs(shift[1]) > (transm_mask.shape)[1] or abs(shift[1]) > (transm_mask.shape)[0] or abs(shift[0]) > (transm_mask.shape)[1] or abs(shift[0]) > (transm_mask.shape)[0]:
        x_beam = x_beam
    elif shift[0] == 0 and shift[1] > 0:
        x_beam[shift[1]:,:] = x_beam[shift[1]:,:] * transm_mask[:-shift[1], :]
    elif shift[0] == 0 and shift[1] < 0:
        x_beam[:shift[1], :] = x_beam[:shift[1], :] * transm_mask[-shift[1]:, :]
    elif shift[0] > 0 and shift[1] == 0:
        x_beam[:, shift[0]:] = x_beam[:, shift[0]:] * transm_mask[:,:-shift[0]]
    elif shift[0] < 0 and shift[1] == 0:
        x_beam[:, :shift[0]] = x_beam[:, :shift[0]] * transm_mask[:, -shift[0]:]
    elif shift[0] > 0 and shift[1] > 0:
        x_beam[shift[1]:,shift[0]:] = x_beam[shift[1]:,shift[0]:] * transm_mask[:-shift[1], :-shift[0]]
    elif shift[0] < 0 and shift[1] < 0:
        x_beam[:shift[1],:shift[0]] = x_beam[:shift[1],:shift[0]] * transm_mask[-shift[1]:, -shift[0]:]
    elif shift[0] < 0 and shift[1] > 0:
        x_beam[shift[1]:,:shift[0]] = x_beam[shift[1]:,:shift[0]] * transm_mask[:-shift[1], -shift[0]:]
    elif shift[0] > 0 and shift[1] < 0:
        x_beam[:shift[1], shift[0]:] = x_beam[:shift[1], shift[0]:] * transm_mask[-shift[1]:, :-shift[0]]
    return x_beam
# fig, ax = plt.subplots()
# cs = ax.pcolormesh(set_experiment_test(x_beam, abs_ball, [100,-10000]))
# cbar = fig.colorbar(cs)
# plt.show()

compare_int_dif_beam_shift = open('compare dif beam shift.txt',"w")
compare_int_dif_beam_shift.write('[y_shift,z_shift]\t sum_dif_int\n')
def compare_int (rns_rng, frm_rng, initial_shift, beam_shift): #shifts are mm
    #compare_table = open('compare_book.txt', "w")
    #compare_table.write('OM\tPHI\tDIF_SUM_I\n')
    tran_mask = abs_ball
    sum_diff_int = 0
    beam_shift[0] = beam_shift[0] / 1000
    beam_shift[1] = beam_shift[1] / 1000
    for run in rns_rng:
        for frm in frm_rng:
            rot_vec = initial_shift.rot_euler('XYZ', [-54.74, 0., float(read_run(run, frm).phi)])
            rot_vec = rot_vec.rot_euler('XYZ', [0., 0., float(read_run(run, frm).om)])
            ball_yz_shift = rot_vec.proj_yz
            shift = (ball_yz_shift - np.array(beam_shift)).tolist()
            shift[0] = shift[0] * 1000  # to microns
            shift[1] = shift[1] * 1000
            shift = [int(round(shift[0])), int(round(shift[1]))]
            diff_int = abs(read_run(run,frm).sum_I - np.sum(set_experiment(x_beam,tran_mask,shift), dtype=np.float32))
            #compare_table.write(f"{read_run(run,frm).om}\t{read_run(run,frm).phi}\t{diff_int}\n")
            #print(diff_int, read_run(run, frm).phi, shift)
            sum_diff_int += diff_int
    print(f"run {run} is done")
    compare_int_dif_beam_shift.write(f"{beam_shift}\t{sum_diff_int}\n")
    #compare_table.close()

#compare_int(range(1,6,1),range(1,360,1), Vec3D(-0.3 ,0. ,0.), [0,0])

for y_shif in range(-25,-14,1):
    print(y_shif)
    for z_shift in range(-25,-24,1):
        compare_int(range(1,6,1),range(30,70,1), Vec3D(-0.3 ,0. ,0.), [y_shif,z_shift])
compare_int_dif_beam_shift.close()
print(time.process_time())

