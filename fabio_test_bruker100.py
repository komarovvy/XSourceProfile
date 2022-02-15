# import numpy as np
# import matplotlib.pyplot as plt
# from scipy import interpolate
# x = np.arange(5, 20)
# y = np.exp(x/3.0)
# f = interpolate.interp1d(x, y)
# x1 = np.arange(6, 19, 0.2)
# y1 = f(x1)
# plt.plot(x, y, 'o', x1, y1, 'x')
# plt.show()

import os
import numpy as np
import matplotlib.pyplot as plt
import pathlib as path


from fabio.bruker100image import bruker100image

from constants import ROI_XY
from constants import SFRM_DIR, SFRM_NAME_PATTERN
#from constants import RUNS, FRAMES
from constants import FRAMES_OF_INTEREST


I = bruker100image()
# print(" *** Read frame: ")
I.read("mo_20211122_StillBall2_01_0001.sfrm")
# print(" *** Show header: ")
I.show_header(0)
# print(" *** Show data: ")
# print(I.data[ROI_XY['yl']:ROI_XY['yh'], ROI_XY['xl']:ROI_XY['xh']])
# print("Sum =", np.sum(I.data[ROI_XY['yl']:ROI_XY['yh'], ROI_XY['xl']:ROI_XY['xh']]))
# print(f'DETTYPE: {I.header["DETTYPE"]}')


# print(" *** Show frame picture: ")
# frm_fig = plt.imshow(np.log(I.data),
#                      extent=[0, I.dim2, 0, I.dim1],
#                      origin="lower")
# frm_fig.axes.set_autoscale_on(False)
# plt.grid()
# plt.scatter(50,70,20,facecolor='none')
# plt.colorbar()




#for run in range(RUNS['min'], RUNS['max']+1):
#    print(' *** run #', run)
#    for frm in range(FRAMES['min'], FRAMES['max']+1, FRAMES['step']):
#        path = os.path.join(SFRM_DIR, SFRM_NAME_PATTERN.format(run=run, frm=frm))
#        print(path)
sfrm_image = bruker100image()

# table_phi_sumI = open("mo_run.txt", "w+")
# table_phi_sumI.write('OM\tPHI\tSUM_I\n')
#
# for run in FRAMES_OF_INTEREST:
#     print(' *** run #', run)
#
#
#     #table_phi_sumI.write(f'run {run}\n')
#     for frm in FRAMES_OF_INTEREST[run]:
#         path = os.path.join(SFRM_NAME_PATTERN.format(run=run, frm=frm))
#         sfrm_image.read(path)
#         sum_I = (np.sum(sfrm_image.data[ROI_XY['yl']:ROI_XY['yh'], ROI_XY['xl']:ROI_XY['xh']]))
#         #print(f"Frm = {frm},  Sum = {sum_I}, \t\t path = {path}")
#         # start value of the changing angle; all 4 ending angles ( 2T, OM, PH, CH); incremente of the changing angle; all 4 starting angles
#         #print(f'{sfrm_image.header["START"]}-{sfrm_image.header["ENDING"]}-{sfrm_image.header["INCREME"]}, {sfrm_image.header["ANGLES"]}')
#         sfrm_angles = sfrm_image.header["ANGLES"].split()
#         #print(f'OM = {sfrm_angles[1]}, PHI (start) = {sfrm_angles[2]}')
#         table_phi_sumI.write(f"{sfrm_angles[1]}\t{sfrm_angles[2]}\t{sum_I}\n")
# table_phi_sumI.close()

class read_run():
    def __init__(self,run, frm):
        path = os.path.join(SFRM_NAME_PATTERN.format(run=run, frm=frm))
        sfrm_image.read(path)
        sum_I = (np.sum(sfrm_image.data[ROI_XY['yl']:ROI_XY['yh'], ROI_XY['xl']:ROI_XY['xh']]))
        sfrm_angles = sfrm_image.header["ANGLES"].split()
        self.om = sfrm_angles[1]
        self.phi = sfrm_angles[2]
        self.sum_I = sum_I


print(read_run(1,30).om)

print(read_run(1,30).phi)
print(float(read_run(1,30).phi))
print(read_run(1,90).sum_I)
