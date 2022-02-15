SFRM_DIR = 'drive:\\path\\'
SFRM_NAME_PATTERN = 'mo_20211122_StillBall2_{run:02d}_{frm:04d}.sfrm'

RUNS = {'min': 1,  'max': 1}
FRAMES = {'min': 1,  'max': 360, 'step': 1}

# key = run number, value = iterator of frames in the run
FRAMES_OF_INTEREST = {r: range(1, 361, 1) for r in range(1,5)}
BM_CENTRE = {'height': 520, 'width':388}
ROI_SIZE = {'height': 18, 'width': 18}
FRAME_DIM = {'height': 1024, 'width':768}
ROI_CENTRE = {'height': FRAME_DIM['height']//2, 'width': FRAME_DIM['width']//2}

ROI_XY = {'xl': BM_CENTRE['width'] - ROI_SIZE['width']//2,
          'xh': BM_CENTRE['width'] + ROI_SIZE['width']//2,
          'yl': BM_CENTRE['height'] - ROI_SIZE['height']//2,
          'yh': BM_CENTRE['height'] + ROI_SIZE['height']//2,
          }