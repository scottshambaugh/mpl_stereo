import matplotlib as mpl
from pathlib import Path
from mpl_stereo import StereoSquare2D

curr_dir = Path(__file__).parent

## User inputs
data_dir = curr_dir
image_left = data_dir / 'gas_tube_holder_left.png'
image_right = data_dir / 'gas_tube_holder_right.png'
output = curr_dir / 'gas_tube_holder_stereo_square.gif'

## Plot the stereo square
data_left = mpl.image.imread(image_left)
data_right = mpl.image.imread(image_right)
stereosquare = StereoSquare2D()
stereosquare.imshow_stereo(data_left, data_right, crop=True)
stereosquare.fig.set_size_inches(6, 6)
stereosquare.fig.set_dpi(320)

## Save to file
stereosquare.wiggle(output)
print(f'Wrote {output}')
