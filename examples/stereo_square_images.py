import matplotlib as mpl
from pathlib import Path
from mpl_stereo import StereoSquare2D

curr_dir = Path(__file__).parent

## User inputs
data_dir = curr_dir / '..' / 'src' / 'mpl_stereo' / 'data'
image_left = data_dir / 'church_left.jpg'
image_right = data_dir / 'church_right.jpg'
output = curr_dir / 'output.gif'

## Plot the stereo square
data_left = mpl.image.imread(image_left)
data_right = mpl.image.imread(image_right)
stereosquare = StereoSquare2D()
stereosquare.imshow_stereo(data_left, data_right, crop=True)
stereosquare.fig.set_size_inches(8, 6)
stereosquare.fig.set_dpi(320)

## Save to file
stereosquare.wiggle(output)
print(f'Wrote {output}')
