# Changelog

## Future Work
### Features & Maintenance:
- Better beginner instructions
- Add on more known 2D plotting functions
- Flesh out testing of 3D plotting
- 3D take into account roll angle
- GH release publish to pypi
- AxesAnaglyph.imshow_stereo() nonlinear color mapping methods

### Known issues:
- Inconsistent coloring across redraws if not speficied

----

## [Unreleased]
### Added    
### Changed    
### Removed    

## [0.6.0] - 2024-01-14
### Added    
* Church full color anaglyph example
* Support for full color image anaglyphs, with methods `'dubois'`, `'photoshop'`, `'photoshop2'`
* Support for colormaps in image anaglyphs
### Changed    
* Default `zscale` is now 1/4 x-axis range, not 1
* Fix anaglyph color ordering

## [0.5.0] - 2023-12-28
### Added    
* Consistent `zlim` handling added, with functions `get_zlim`, `set_zlim`, `autoscale_z`, `_calc_bounding_zlim`
* Method `redraw` for side-by-side axes
### Changed    
* Default `zscale` fixed
* `z_scale` and `z_zero` renamed to `zscale` and `zzero`
* `AxesStereo` class renamed to `AxesStereoSideBySide`
* `AxesStereo2DBase` abstract base class added

## [0.4.1] - 2023-12-24
### Changed    
* Fix anaglyph color ordering

## [0.4.0] - 2023-12-23
### Added    
* `AxesAnaglyph.imshow_stereo()` for generating anaglyphs from existing image data
* `mpl_stereo.example_data` module
* `AxesStereo2D` and `AxesStereo3D` keep track of their artists in `.artists_left` and `artists_right`
* Derivation of geometry
### Changed    
* Minimal whitespace between AxesStereo2D subplots

## [0.3.1] - 2023-12-21
### Added    
* Project logo
* Support for 2D plotting methods: `text`
### Changed    
* Docs improvements
* Fix bug when no variation in z data
### Removed    

## [0.3.0] - 2023-12-18
### Added    
* Animation example for 2D stereograms
* `z_zero` keyword argument to 2D plotting methods to set focal plane
### Changed    
* Allow passing in existing axes to object creation
* Rename `focal_plane` to `eye_balance`
* `z_scale` is now a keyword argument to 2D plotting methods to set z-axis scaling, fixed units to match x-axis
* Fixed 3D plotting z scaling
### Removed    
* `z_scale` argument from AxesStereo3D

## [0.2.0] - 2023-12-17
### Added    
* AxesAnaglyph for red-cyan anaglyphs
* Publishing instructions
### Changed    
* Only apply inaccurate axis label transparency to x axis labels
* 3D plots now share their view while rotating, with offset applied
* Everywhere applicable, pass-through methods now return a tuple of (res_left, res_right) for the results from both axes calls

## [0.1.2] - 2023-12-15
### Added
* Tests
### Changed 
* 2D scatter plots have their z values sorted to not improperly occlude   

## [0.1.1] - 2023-12-15
### Added
* README warning about data shifting in 2D plots

## [0.1.0] - 2023-12-15
### Added
* Initial release!
