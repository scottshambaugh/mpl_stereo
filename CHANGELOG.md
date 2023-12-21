# Changelog

## Future Work
### Features & Maintenance:
- Better beginner instructions
- Add on more known 2D plotting functions
- Flesh out testing of 3D plotting
- 3D take into account roll angle
- Document derivation of offsets
- GH release publish to pypi
- Raster anaglyphs, see [example in sunpy](https://docs.sunpy.org/en/stable/generated/gallery/showcase/stereoscopic_3d.html#sphx-glr-generated-gallery-showcase-stereoscopic-3d-py)

----

## [Unreleased]
### Added    
### Changed    
### Removed    

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
