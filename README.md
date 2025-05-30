# panlib
The Pancine library consists of over 30 headers and 75 program modules,
with a mixture of C and C++ interfaces.
Many of these routines depend on the Radiance library held in ray/src/common,
which must be included with all builds as a header directory, and linked to librtrad.a
to create a binary.  Most tools further depend on image i/o libraries for JPEG, TIFF,
and OpenEXR.  (Radiance HDR and BMP are included in -lrtrad)
Support for image i/o is decided and set at the main() call-level, offering maximum
flexibility for builds on system with missing libraries or limited i/o requirements.
(E.g., a tool that reads and analyzes images may have no need for image writers at all.)
In general, it is best
to link all libraries statically, as there are too many different versions floating
about and they disagree on what features to support and the details of how to support them.

The C interface prototypes may be found in "pimage.h" and headers it includes.
These functions provide basic image format support (reading and writing) along with
general routines for copying, linking and unlinking image subareas, resampling, colorspace
conversion, histogram computation and matching, weighted averaging, blurring,
image pyramids for fast resampling, general filtering and convolution, warping,
rotation, dilation and erosion, alpha map manipulation, image blending, white-balancing,
bilaterial filtering and denoising, and blending panormas.

The C++ interface prototypes are spread between "panimage.h" which describes the
basic image class, and specialty classes and APIs for various operations supported.
The PanImage class itself encapsulates much of the functionality described in the
C API above, with more convenient reading, writing, and pixel accessors.  Simple
operations such as adding and subtracting images and pixels are supported with
arithmetic operator overloading for this class and the PixelVal type defined in
the "imagio.h" header included by "pimage.h", which in turn is included in
"panimage.h".  None of the types or code in this library use or depend on STL.

Some headers and the functionality they support are listed alphabetically below:

ablockarray.h - efficient block allocator/deallocator template for large arrays of identically-sized chunks

astack.h - block-allocated stack template

cache.h - memory cache handling used throughout panlib to manage RAM footprint

dbase.h - extensive database manager used by Photosphere and some other C++ code

dbhlink.h - header manager used by the database routines (and others)

exif.h - EXIF header i/o

gaussjord.h - matrix equation solver template for square and overdetermined systems of linear equations

pancine.h - inclusive header used by Photosphere routines

panimage.h - described above

pdispimg.h - sophisticated and confusing routines for managing image display in Photosphere

pdraw.h - useful set of routines for drawing lines, circles, and polygons in bitmaps and images

pfeatures.h - a feature-matching and alignment warping method

phdrimg.h - API for merging LDR images into HDR, including alignment, ghost- and flare-removal

photophile.h - classes and APIs specifically tailored to Photosphere

ppano.h - class for joining panoramas

pstrings.h - generally useful routines and classes for char * strings and groups of sorted strings

pthumb.h - Photosphere thumbnail creator and memory/disk cache manager

textform.h - text drawing base class

textmap.h - text class for drawing into bitmaps

tiffin.h - basic TIFF metadata input class

If you have specific questions about the library or its use, please write to Greg Ward at
<GJWard@lbl.gov>.  If you have a bug to report or a contribution to make, you probably
know how to use github better than he does, so go for it.
