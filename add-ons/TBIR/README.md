# Template-Based Image Reconstruction from Sparse Tomographic Data

## What is it?

This code implements indirect image registration based on [LagLDDMM](https://github.com/C4IR/FAIR.m/tree/master/add-ons/LagLDDMM).

This is a MATLAB implementation of the approach described in:

    Lukas F. Lang, Sebastian Neumayer, Ozan Öktem, Carola-Bibiane Schönlieb. Template-Based Image Reconstruction from Sparse Tomographic Data, 2018.

If you use this software in your work please cite the abovementioned paper in any resulting publication.

BibTeX:

    @techreport{LanNeuOktScho18_report,
      author     = {Lang, L.~F. and Neumayer, S. and Öktem, O. and Schönlieb, C.-B.},
      title      = {Template-Based Image Reconstruction from Sparse Tomographic Data},
      number     = {arXiv:1810.08596},
      numpages   = {26},
      type       = {Preprint on ArXiv},
      url        = {https://arxiv.org/abs/1810.08596},
      year       = {2018}
    }

See https://arxiv.org/abs/1810.08596 for the preprint.

## License & Disclaimer

Copyright 2019 Lukas F. Lang and Sebastian Neumayer

This file is part of TBIR. TBIR is free software: you can redistribute it
and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

TBIR is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
details.

You should have received a copy of the GNU General Public License along
with TBIR.  If not, see <http://www.gnu.org/licenses/>.

For the full license statement see the file LICENSE.

## Requirements

This software was originally written for and tested with MATLAB R2017b.

The following libraries are required:

astra-toolbox: Matrix-free implementation of the Radon transform.
GitHub: https://github.com/astra-toolbox/astra-toolbox
URL: https://www.astra-toolbox.com/
Version used: 3d07f5b

>> git clone https://github.com/astra-toolbox/astra-toolbox

In order to work with the abovementioned version type:

>> cd astra-toolbox
>> git checkout 3d07f5b

% TODO: Add information reg. compiling astra-toolbox, setting up path, GPU, etc.

## Usage

% TODO: Add info that only matrix-free operators are supported at the moment.

In order to use GPU support MATLAB may be required to be started with:

>> LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 \
>> LD_LIBRARY_PATH={LOCATION_OF_ASTRA}/astra/lib \
>> matlab-r2017b

Make sure to download ASTRA Toolbox and to set the path in TBIRstartup.m 
properly. Then simply run the FAIR startup script (FAIRstartup.m) and
subsequently the TBIR startup script.

To run the test cases execute

>> runtests('test')

% TODO: script runtests.sh
% TODO: If no GPU available all tests will pass but 3D Radon tests will throw warnings.

The figures in the paper were created with the following scripts:

- paperfigures.m

## Acknowledgements

We gratefully acknowledge the support of NVIDIA Corporation with the
donation of the Quadro P6000 GPU used for this research.

The artificial brain phantom is taken from:

M. Guerquin-Kern, L. Lejeune, K. P. Pruessmann, and M. Unser. Realistic Analytical Phantoms for Parallel Magnetic Resonance Imaging, IEEE Transactions on Medical Imaging, vol. 31, no. 3, pp. 626-636, March 2012.
