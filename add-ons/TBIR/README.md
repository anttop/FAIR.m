# Template-Based Image Reconstruction from Sparse Tomographic Data

## What is it?

This code implements indirect image registration based on [LagLDDMM](https://github.com/C4IR/FAIR.m/tree/master/add-ons/LagLDDMM).

![Template image](/add-ons/TBIR/images/Brain_source.png)
![Unknown image](/add-ons/TBIR/images/Brain_target.png)
![Measurements](/add-ons/TBIR/images/Brain_sino.png)

It is an extension to [FAIR](https://github.com/C4IR/FAIR.m) as described in:

    Lukas F. Lang, Sebastian Neumayer, Ozan Öktem, Carola-Bibiane Schönlieb. Template-Based Image Reconstruction from Sparse Tomographic Data, 2018.

See https://arxiv.org/abs/1810.08596 for the preprint.

If you use this software in your work please cite the abovementioned paper in any resulting publication:

    @techreport{LanNeuOktScho18_report,
      author     = {Lang, L.~F. and Neumayer, S. and Öktem, O. and Schönlieb, C.-B.},
      title      = {Template-Based Image Reconstruction from Sparse Tomographic Data},
      number     = {arXiv:1810.08596},
      numpages   = {26},
      type       = {Preprint on ArXiv},
      url        = {https://arxiv.org/abs/1810.08596},
      year       = {2018}
    }

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

This software was originally written for and tested with MATLAB 2017b.

The following additional libraries are required:

astra-toolbox: Matrix-free implementation of the Radon transform. <br />
GitHub: https://github.com/astra-toolbox/astra-toolbox <br />
URL: https://www.astra-toolbox.com/ <br />
Version used: 3d07f5b <br />

In order to compile it clone the library with

    git clone https://github.com/astra-toolbox/astra-toolbox.git

In order to work with the abovementioned version type:

    cd astra-toolbox
    git checkout 3d07f5b

See https://www.astra-toolbox.com/docs/install.html#for-matlab for instructions.

## Usage

At the moment, only matrix-free operators are supported by this extension.
Moreover, the implementation is restricted to square/cubic geometries.

In order to use GPU support, MATLAB may be required to be started with:

    LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6 \
    LD_LIBRARY_PATH=/home/ll542/store/git/astra/lib \
    matlab-r2017b

To run the code set the path in TBIRstartup.m to where the ASTRA library is located.
Then simply run the startup scripts:

    run('FAIRstartup.m')
    run('TBIRstartup.m')

To run the test cases, execute

    runtests('add-ons/TBIR/tests')

We also provide a script: 

    sh runtests.sh

If no GPU is available all tests will pass but 3D Radon tests will throw warnings.

The figures in the paper were created using the example scripts.
To run all these examples you can use:

    sh runexamples.sh

Examples 6-8 require X-ray tomography data to be downloaded beforehand (see the scripts).
Alternatively, run

    cd add-ons/TBIR/
    downloaddata.sh

to download the data to the folder 'data'.

To reproduce the results for the metamorphosis approach the following code is required:

GitHub: https://github.com/bgris/IndirectMatchingMetamorphosis <br />
Version used: ae0b510 <br />

Results can be reproduced by 

1. Get and install Anaconda from https://www.anaconda.com/.

2. Create environment  

    conda create -c odlgroup -n odl-py35 python=3.5 odl matplotlib pytest scikit-image spyder

2. Activate environment

    source activate odl-py35

3. Install ASTRA toolbox

    conda install -c astra-toolbox astra-toolbox

4. Clone IndirectMatchingMetamorphosis implementation to some location

    git clone https://github.com/bgris/IndirectMatchingMetamorphosis.git

In order to work with the abovementioned version type:

    cd IndirectMatchingMetamorphosis
    git checkout ae0b510

5. Place the files

* metamorphosis_brain.py
* images/Brain_source.png
* images/Brain_target.png

in the directory 'IndirectMatchingMetamorphosis' and adjust the path in
metamorphosis_brain.py.

1. Run

    python3 metamorphosis_brain.py

Results can then be found in the 'results' folder.

## Acknowledgements

Lukas F. Lang and Carola-Bibiane Schönlieb acknowledge support from the Leverhulme Trust project "Breaking the non-convexity barrier", the EPSRC grant EP/M00483X/1, the EPSRC Centre Nr.\ EP/N014588/1, the RISE projects ChiPS and NoMADS, the Cantab Capital Institute for the Mathematics of Information, and the Alan Turing Institute.
Sebastian Neumayer is funded by the German Research Foundation (DFG) within the Research Training  Group  1932,  project  area  P3.
Ozan Öktem is supported by the Swedish Foundation of Strategic Research, grant AM13-0049.
We gratefully acknowledge the support of NVIDIA Corporation with the donation of the Quadro P6000 GPU used for this research.