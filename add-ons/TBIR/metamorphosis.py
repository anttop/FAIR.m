#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Copyright 2019 Lukas F. Lang and Sebastian Neumayer

This file is part of TBIR.

    TBIR is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    TBIR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TBIR.  If not, see <http://www.gnu.org/licenses/>.

This script reproduces the result obtained by the approach in:

Gris, B. and Chong, C. and Öktem, O.: Image reconstruction through metamorphosis, 2018.
URL: https://arxiv.org/abs/1806.01225
Code: https://github.com/bgris/IndirectMatchingMetamorphosis

"""
import odl
import numpy as np
import imageio
import os
from skimage.measure import compare_ssim as ssim
##%%
namepath= '/home/ll542/store/git/IndirectMatchingMetamorphosis/'

## The parameter for kernel function
sigma = 2.0
name_sigma=str(int(sigma))

niter=200
epsV=0.02
epsZ=0.0002
## Give regularization parameter
lamb = 1e-5
name_lamb='1e_' + str(-int(np.log(lamb)/np.log(10)))
tau = 1
name_tau='1e_' + str(-int(np.log(tau)/np.log(10)))

# Give the number of time points
time_itvs = 20
nb_time_point_int=time_itvs

name_val_template = 'brain-T.png'
name_val = 'brain-R.png'
num_angles = 6
maxiangle = 'pi'
max_angle = 5*np.pi/6
noise_level = 0.0
noi = '0'
min_angle = np.pi/2

name_exp = 'Brain'
path_data = namepath
path_result_init = namepath + '/results/'
path_result = path_result_init + name_exp
if not os.path.exists(path_result):
    os.makedirs(path_result)

# Discrete reconstruction space: discretized functions on the rectangle
name_ground_truth = path_data + name_val
ground_truth = imageio.imread(name_ground_truth)/255
print(np.sum(np.sum(ground_truth)))
rec_space = odl.uniform_discr(
    min_pt=[-16, -16], max_pt=[16, 16], shape = ground_truth.shape,
    dtype='float32', interp='linear')
ground_truth = rec_space.element(ground_truth)

name_template = path_data + name_val_template
template = rec_space.element(imageio.imread(name_template)/255)
print(np.sum(np.sum(template)))

# Give kernel function
def kernel(x):
    scaled = [xi ** 2 / (2 * sigma ** 2) for xi in x]
    return np.exp(-sum(scaled))

## Create forward operator
## Create the uniformly distributed directions
angle_partition = odl.uniform_partition(min_angle, max_angle, num_angles,
                                    nodes_on_bdry=[(True, True)])

## Create 2-D projection domain
## The length should be 1.5 times of that of the reconstruction space
detector_partition = odl.uniform_partition(-24, 24, int(round(rec_space.shape[0]*np.sqrt(2))))

## Create 2-D parallel projection geometry
geometry = odl.tomo.Parallel2dGeometry(angle_partition, detector_partition)

## Ray transform aka forward projection. We use ASTRA CUDA backend.
forward_op = odl.tomo.RayTransform(rec_space, geometry, impl='astra_cpu')


## load data
data_load = forward_op(rec_space.element(ground_truth))

mini= 0
maxi = 1

# show noise-free data
ground_truth.show('noise-free data', clim=[mini, maxi])
# show noisy data
data_load.show('noisy data')

data=[data_load]
#data=[proj_data]
data_time_points=np.array([1])
forward_operators=[forward_op]
Norm=odl.solvers.L2NormSquared(forward_op.range)
Norm_list = [Norm]


##%% Define energy operator
import Metamorphosis as meta
functional=meta.TemporalAttachmentMetamorphosisGeom(nb_time_point_int,
                            lamb,tau,template ,data,
                            data_time_points, forward_operators,Norm_list, kernel,
                            domain=None)


##%% Gradient descent
X_init=functional.domain.zero()

##%%
import Optimizer as opt
X_final = opt.GradientDescent(niter, epsV, epsZ, functional, X_init)

##%% Compute estimated trajectory
image_list_data=functional.ComputeMetamorphosis(X_final[0],X_final[1])


image_list=functional.ComputeMetamorphosisListInt(X_final[0],X_final[1])

deformation_evo=meta.ShootTemplateFromVectorFields(X_final[0], template)

zeta_transp=meta.ShootSourceTermBackwardlist(X_final[0],X_final[1])

image_evol=meta.IntegrateTemplateEvol(template,zeta_transp,0,nb_time_point_int)

ssim_ref = ssim(template, ground_truth)
ssim_res = ssim(deformation_evo[-1], ground_truth)
display(ssim_res)

##%% save results
for i in range(nb_time_point_int + 1):
    tmp = np.array(image_list[i])
    tmp = np.uint8(255*(tmp - np.min(tmp))/np.ptp(tmp))
    imageio.imwrite(path_result + '/Image_t_' + str(i) + '.png', tmp)
    tmp = np.array(image_evol[i])
    tmp = np.uint8(255*(tmp - np.min(tmp))/np.ptp(tmp))
    imageio.imwrite(path_result + '/TemplatePart_t_' + str(i) + '.png', tmp)
    tmp = np.array(deformation_evo[i])
    tmp = np.uint8(255*(tmp - np.min(tmp))/np.ptp(tmp))
    imageio.imwrite(path_result + '/DeformationPart_t_' + str(i) + '.png', tmp)
    tmp = np.array(zeta_transp[i])
    tmp = np.uint8(255*(tmp - np.min(tmp))/np.ptp(tmp))
    imageio.imwrite(path_result + '/SourcePart_t_' + str(i) + '.png', tmp)
#
