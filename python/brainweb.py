# -*- coding: utf-8 -*-
"""
Created on Fri May 17 16:34:30 2019

@author: abm15
"""
import numpy as np
import os
from tqdm.auto import trange


def PETbrainWebPhantom(
      phanPath, phantom_number=0,voxel_size=None,image_size=None, num_lesions=10,
      lesion_size_mm=[2, 10], pet_lesion=False, t1_lesion=False, t2_lesion=False,
      hot_cold_ratio=0.5):

    if type(phantom_number)==list:
         pet = np.zeros([len(phantom_number),]+image_size)
         mumap = 0*pet
         t1 = 0*pet
         t2 = 0*pet
         for i in range(len(phantom_number)):
             pet[i,:,:,:], mumap[i,:,:,:], t1[i,:,:,:], t2[i,:,:,:] = PETbrainWebPhantom(phanPath, phantom_number[i], \
                voxel_size,image_size, num_lesions, lesion_size_mm, pet_lesion, t1_lesion, t2_lesion,hot_cold_ratio)
         return pet, mumap, t1, t2

    if voxel_size is None:
        voxel_size = [2.08625, 2.08625, 2.03125]
    if image_size is None:
        image_size = [344,344,127]
    #filename='D:\\pyTorch\\brainweb_20_raws\\subject_04.raws'
    filename = download_brain_web(phanPath, phantom_number)
    if filename.endswith('.gz'):
         import gzip
         file = gzip.open(filename, "r")
         phantom = np.frombuffer(file.read(), dtype='uint16').copy()
    else:
          phantom = np.fromfile(filename, dtype='uint16')
    phantom = phantom.reshape([362, 434, 362]).transpose(1,2,0)
    phantom =phantom[::-1, :, :]

    # PHANTOM PARAMETER
    indicesCsf = phantom == 16
    indicesWhiteMatter = phantom == 48
    indicesGrayMatter = phantom == 32
#    indicesFat = phantom == 64
#    indicesMuscleSkin = phantom == 80
    indicesSkin = phantom == 96
    indicesSkull = phantom == 112
#    indicesGliaMatter = phantom == 128
#    indicesConnectivity = phantom == 144
    indicesMarrow = phantom == 177
    indicesDura = phantom == 161
    indicesBone = indicesSkull | indicesMarrow | indicesDura
    indicesAir  = phantom ==0

    # 0=Background, 1=CSF, 2=Gray Matter, 3=White Matter, 4=Fat, 5=Muscle, 6=Muscle/Skin, 7=Skull, 8=vessels, 9=around fat, 10 =dura matter, 11=bone marrow
    mumap = np.zeros(phantom.shape,dtype='float')
    mu_bone_1_cm = 0.13;
    mu_tissue_1_cm = 0.0975;
    mumap[phantom >0] = mu_tissue_1_cm
    mumap[indicesBone] = mu_bone_1_cm

    #TRANSFORM THE ATANOMY INTO PET SIGNALS
    whiteMatterAct = 32
    grayMatterAct = 96
    skinAct = 16;
    pet = phantom;
    pet[indicesWhiteMatter] = whiteMatterAct
    pet[indicesGrayMatter] = grayMatterAct
    pet[indicesSkin] = skinAct
    pet[~indicesWhiteMatter &  ~indicesGrayMatter & ~indicesSkin] = skinAct/2
    pet[indicesAir] = 0

    # T1
    t1 = 0*phantom;
    whiteMatterT1 = 154
    grayMatterT1 = 106
    skinT1 = 92
    skullT1 = 48
    marrowT1 = 180
    duraT1 = 48
    csfT2 = 48
    t1[indicesWhiteMatter] = whiteMatterT1
    t1[indicesGrayMatter] = grayMatterT1
    t1[indicesSkin] = skinT1
    t1[~indicesWhiteMatter & ~indicesGrayMatter & ~indicesSkin & ~indicesBone] = 0
    t1[indicesSkull] = skullT1
    t1[indicesMarrow] = marrowT1
    t1[indicesBone] = duraT1
    t1[indicesCsf] = csfT2
    # T2
    t2 = 0*phantom;
    whiteMatterT2 = 70;
    grayMatterT2 = 100;
    skinT2 = 70;
    skullT2 = 100;
    marrowT2 = 250;
    csfT2 = 250;
    duraT2 = 200;
    t2[indicesWhiteMatter] = whiteMatterT2
    t2[indicesGrayMatter] = grayMatterT2
    t2[indicesSkin] = skinT2
    t2[~indicesWhiteMatter & ~indicesGrayMatter & ~indicesSkin & ~indicesBone] = 0
    t2[indicesCsf] = csfT2
    t2[indicesSkull] = skullT2
    t2[indicesMarrow] = marrowT2
    t2[indicesBone] = duraT2

    if pet_lesion:
         lesion_pet = random_lesion(indicesWhiteMatter, num_lesions,lesion_size_mm)
         lesion_values = np.zeros(num_lesions)
         indx = list(range(num_lesions))
         np.random.shuffle(indx)
         split = int(np.floor(num_lesions*(hot_cold_ratio)))
         cold_idx,hot_idx = indx[split:],indx[:split]
         lesion_values[hot_idx] = grayMatterAct*1.5
         lesion_values[cold_idx] = whiteMatterAct*0.5
         for le in range(num_lesions):
              pet[lesion_pet[:,:,:,le]] = lesion_values[le]
    if t1_lesion:
         lesion_t1 = random_lesion(indicesWhiteMatter, num_lesions,lesion_size_mm)
         lesion_values = np.zeros(num_lesions)
         indx = list(range(num_lesions))
         np.random.shuffle(indx)
         split = int(np.floor(num_lesions*(hot_cold_ratio)))
         cold_idx,hot_idx = indx[split:],indx[:split]
         lesion_values[hot_idx] = whiteMatterT1*1.5
         lesion_values[cold_idx] = grayMatterT1*0.8
         for le in range(num_lesions):
              t1[lesion_t1[:,:,:,le]] = lesion_values[le]
    if t2_lesion:
         lesion_t2 = random_lesion(indicesWhiteMatter, num_lesions,lesion_size_mm)
         lesion_values = np.zeros(num_lesions)
         indx = list(range(num_lesions))
         np.random.shuffle(indx)
         split = int(np.floor(num_lesions*(hot_cold_ratio)))
         cold_idx,hot_idx = indx[split:],indx[:split]
         lesion_values[hot_idx] = whiteMatterT1*1.5
         lesion_values[cold_idx] = grayMatterT1*0.8
         for le in range(num_lesions):
              t2[lesion_t2[:,:,:,le]] = lesion_values[le]

    pet = regrid(pet,[0.5,0.5,0.5],voxel_size)
    pet = zero_pad(pet,image_size)
    pet[pet<0]=0

    mumap = regrid(mumap,[0.5,0.5,0.5],voxel_size)
    mumap = zero_pad(mumap,image_size)
    mumap[mumap<0]=0

    t1 = regrid(t1,[0.5,0.5,0.5],voxel_size)
    t1 = zero_pad(t1,image_size)
    t1[t1<0]=0

    t2 = regrid(t2,[0.5,0.5,0.5],voxel_size)
    t2 = zero_pad(t2,image_size)
    t2[t2<0]=0

    return pet, mumap, t1, t2

def random_lesion(img, num_lesions,lesion_size_mm):
     import random

     #img : indicesWhiteMatter
     idx = np.array(np.nonzero(img.flatten('F').reshape(-1,))).reshape(-1,)
     i = random.sample(range(0, idx.shape[0]), num_lesions)
     i,j,k=col2ijk(idx[i],img.shape[0],img.shape[1],img.shape[2])

     x = np.arange(img.shape[0])
     y = np.arange(img.shape[1])
     z = np.arange(img.shape[2])
     xx, yy,zz = np.meshgrid(x, y,z,indexing='ij')

     r = lesion_size_mm[0]/0.5+(lesion_size_mm[1]-lesion_size_mm[0])/0.5**np.random.rand(num_lesions)
     lesions = np.zeros((img.shape[0],img.shape[1],img.shape[2],num_lesions))

     for le in trange(num_lesions, unit="lesion"):
          lesions[:,:,:,le] = (((xx-i[le])**2+(yy-j[le])**2+(zz-k[le])**2)<r[le]**2).astype('float')

     return lesions>0

def col2ijk(m,Nx,Ny,Nk):
     n = Nx*Ny
     m+=1
     if np.max(m) > n**2:
        raise ValueError("m is greater than the max number of elements")
     k = np.ceil(m/n)
     temp = (m -(k-1)*n)
     j = np.ceil(temp/Nx)
     i = temp -(j-1)*Nx
     return i.astype(int)-1,j.astype(int)-1,k.astype(int)-1

def regrid(img, voxel_size, new_voxel_size, method='linear'):
    # img (numpy.ndarray)
    # method: linear/ nearest
    from scipy.interpolate import RegularGridInterpolator

    new_shape = [int(np.round(img.shape[k]*voxel_size[k]/new_voxel_size[k])) for k in range(np.ndim(img))]
    if np.ndim(img)==3:
        x, y, z = [voxel_size[k] * np.arange(img.shape[k]) for k in range(3)]
        f = RegularGridInterpolator((x, y, z), img, method=method)
        x, y, z = [np.linspace(0, voxel_size[k] * (img.shape[k]-1), new_shape[k]) for k in range(3)]
        new_grid = np.array(np.meshgrid(x, y, z, indexing='ij'))
        new_grid = np.moveaxis(new_grid, (0, 1, 2, 3), (3, 0, 1, 2))
    else:
        x, y = [voxel_size[k] * np.arange(img.shape[k]) for k in range(2)]
        f = RegularGridInterpolator((x, y), img, method=method)
        x, y = [np.linspace(0, voxel_size[k] * (img.shape[k]-1), new_shape[k]) for k in range(2)]
        new_grid = np.array(np.meshgrid(x, y, indexing='ij'))
        new_grid = np.moveaxis(new_grid, (0, 1, 2), (2, 0, 1))

    new_img = f(new_grid)
    return new_img

def zero_pad(x,new_size):
    #new_size -->tuple
    def idxs(m,q):
        if m <q:
            i = (q-m)/2
            ii=[int(i),int(q-i)]
        else:
            ii=[0,m]
        return ii

    X = np.zeros(new_size,dtype=x.dtype)

    if np.ndim(x)==3 and x.shape[2]>1:
        i,j,k = [idxs(x.shape[k],new_size[k]) for k in range(3)]
        X[i[0]:i[1],j[0]:j[1],k[0]:k[1]] = x
    else:
        i,j = [idxs(x.shape[k],new_size[k]) for k in range(2)]
        X[i[0]:i[1],j[0]:j[1]] = x
    return X

def download_brain_web(phanPath, phantom_number = 0, download_all = False):
     # return file name of phantom_number (0:19), if dosn't exist download it
     links=[
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject04_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D',
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject05_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D',
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject06_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D',
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject18_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D',
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject20_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D',
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject38_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D',
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject41_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D',
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject42_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D',
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject43_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D',
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject44_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D',
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject45_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D',
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject46_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D',
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject47_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D',
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject48_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D',
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject49_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D',
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject50_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D',
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject51_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D',
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject52_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D',
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject53_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D',
     'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?do_download_alias=subject54_crisp&format_value=raw_short&zip_value=gnuzip&download_for_real=%5BStart+download%21%5D']
     if download_all:
          for i in range(len(links)):
               download_brain_web(phanPath, phantom_number = i)
          return
     if phantom_number>19:
          raise ValueError("Choose a phantom number in [0, 19]")
     flname = phanPath + 'brainWeb_subject_'+ str(phantom_number)+'.raws.gz'
     if not os.path.isfile(flname):
        import urllib
        urllib.request.urlretrieve(links[phantom_number], flname)

     return flname

def imRotation(img,angle,num_rand_rotations=0):

     # example: imRotation(img,15), imRotation(img,15,5), imRotation(img,[3,45,-10])
    from scipy.ndimage import rotate

    if num_rand_rotations>0:
        # take angle as an interval
        Angle = 2*angle*np.random.rand(num_rand_rotations)-angle
        if num_rand_rotations>1:
             Angle[0] = 0 # to include no rotation
    else:
         if np.isscalar(angle):
              Angle = [angle]
         else:
              Angle = angle

    imgr = np.zeros((len(Angle),) + img.shape,dtype=img.dtype)
    for a in range(len(Angle)):
         if Angle[a]!=0:
             if np.ndim(img)==3:
                 for i in range(img.shape[2]):
                     imgr[a,:,:,i] = rotate(img[:,:,i],Angle[a],reshape=False,order=1)
             else:
                 imgr[a,:,:] = rotate(img,Angle[a],reshape=False,order=1)
         else:
              if np.ndim(img)==3:
                   imgr[a,:,:,:] = img
              else:
                   imgr[a,:,:] = img
    return np.squeeze(imgr)
