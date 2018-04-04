#!/usr/bin/env python
"""Usage:
  show [options] <rawbinary>...

Options:
  --log=<lvl>  : WARN(ING)|[default: INFO]|DEBUG
  --cmap=<C>   : [default: hot]
  --noise=<n>  : csv float [default: 0.0]
  --vmin=<v>   : [default: each]|auto|none|<float>
  --vmax=<v>   : [default: each]|auto|none|<float>
  --dpi=<D>    : [default: 120:int]
  --outres=<O>  : [default: mMR]|MR
  --save       : bool
  --sigma=<s>  : [default: 1:float]

Examples:
  show --sigma 1.5 --noise 0.5 --outres MR --save \\
    brainweb.raws/subject_??.raws
"""
from __future__ import division
import numpy as np
from caspyr.scripts.showslice import MainApp
from skimage.transform import resize
from skimage.filters import gaussian
from argopt import argopt
from tqdm import tqdm
import logging

__author__ = "Casper da Costa-Luis"
__date__ = "2018-04"
__licence__ = "MPLv2.0"
__requires__ = ["caspyr>=1.0.0"]


class Act(object):
  """careful: occasionally other bits may be set"""
  background, csf, greyMatter, whiteMatter, fat, muscle, skin, skull, vessels,\
      aroundFat, dura, marrow\
      = [i << 4 for i in range(12)]
  bone = skull | marrow | dura

  @classmethod
  def indices(cls, im, attr):
    if attr == "bone":
      return (cls.indices(im, "skull") +
              cls.indices(im, "marrow") +
              cls.indices(im, "dura") > 0)
    return abs(im - getattr(cls, attr)) < 1


class Pet(Act):
  whiteMatter = 32
  greyMatter = whiteMatter * 4
  skin = whiteMatter // 2
  attrs = ["whiteMatter", "greyMatter", "skin"]


class T1(Act):
  whiteMatter = 154
  greyMatter = 106
  skin = 92
  skull = 48
  marrow = 180
  bone = 48
  csf = 48
  attrs = ["whiteMatter", "greyMatter", "skin", "skull", "marrow", "bone",
           "csf"]


class T2(T1):
  whiteMatter = 70
  greyMatter = 100
  skin = 70
  skull = 100
  marrow = 250
  csf = 250
  bone = 200


mu_bone_1_cm = 0.13
mu_tissue_1_cm = 0.0975


class Res(object):
  mMR = np.array([2.0312, 2.0863, 2.0863])
  MR = np.array([1.0, 1.0, 1.0])
  brainweb = np.array([0.5, 0.5, 0.5])


class Shape(object):
  mMR = np.array([127, 344, 344])
  MR = mMR * Res.mMR / Res.MR
  brainweb = mMR * Res.mMR / Res.brainweb


def getRaw(fname):
  """z, y, x"""
  return np.fromfile(fname, dtype=np.uint16).reshape((362, 434, 362))


def noise(im, n, warn_zero=True, sigma=1):
  """
  @param n  : float, noise fraction (0, inf)
  @param sigma  : float, smoothing of noise component
  @return[out] im  : array like im with +-n *100%im noise added
  """
  log = logging.getLogger(__name__)
  if n < 0:
    raise ValueError("Noise must be positive")
  elif n == 0:
    if warn_zero:
      log.warn("zero noise")
    return im
  r = gaussian(np.random.rand(*im.shape), sigma=sigma, multichannel=False)
  return im * (1 + n * (2 * r - 1))


def toPetMmr(im, pad=True, dtype=np.float32, outres="mMR"):
  log = logging.getLogger(__name__)

  out_res = getattr(Res, outres)
  out_shape = getattr(Shape, outres)

  # PET
  # res = np.zeros(im.shape, dtype=dtype)
  res = np.zeros_like(im)
  for attr in Pet.attrs:
    log.debug("PET:%s:%d" % (attr, getattr(Pet, attr)))
    res[Act.indices(im, attr)] = getattr(Pet, attr)

  # muMap
  muMap = np.zeros(im.shape, dtype=dtype)
  muMap[im != 0] = mu_tissue_1_cm
  muMap[Act.indices(im, "bone")] = mu_bone_1_cm

  # MR
  # t1 = np.zeros(im.shape, dtype=dtype)
  t1 = np.zeros_like(im)
  for attr in T1.attrs:
    log.debug("T1:%s:%d" % (attr, getattr(T1, attr)))
    t1[Act.indices(im, attr)] = getattr(T1, attr)
  # t2 = np.zeros(im.shape, dtype=dtype)
  t2 = np.zeros_like(im)
  for attr in T2.attrs:
    log.debug("T2:%s:%d" % (attr, getattr(T2, attr)))
    t2[Act.indices(im, attr)] = getattr(T2, attr)

  # resize
  new_shape = np.rint(np.asarray(im.shape) * Res.brainweb / out_res)
  padLR, padR = divmod((np.array(out_shape) - new_shape), 2)

  def resizeToMmr(arr):
    # oldMax = arr.max()
    # arr = arr.astype(np.float32)
    # arr /= arr.max()
    arr = resize(arr, new_shape,
                 order=1, mode="constant", anti_aliasing=False)
    if pad:
      arr = np.pad(arr, [(p, p + r) for (p, r)
                         in zip(padLR.astype(int), padR.astype(int))],
                   mode="constant")
    if arr.dtype == np.uint16:
      return np.asarray(arr, dtype=np.float32) * np.float32(2 ** 16)
    return arr.astype(dtype)

  return [resizeToMmr(i) for i in [res, muMap, t1, t2]]


def matify(mat, dtype=np.float32, transpose=None):
  """@param transpose  : tuple<int>, (default: range(mat.ndim)[::-1])"""
  if transpose is None:
    transpose = tuple(range(mat.ndim)[::-1])
  return mat.transpose(transpose).astype(dtype)


def main(args):
  log = logging.getLogger(__name__)
  np.random.seed(2 ** 30)

  noiseLevels = map(float, args.noise.split(','))
  # sigma = float(args.sigma)
  resolution = getattr(Res, args.outres)
  log.debug(resolution)
  dat = [noise(im, 0 if i == 1 else n,  # no noise for muMap
               warn_zero=False, sigma=args.sigma)
         for f in tqdm(args.rawbinary, unit="file")
         for n in tqdm(noiseLevels, desc="noise", leave=False,
                       disable=len(noiseLevels) < 2)
         for i, im in enumerate(toPetMmr(getRaw(f), pad=args.save,
                                         outres=args.outres))]

  if args.save:
    from hdf5storage import write as savemat
    from caspyr._script_showslice import stripEnd
    i = 0
    for f in args.rawbinary:
      for n in noiseLevels:
        fSave = (stripEnd(f, ".raws") +
                 "_sigma%.3g_noise%.3g_%s.mat" % (args.sigma, n, args.outres))
        log.info("save:" + fSave)
        savemat({u"MultiMaps_Ref":
                 {u"PET": matify(dat[i]),
                  u"uMap": matify(dat[i + 1]),
                  u"T1": matify(dat[i + 2]),
                  u"T2": matify(dat[i + 3])}},
                '.', fSave,
                matlab_compatible=True, truncate_existing=True)
        i += 4
    return

  app = MainApp(args, dat=dat)
  args.ax_set = [dict(title="PET"), dict(title="mm"),
                 dict(title="T1"), dict(title="T2")]
  args.ax_set *= len(dat) // 4
  app.mainloop()


if __name__ == "__main__":
  args = argopt(__doc__).parse_args()
  args.parfile = "rawShow"
  # args.dpi = int(args.dpi)
  logging.basicConfig(level=getattr(logging, args.log, logging.INFO))
  log = logging.getLogger(__name__)
  # log.debug(args.__dict__['d'])
  log.debug(args)
  main(args)
