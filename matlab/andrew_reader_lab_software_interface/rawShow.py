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
    brainweb.raws/subject_??.bin
"""
from __future__ import division
import numpy as np
from caspyr.scripts.showslice import MainApp
from argopt import argopt
from tqdm import tqdm
import logging

from brainweb.utils import Res, getRaw, noise, toPetMmr, matify

__author__ = "Casper da Costa-Luis"
__date__ = "2018-04"
__licence__ = "MPLv2.0"
__requires__ = ["caspyr>=1.0.0"]


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

  args.out = None
  args.medfilt = 0
  args.d3 = False
  args.d4 = False
  args.anim = False
  args.n = -1
  app = MainApp(args, dat=dat)
  args.ax_set = [dict(title="PET"), dict(title="uMap", cmap="bone"),
                 dict(title="T1", cmap="Greys_r"), dict(title="T2", cmap="Greys_r")]
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
