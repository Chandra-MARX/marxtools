from tempfile import TemporaryDirectory
import subprocess
import shutil
from glob import glob
from os import path

import numpy as np
from astropy.table import Table

from ..io import read_marx_file

'''
Function works (checked by hand), but test does not.
Comparing to marx2fits is a good idea in principle,
but marx2fits does other reformatting (e.g. change PHA to integer)
that makes is very general test fail.
Need to modify test to test columns that are known to work or
mimick all marx2fits effects here.
That would be a good test for consistency of marx2fits, too.


Also, fix TemporaryDirectory. Not sure what fails there.
'''

def test_read_file():
        '''Test that binary data is read correctly'''
        #with TemporaryDirectory() as d:
        d = '.'
        process = subprocess.Popen(['which', 'marx'], stdout=subprocess.PIPE)
        out, err = process.communicate()
        out = out.decode()
        shutil.copyfile(path.join(path.dirname(path.dirname(out.strip())),
                                  'share', 'marx','pfiles', 'marx.par'),
                        path.join(d, 'marx.par'))
        subprocess.call('marx')
        subprocess.call(['marx2fits', path.join(d, 'point'),
                         path.join(d, 'point.fits')])
        evt = Table.read(path.join(d, 'point.fits'), hdu=1)
        datfiles = glob(path.join(d, 'point', '*.dat'))
        for d in datfiles:
            data, name = read_marx_file(d)
            assert np.allclose(evt[name], data)
        assert len(datfiles) == len(evt.columns)
