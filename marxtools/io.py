from glob import glob
from os import path
from struct import unpack
from warnings import warn

import numpy as np


__all__ = ['read_marx_file',
           'marx_bin_to_table',
           ]

marxtypes = {b'A': np.int8,  # byte order not relevant for 1 byte
             b'I': np.dtype('>i2'),
             b'J': np.dtype('>i4'),
             b'E': np.dtype('>f4'),
             b'D': np.dtype('>f8')}


def read_marx_file(filename):
    '''Read binary marx file into numpy array

    Parameters
    ----------
    filename : string
        Filename and path to marz binary file

    Returns
    -------
    data : np.ndarray
        Data read from file
    colname : string
        Name of data column
    '''
    with open(filename, mode='rb') as file:
        fileContent = file.read()
    m1, m2, m3, m4, marxtype, colname, nrow, ncol, reserved = unpack('>BBBBc15siii', fileContent[:32])
    if not (m1, m2, m3, m4) == (131, 19, 137, 141):
        raise IOError('{} is not a valid marx file.'.format(filename))
    data = np.frombuffer(fileContent, dtype=marxtypes[marxtype], offset=32)
    if ncol > 0:
        data = data.reshape(nrow, ncol)
    return data, colname.rstrip(b'\x00').decode()


parconverter = {
    'b' : lambda x : 'yes' == x.strip(),
    'f' : lambda x : x,
    'i' : int,
    'r' : float,
    's' : lambda x : x,
    }


def marx_bin_to_table(dirname, cols=None):
    '''Read binary marx simulation output into `astropy.table.Table`.

    This routine is mainly useful for marx development. When using marx
    in simulating Chandra data, it is recommended to use
    `marx2fits` to convert marx output to fits files and process those.

    Missing datafiles are skipped with a warning message. It also read the
    ``marx.par`` file and places the as-run parameters in the meta
    information of the returned table

    TODO: (currently all of them as string,
    but that could be made more clever easily by using the format information in
    the .par file - just need )

    Parameters
    ----------
    dirname : string
        Directory name and path to marx native simulation output
    cols : None or list of strings
        If ``None``, all of marx output vectors will be read, otherwise
        this can be a list of strings, e.g. ``['xpos', 'energy']`` will
        read just the files ``xpos.dat`` and ``energy.dat``.

    Returns
    -------
    tab : `astropy.table.Table`
        Table with photon data
    '''
    from astropy.table import Table
    tab = Table()
    if cols is None:
        cols = glob(path.join(dirname, '*.dat'))
    for dat in cols:
        try:
            data, colname = read_marx_file(path.join(dirname, dat + '.dat'))
        except IOError as e:
            warn('Skipping column ' + dat + 'because: ' + e.message)
        tab[colname] = data

    marxpartab = Table.read(path.join(dirname, 'marx.par'),
                            format='ascii.no_header')
    for r in marxpartab:
        tab.meta[r['col1']] = parconverter[r['col2']](r['col4'])
    return tab
