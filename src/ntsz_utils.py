import os
import numpy as np
from pathlib import Path
from astropy.io import ascii
from astropy.io import fits


def get_data_path(subfolder: str = None):
    env_data_path = os.getenv("NTSZ_DATA_DIR")
    if env_data_path:
        base_path = Path(env_data_path).expanduser().resolve()
    else:
        base_path = Path(__file__).resolve().parent
    if subfolder:
        return base_path / subfolder
    return base_path


def get_path(name, path):
    """A function to get path to a file with a name
    located in any of the subfolders in the path specified by user.

    Parameters
    ----------
    name: str
            Name of the file
    path: str
            Bash variable assigning the path to directory containing file.
    Returns
    ---------
    Path to file with filename <name>.
    """
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)


def get_planck_freq_index(frequency):
    frequency_index = {'30': 0, '44': 1, '70': 2, '100': 3, '143': 4, '217': 5,
                       '353': 6, '545': 7, '857': 8}
    return frequency_index[str(frequency)]


def get_planck_bandpass(frequency):
    """ A function to obtain the frequency and transmission columns from the
    Planck IMO which are used to apply bandpass corrections.

    Parameter
    ---------
    frequency: int
               Nominal frequency of a Planck frequency band in GHz. Choose
               from 30, 44, 70 for LFI or 100, 143, 217, 353, 545, 857 for HFI.
    Returns
    --------
    frequencies: A list of frequencies in GHz.
    transmission: A list of corresponding transmission.
    """

    dir_path = Path(__file__).resolve().parent
    bandpass_path = get_path(name="planck_bandpass.fits", path=dir_path)
    ind = get_planck_freq_index(frequency)
    hdu = fits.open(bandpass_path)
    freq = hdu[ind+1].data['Frequency']
    trans = hdu[ind+1].data['Transmission']
    return freq, trans


# Function to write data into an output file.
def writefile(filename, data, names=None, overwrite=True):
    """A function to write data into a file.

    Parameters
    ----------
    filename: str
            Name of file for data to be written to.
    data
            Data to be written to.
    names: str, optional
            Name of columns to be written.
    overwrite: bool, optional
            Whether to overwrite an already existing file with same filename.
            Default is True.
    """
    ascii.write(
        data,
        filename,
        format="fixed_width",
        delimiter_pad=" ",
        delimiter=None,
        fast_writer=False,
        overwrite=overwrite,
        names=names,
    )
    return None


# Function to create fits file
def create_fits(file_name, data):
    """A function to A function to write data into a fits file.

    Parameters
    ----------
    file_name: str
            Name of file for data to be written to.
    data
            Data to be written to.
    """
    hdu = fits.PrimaryHDU()
    hdu.data = np.array(data, dtype=np.float32)
    hdu.writeto(file_name, overwrite=True)
    return None


# Function to read the above created fits file
def read_fits(file_name):
    """A function to read extracted cluster patches saved as fits file.

    Parameters
    ----------
    file_name: str
            Path to the fits file.
    Returns
    -------
    data: float array
            Returns the data in the fits file as a tuple or array
    """
    hdul = fits.open(file_name)
    map = hdul[0].data
    hdul.close()
    return map
