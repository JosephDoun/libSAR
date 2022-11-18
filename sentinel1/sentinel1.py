from   glob   import glob
from   osgeo  import gdal, gdal_array
from   abc    import ABC, abstractmethod
from ..base   import SARImage, XMLMetadata
from   typing import List

import os
import numpy  as     np


class S1SARImage(SARImage):
    "SENTINEL-1 SAR image parser class."
    
    __NSWATHS = {"WV": 2,
                 "IW": 3,
                 "EW": 5,
                 "SM": 6}
    
    def __init__(self, SAFE: str):
        data: list = SAFE.replace('.', '_').split('_')
        '' in data and data.remove('')
        
        (self.platform,
         self.mode,
         self.product,
         
         (  
            self.p_level,
            self.p_class,
            __sod,
            __pol
         ),
         
         self.start_time,
         self.stop_time,
         self.abs_orbit,
         self.mission_data_take_id,
         self.unique_id,
         _) = data
        
        self._bands      = self._POLARISATIONS[__sod][__pol]
        self._num_swaths = self.__NSWATHS[self.mode]

    @property
    def bands(self):
        return self._bands


class SLC(S1SARImage):
    def __init__(self, SAFE: str):
        "SAFE: Path to SAFE directory"
        super().__init__(SAFE)
        self. _SAFE   = SAFE
        self.__swaths = [
                SubSwath(i, self) for i in range(1, self._num_swaths+1)
                ]

    def __getitem__(self, idx):
        "Return desired SubSwath."
        return self.__swaths[idx]
    
    def __len__(self):
        return len(self.__swaths)

    @property
    def bands(self):
        return self._bands

    def __matmul__(self, other: "SLC"):
        ...


class SubSwath:
    def __init__(self, i: int, parent: SLC):
        self.__parent = parent
        self._iw      = i
        self._bands   = [Band(band, parent, i) for band in self.__parent._bands]
    
    @property
    def bands(self):
        return self._bands

    def __getitem__(self, idx):
        "Return desired band."
        return self._bands[idx]


class Band:
    def __init__(self, band: str, parent: SLC, _iw: int):
        base =(f'{parent.platform.lower()}-'
               f'{parent.mode.lower()}{_iw}-'
               f'{parent.product.lower()}-'
               f'{band.lower()}-*')
        self._measurement = Measurement(
                glob(os.path.join(parent._SAFE, 'measurement', base))[0]
                )
        self._annotation  = Annotation(
                glob(os.path.join(parent._SAFE, 'annotation', base))[0]
                )
        self.__bursts     = [Burst(i, self) for i
                             in range(self._annotation._num_bursts)]
    
    def __getitem__(self, idx):
        "Return desired burst of Band."
        return self.__bursts[idx]


class Measurement:
    def __init__(self, path: str):
        self._path: str = path
        self._file: str = os.path.split(path)[-1]
        self._band: str = self._file.split("-")[3]


class Annotation(XMLMetadata):
    def __init__(self, path: str):
        super().__init__(path)
        self._file       = os.path.split(path)[-1]
        self._num_bursts = len(self.swathTiming.burstList)
        self._linespb    = int(self.swathTiming.linesPerBurst.text)        
        self._sampspb    = int(self.swathTiming.samplesPerBurst.text)


class Burst:
    def __init__(self, i: int, parent: Band):
        self.burst_info = parent._annotation.swathTiming.burstList[f'burst_{i}']
        self.__path     = parent._measurement._path
        self.__lpb      = parent._annotation._linespb
        self.__spb      = parent._annotation._sampspb
        self.__id       = self.burst_info.burstId
        self.__line     = self.__lpb * i
        self.__fsample  = np.fromstring(self.burst_info.firstValidSample.text,
                                        sep=' ',
                                        dtype=int)
        self.__lsample  = np.fromstring(self.burst_info.lastValidSample.text,
                                        sep=' ',
                                        dtype=int)
        self._vstart: int = (self.__fsample > 0).argmax().item()
        self._vend  : int = (self.__fsample[self._vstart:] < 0).argmax().item() + self._vstart
        self._hstart: int = self.__fsample.max().item()
        self._hend  : int = self.__lsample.max().item()

        self._vstart += self.__line
        self._vend   += self.__line

    @property
    def array(self):
        return gdal_array.LoadFile(self.__path,
                                   xoff=self._hstart,
                                   yoff=self._vstart,
                                   xsize=self._hend-self._hstart,
                                   ysize=self._vend-self._vstart)

    @property
    def amplitude(self):
        amp = np.abs(self.array) + 1
        amp = 10*np.log10(amp)
        return amp

    @property
    def phase(self):
        return np.angle(self.array)

