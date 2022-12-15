from   glob   import glob
from   osgeo  import gdal, gdal_array
from   abc    import ABC, abstractmethod
from ..base   import SARImage, XMLMetadata
from   typing import List
from   copy   import copy

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
        
        for i, band in enumerate(self.__parent._bands):
            self.__dict__[band] = self._bands[i]

    @property
    def bands(self):
        return self._bands

    def __getitem__(self, idx):
        "Return desired band."
        return self._bands[idx]

    def __repr__(self):
        return f"<{type(self).__name__} {self._iw} object>"


class Band:
    def __init__(self, band: str, parent: SLC, _iw: int):
        base =(f'{parent.platform.lower()}-'
               f'{parent.mode.lower()}{_iw}-'
               f'{parent.product.lower()}-'
               f'{band.lower()}-*')
        self.__name: str  = band

        self._measurement = Measurement(glob(os.path.join(parent._SAFE,
                                                          'measurement',
                                                          base))[0])
        self._annotation  = Annotation(glob(os.path.join(parent._SAFE,
                                                         'annotation',
                                                         base))[0])
        self.__bursts     = [Burst(       i, self           ) for i
                             in range(self._annotation._num_bursts)]
    
    def __getitem__(self, idx):
        "Return desired bursts of Band."
        if not isinstance(idx, slice):
            idx = slice(idx, idx + 1)
        return sum( self.__bursts[idx] )

    def __compose__(self, *args):
        "Compositor method for multiple bursts."

    def __repr__(self):
        return f"<{type(self).__name__} {self.__name} object>"


class Measurement:
    def __init__(self, path: str):
        self._path: str = path
        self._file: str = os.path.split(path)[-1]
        self._band: str = self._file.split("-")[3]

        self._ds  : gdal.Dataset = gdal.Open(path)


class Annotation(XMLMetadata):
    def __init__(self, path: str):
        super().__init__(path)
        self._file       = os.path.split(path)[-1]
        self._num_bursts = len(self.swathTiming.burstList)
        self._linespb    = int(self.swathTiming.linesPerBurst.text)        
        self._sampspb    = int(self.swathTiming.samplesPerBurst.text)
        print(self.generalAnnotation,
              self.geolocationGrid.geolocationGridPointList.geolocationGridPoint_0.line,
              self.geolocationGrid.geolocationGridPointList.geolocationGridPoint_0.pixel)

class Burst:
    "Class representing a S1 TOPSAR burst and all its attributes."
    def __init__(self, i: int, parent: Band):
        self.burst_info = parent._annotation.swathTiming.burstList[f'burst_{i}']
        self.__i        = {i}
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
        
        self._hstart: int = self.__fsample.max().item()
        self._hend  : int = self.__lsample.max().item()
        self._vstart: int = (self.__fsample > 0).argmax().item()
        self._vend  : int = (
                self.__fsample[self._vstart:] < 0
                ).argmax().item() + self._vstart

        self._vstart += self.__line
        self._vend   += self.__line

        self.__array_coords = {(self._hstart,
                                self._vstart,
                                self._hend-self._hstart,
                                self._vend-self._vstart)}
        
        self.ds = self._measurement._ds
        
        # More debugging stuff.
        print(self.__array_coords)

    @property
    def array(self):
        zeros = np.zeros((sum([x[3] for x in self.__array_coords]),
                          min([x[2] for x in self.__array_coords])),
                         dtype=np.complex64)

        _coords = np.array(sorted(self.__array_coords,
                                  key=lambda x: x[1]))
        
        i, j, p    = 0, 0, 0
        for k, coords in enumerate(_coords):
            j         += coords[3]
            coords[2]  = zeros.shape[-1]
            zeros[i:j] = gdal_array.LoadFile(self.__path, *coords.tolist())
            i         += coords[3]

        return zeros
    
    @property
    def amplitude(self):
        amp = np.abs(self.array) + 1
        amp = 10*np.log10(amp)
        return amp

    @property
    def phase(self):
        return np.angle(self.array)

    def __geolocate__(self):
        ...

    def __add__(self, other: "Burst"):
        assert self.__path == other.__path, "You can only merge bursts "
        "of the same swath and band."
        if not self == other:
            
            _copy = copy(self)
            _copy.__array_coords = _copy.__array_coords.union(
                   
                    other.__array_coords
                   
                    )
            _copy.__i = _copy.__i.union( other.__i )

        return _copy

    def __radd__(self, other):
        return self

    def __repr__(self):
        return f"<{type(self).__name__} {sorted(self.__i)} object>"



