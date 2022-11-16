from   glob   import glob
from   osgeo  import gdal, gdal_array
from   abc    import ABC, abstractmethod
from ..base   import SARImage, XMLMetadata
from   typing import List

import os


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

    @bands.setter
    def bands(self, bands: List):
        self._bands = bands


class SLC(S1SARImage):
    def __init__(self, SAFE: str):
        "SAFE: Path to SAFE directory"
        super().__init__(SAFE)
        self._SAFE    = SAFE
        self.__swaths = [
                SubSwath(i, self) for i in range(1, self._num_swaths+1)
                ]

    def __getitem__(self, idx):
        pass
    
    @property
    def bands(self):
        return self._bands


class SubSwath(SLC):
    def __init__(self, i: int, parent: SLC):
        self.__parent = parent
        self._iw      = i
        self._bands   = [Band(band, parent, i) for band in self.__parent._bands]
    
    @property
    def bands(self):
        return self._bands


class Band(SubSwath):
    def __init__(self, band: str, parent: SLC, _iw: int):
        base =(f'{parent.platform.lower()}-'
               f'{parent.mode.lower()}{_iw}-'
               f'{parent.product.lower()}-'
               f'{band.lower()}-*')
        self.__measurement = Measurement(
                glob(os.path.join(parent._SAFE, 'measurement', base))[0]
                )
        self.__annotation  = Annotation(
                glob(os.path.join(parent._SAFE, 'annotation', base))[0]
                ) 


class Measurement:
    def __init__(self, path: str):
        self.__file = os.path.split(path)[-1]
        print(path)


class Annotation(XMLMetadata):
    def __init__(self, path: str):
        super().__init__(path)
        self.__file = os.path.split(path)[-1]


