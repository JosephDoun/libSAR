from   glob   import glob
from   osgeo  import gdal, gdal_array, osr
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
        assert SAFE.endswith('.SAFE'), "Constructor expects SAFE directory."
        SAFE = os.path.split(SAFE)[-1]
        
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
        "Return desired bursts of Band object."
        if not isinstance(idx, slice):
            idx = slice(idx, idx + 1 or None)
        returnable = self.__bursts[idx]
        return returnable[0] if len(returnable) == 1 else BurstGroup(returnable)

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
        self._dt         = float(self.imageAnnotation
                                     .imageInformation
                                     .azimuthTimeInterval
                                     .text)
        # Debugging.
        print(self.generalAnnotation,
              self.geolocationGrid.geolocationGridPointList.geolocationGridPoint_0.line,
              self.geolocationGrid.geolocationGridPointList.geolocationGridPoint_0.pixel)

class Burst:
    "Class representing a S1 TOPSAR burst and all its attributes."
    def __init__(self, i: int, parent: Band):
        self.burst_info = parent._annotation.swathTiming.burstList[f'burst_{i}']
        self._i        = {i}
        self._path     = parent._measurement._path
        self._lpb      = parent._annotation._linespb
        self._spb      = parent._annotation._sampspb
        self._id       = self.burst_info.burstId
        self._line     = self._lpb * i
        self._fsample  = np.fromstring(self.burst_info.firstValidSample.text,
                                        sep=' ',
                                        dtype=int)
        self._lsample  = np.fromstring(self.burst_info.lastValidSample.text,
                                        sep=' ',
                                        dtype=int)
        
        self._dt       = parent._annotation._dt
        self._t        = float(self.burst_info.azimuthAnxTime.text)
        self._atimes   = self._t + np.arange(self._lpb) * self._dt
        
        self.t0 = self._t
        self._to_time  = lambda l: self._t + l * self._dt
        self._to_line  = lambda t: int((t - self.t0) * 1 / self._dt + .1)
        
        self._hstart: int = self._fsample.max().item()
        self._hend  : int = self._lsample.max().item()
        self._vstart: int = (self._fsample > 0).argmax().item()
        self._vend  : int = (
                # Comment on this.
                # Find index of last valid (non negative) value.
                self._fsample[self._vstart:] < 0
                ).argmax().item() + self._vstart

        # Comment on this part.
        self._vstart += self._line
        self._vend   += self._line

        # Each burst should be explicitly bound to its super-objects.
        # Higher objects should be available for access recursively.
        # self._band = parent.__name
        
        self._src_coords = (self._hstart,
                            self._vstart,
                            self._hend-self._hstart,
                            self._vend-self._vstart,)
        
        ds = parent._measurement._ds
                
        # Assumption: Each burst has two rows of GCPs
        # associated with it, which are independent.
        GCPs = ds.GetGCPs()[:]#[i*21:(i+2)*21]
        # self.GeoTransform = gdal.GCPsToGeoTransform(self.GCPs)
        # self.InvGeoTransform = gdal.InvGeoTransform(self.GeoTransform)
        
        # self.ds = self.__build_ds(GCPs)
        
        # More debugging stuff.
        print(self._src_coords)

    def __build_ds(self, gcps: List[gdal.GCP]):
        """
        Create a Virtual Dataset that belongs to the burst.
        
        # TODO
        # REMOVE.
        
        # CHANGE THIS FUNCTION TO HANDLE GEOTRANSFORMATION.
        # i.e. Find origin.
        # TODO
        """
        # ds: gdal.Dataset = gdal.Translate('',
        #                                   srcDS=self.__path,
        #                                   options=gdal.TranslateOptions(
        #                                       srcWin=(0, self.__line,
        #                                               self.__spb, self.__lpb),
        #                                       format='VRT',
        #                                       noData=0,
        #                                   ))
        # ref = osr.SpatialReference()
        # ref.ImportFromEPSG(4326)
        # ds.SetGCPs(gcps, ref)
        
        # ds = gdal.Warp('', ds, options=gdal.WarpOptions(
        #     dstSRS='EPSG:4326',
        #     format='VRT'
        #     ))
        # gdal.BuildVRTOptions()
        return None
    
    @property
    def array(self):
        return gdal_array.LoadFile(self._path, *self._src_coords)
    
    @property
    def amplitude(self):
        amp = np.abs(self.array) + 1
        amp = 10*np.log10(amp)
        return amp

    @property
    def phase(self):
        return np.angle(self.array)

    def __repr__(self):
        return f"<{type(self).__name__} {sorted(self._i)} object>"


class BurstGroup:
    def __init__(self, bursts: List[Burst]):
        self.__bursts     = sorted(bursts, key=lambda x: x._i)
        self.__i          = [burst._i for burst in bursts]
        self.__src_coords = [burst._src_coords for burst in bursts]
        self.__overlaps   = [
            0,
            *[self.__get_overlap(
                bursts[i-1],
                bursts[i]
                ) for i in range(1, len(bursts))]
        ]
        self.__shape     = (
            # Calculate dimensions of debursted array.
            
            # Max line minus min line minus overlaps.
            self.__src_coords[-1][1] -
            self.__src_coords[0][1]  -
            sum(self.__overlaps),
            
            # Maximum width of all bursts.
            min(burst._src_coords[2] for burst in bursts) -
            max(burst._src_coords[0] for burst in bursts)
        )

    @property
    def array(self):
        array = np.zeros((self.__shape))
        for burst in self.__bursts:
            pass
        return array
    
    def __get_overlap(self, x: Burst, y: Burst) -> int:
        """
        Description: Index azimuth times of Bursts and...
        """
        
        overlap = (x._atimes[x._src_coords[1] + x._src_coords[3] - x._line] -
                   y._atimes[y._src_coords[1] - y._line])
        
        # Following lines might have to change
        # from floor division to rounding.
        overlap //= x._dt
        overlap  += 1
        
        return int(overlap)
    
    def __repr__(self):
        return f"<{type(self).__name__} {self.__i} object>"
