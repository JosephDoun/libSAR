from   glob    import glob
from   osgeo   import gdal, gdal_array, osr
from   abc     import ABC, abstractmethod
from ..shared  import SARImage, XMLMetadata
from   typing  import List

import os
import numpy  as     np


class S1SARImage(SARImage):
    "SENTINEL-1 SAR image parser class."
    
    __NSWATHS = {"WV": 2,
                 "IW": 3,
                 "EW": 5,
                 "SM": 6}
    
    def __init__(self, SAFE: str):
        super().__init__(SAFE)
        self.__isValidSAFE(SAFE)
        SAFE = os.path.split(SAFE)[-1]
        
        metadata: list = SAFE.replace('.', '_').split('_')
        '' in metadata and metadata.remove('')
        
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
         _) = metadata
        
        self._bands      = self._POLARISATIONS[__sod][__pol]
        self._num_swaths = self.__NSWATHS[self.mode]

    @property
    def bands(self):
        return self._bands
    
    def __isValidSAFE(self, SAFE: str):
        "Validate SAFE directory path."
        assert SAFE.endswith('.SAFE'), "Constructor expects SAFE directory."
        

class SLC(S1SARImage):
    def __init__(self, SAFE: str):
        "SAFE: Path to SAFE directory"
        super().__init__(SAFE)
        self. _SAFE  = SAFE
        self.__bands = [Band(band, self) for band in self.bands]
        
        for i, band in enumerate(self.bands):
            self.__dict__[band] = self.__bands[i]

    def __getitem__(self, idx):
        "Return desired band."
        return self.__bands[idx]
    
    def __len__(self):
        return len(self.__swaths)

    @property
    def bands(self):
        return self._bands

    def __matmul__(self, other: "SLC"):
        ...


class Band:
    def __init__(self, band: str, parent: SLC):
        
        self._name: str  = band
        self.__swaths: SubSwath    = [SubSwath(self, nswath, parent) for nswath
                                      in range(1, parent._num_swaths + 1)]
        self.__merger: SwathMerger = SwathMerger()
    
    def __getitem__(self, idx):
        "Return desired SubSwath."
        return self.__swaths[idx]

    def __assemble__(self):
        pass
    
    def __repr__(self):
        return f"<{type(self).__name__} {self._name} object>"


class SubSwath:
    def __init__(self, band: Band, i: int, slc: SLC):
        self._band = band._name
        self._slc  = slc
        self._iw    = i
        
        base =(f'{slc.platform.lower()}-'
               f'{slc.mode.lower()}{i}-'
               f'{slc.product.lower()}-'
               f'{band._name.lower()}-*')
        
        self._measurement = Measurement(glob(os.path.join(slc._SAFE,
                                                          'measurement',
                                                          base))[0])
        self._annotation  = Annotation(glob(os.path.join(slc._SAFE,
                                                         'annotation',
                                                         base))[0])
        self.__bursts     = [Burst(       i, self           ) for i
                             in range(self._annotation._num_bursts)]

    @property
    def bands(self):
        return self._bands

    def __getitem__(self, idx):
        "Return desired bursts of Band object."
        if not isinstance(idx, slice):
            idx = slice(idx, idx + 1 or None)
        returnable = self.__bursts[idx]
        return returnable[0] if len(returnable) == 1 else BurstGroup(returnable)

    def __repr__(self):
        return f"<{type(self).__name__} {self._iw} object>"


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
        # print(self.generalAnnotation,
        #       self.geolocationGrid.geolocationGridPointList.geolocationGridPoint_0.line,
        #       self.geolocationGrid.geolocationGridPointList.geolocationGridPoint_0.pixel)

class Burst:
    "Class representing a S1 TOPSAR burst and all its attributes."
    def __init__(self, i: int, swath: SubSwath):
        self.burst_info = swath._annotation.swathTiming.burstList[f'burst_{i}']
        self.__swath     = swath
        self.__band      = swath._band
        self.__i         = i
        self.__path      = swath._measurement._path
        self.__lpb       = swath._annotation._linespb
        self.__spb       = swath._annotation._sampspb
        self.__id        = self.burst_info.burstId
        self.__line      = self.__lpb * i
        self.__fsample   = np.fromstring(self.burst_info.firstValidSample.text,
                                        sep=' ',
                                        dtype=int)
        self.__lsample   = np.fromstring(self.burst_info.lastValidSample.text,
                                        sep=' ',
                                        dtype=int)
        
        self.__dt        = swath._annotation._dt
        self.__t         = float(self.burst_info.azimuthAnxTime.text)
        self.__atimes    = self.__t + np.arange(self.__lpb) * self.__dt
        
        self.__hstart: int = self.__fsample.max().item()
        self.__hend  : int = self.__lsample.max().item()
        self.__vstart: int = (self.__fsample > 0).argmax().item()
        self.__vend  : int = (
                # Comment on this.
                # Find index of last valid (non negative) value.
                self.__fsample[self.__vstart:] < 0
                ).argmax().item() + self.__vstart

        # Comment on this part.
        self.__vstart += self.__line
        self.__vend   += self.__line

        # Each burst should be explicitly bound to its super-objects.
        # Higher objects should be available for access recursively.
        # self._band = parent.__name
        
        # The array coordinates of the burst.
        self.__src_coords = (self.__hstart,
                             self.__vstart,
                             self.__hend-self.__hstart,
                             self.__vend-self.__vstart,)
        
        # The dataset of the corresponding Swath - Band pair.
        ds = swath._measurement._ds
                
        # Comment on this from docs/sources.
        # Experimental.
        self.GCPs = ds.GetGCPs()[i*21:(i+1)*21] if ds else []
    
    @property
    def array(self):
        "Fetch the burst out of the swath array."
        return gdal_array.LoadFile(self.__path, *self.__src_coords)
    
    @property
    def amplitude(self):
        "Return the waves' amplitude in log scale."
        amp = np.abs(self.array) + 1
        amp = 10*np.log10(amp)
        return amp

    @property
    def phase(self):
        "Return the waves' phase."
        return np.angle(self.array)

    def __repr__(self):
        return f"<{type(self).__name__} {self.__i} object>"


from  .assembly import Deburster, SwathMerger


class BurstGroup(Burst):
    def __init__(self, bursts: List[Burst]):
        self.__bursts   = []
        self.__i        = []
        self.__width    = 0
        self.__height   = 0
        
        for burst in bursts:
            self.__bursts.append(burst)
            self.__i.append(burst._Burst__i)
            # Sum up all the valid burst heights ignoring
            # overlaps. Subtract overlaps later.
            self.__height += burst._Burst__src_coords[-1]
            
            # Update to the minimum width.
            if burst._Burst__src_coords[2] < self.__width or not self.__width:
                self.__width = burst._Burst__src_coords[2]
        
        self.__deburst   = Deburster(bursts)
        self.__shape     = (
            
            # Sum of burst heights minus overlaps.
            self.__height - sum(self.__deburst.overlaps),
            
            # Minimum of all burst widths.
            self.__width
        )
        
        self.__mod_gcps()

    @property
    def array(self):
        array = np.zeros((self.__shape), dtype=np.complex64)
        
        for p, barray in self.__deburst:      
            # Index destination array.
            array[# Starting from 0,
                  # Number of burst lines minus overlap
                  # with previous burst.
                  p:barray.shape[0] + p,
                  # Fixed width.
                   :self.__shape[1]] = barray[:,
                                              # Ensure correct width.
                                              :self.__shape[1]]
        return array
    
    def __mod_gcps(self):
        "Experimental. TO BE REMOVED."
        self.GCPs = []
        for i in range(1, len(self.__bursts)):
            gcps = self.__bursts[i].GCPs
            OL   = self.__bursts[i-1]._Burst__atimes[-1] - self.__bursts[i]._Burst__atimes[0]
            OL //= self.__bursts[i]._Burst__dt
            OL   = int(OL + 1)
            for gcp in gcps:
                gcp.GCPLine  -= OL
                gcp.GCPPixel -= self.__bursts[i-1]._Burst__src_coords[0]
                
                self.GCPs.append(gcp) if self.__bursts[i-1]._Burst__src_coords[-2] >= gcp.GCPPixel >= 0 else None
    
    def save(self, filename: str):
        "Persistance method. Quick implementation for debugging and other uses."
        "# TO BE REMOVED."
        driver: gdal.Driver  = gdal.GetDriverByName('GTiff')
        ds    : gdal.Dataset = driver.Create(filename,
                                             xsize=self.__shape[1],
                                             ysize=self.__shape[0],
                                             bands=1,
                                             eType=gdal.GDT_Float32)
        ds.WriteArray(self.amplitude)
        ref = osr.SpatialReference()
        ref.ImportFromEPSG(4326)
        ds.SetGCPs(self.GCPs, ref)
        ds.FlushCache()
    
    def __repr__(self):
        return f"<{type(self).__name__} {self.__i} object>"


class GRD(S1SARImage):
    pass
