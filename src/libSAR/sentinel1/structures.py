from   glob        import glob
from   osgeo       import gdal, gdal_array, osr
from   abc         import ABC, abstractmethod
from ..shared      import SARImage, XMLMetadata, InvalidInputError
from   typing      import List

import os
import numpy  as     np


class SAFEDirectory:
    """
    Class for the representation and parsing of SAFE directories.
    
    Responsible for directory validity checks.
    """
    
    assertions = lambda SAFE: all([
        os.path.exists(os.path.join(SAFE, 'annotation')),
        os.path.exists(os.path.join(SAFE, 'measurement')),
        os.path.exists(os.path.join(SAFE,'manifest.safe')),
        os.path.exists(os.path.join(SAFE, 'preview')),
        # os.path.exists(os.path.join(SAFE, os.path.split(SAFE)[-1])),
    ])
    
    def __init__(self, SAFE: str) -> None:
        try:
            "Parsing."
            assert SAFEDirectory.assertions(SAFE)
            self._SAFE: str     = SAFE
            self.manifest: Manifest = Manifest(os.path.join(SAFE,
                                                            'manifest.safe'))
            self._measurements: List[str] = glob(os.path.join(SAFE,
                                                             'measurement/*.tif?'))
            self._annotations: List[str]  = glob(os.path.join(SAFE,
                                                             'annotation/*.xml'))
        except AssertionError:
            "Fail."
            raise InvalidInputError('Does not appear to be a valid SAFE directory.')


class Measurement:
    def __init__(self, path: str):
        self.measurement: str = path
        
    @property
    def ds(self) -> gdal.Dataset:
        return gdal.Open(self.measurement, gdal.GA_ReadOnly)


class Manifest(XMLMetadata):
    def __init__(self, path: str):
        super().__init__(path)
        """
        # TODO
        # Should be populated with
        # interior manifest info.
        """


class Annotation(XMLMetadata):
    def __init__(self, path: str):
        super().__init__(path)
        self._num_bursts = len  (self.swathTiming.burstList)
        self._linespb    = int  (self.swathTiming.linesPerBurst.text)        
        self._samplespb  = int  (self.swathTiming.samplesPerBurst.text)
        self._dt         = float(self.imageAnnotation
                                     .imageInformation
                                     .azimuthTimeInterval
                                     .text)
        

class S1SARImage(SARImage, SAFEDirectory):
    "Sentinel-1 SAR image parser class."
    
    __NSWATHS = {"WV": 2,
                 "IW": 3,
                 "EW": 5,
                 "SM": 6}
    
    def __init__(self, SAFE: str):
        SARImage     .__init__(self, SAFE)
        SAFEDirectory.__init__(self, SAFE)
        
        # Isolate filename and format for parsing.
        SAFE           = os.path.split(SAFE)[-1]
        metadata: list = SAFE.replace('.', '_').split('_')
        
        # Remove empty string.
        '' in metadata and metadata.remove('')
        
        (self.platform,
         self.mode,
         self.product,
         
         (  
            self.p_level,
            self.p_class,
            # Single or Double Polarisation.
            __sod,
            # Polarisation.
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
        

class SLC(S1SARImage):
    """
    Single-Look Complex SAR Image class implemente for
    the Sentinel-1 platforms.
    """
    def __init__(self, SAFE: str):
        "SAFE: Path to SAFE directory"
        super().__init__(SAFE)
        
        # Build bands.
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
    "class Band: Contains multiple Swaths."
    def __init__(self, band: str, parent: SLC):
        self.name:      str  = band
        self._slc:      SLC  = parent
        self.__swaths: Swath = [
            
            # Build Swaths.
            Swath(swathn, self) for swathn
                
            in range(1, parent._num_swaths + 1)
            
        ]
        self.__merger: SwathMerger = SwathMerger()
    
    def __getitem__(self, idx):
        "Return desired Swath."
        return self.__swaths[idx]

    def __assemble__(self):
        pass
    
    def __repr__(self):
        return f"<{type(self).__name__} {self.name} object>"


class Swath(Measurement, Annotation):
    """
    class Swath:
        Associated with a band, measurement and annotation object.
        It contains multiple Burst objects.
    """
    def __init__(self, i: int, parent: Band):
        self._band = parent.name
        self.__i    = i

        base = (f'{parent._slc.platform}-'
                f'{parent._slc.mode}{i}-'
                f'{parent._slc.product}-'
                f'{parent.name}-'
                
                # Start/stop times are different for each swath.
                # f'{parent._slc.start_time}-'
                # f'{parent._slc.stop_time}-'
                # Times can be used to derive order.
                
                '*'
                
                f'{parent._slc.abs_orbit}-'
                f'{parent._slc.mission_data_take_id}-'
                
                # ddd.ext 3 digit incremental number
                # plus file extension.
                '*'
                .lower())

        Measurement.__init__(self, glob(os.path.join(parent._slc._SAFE,
                                                          'measurement',
                                                                base))[0])
        
        Annotation .__init__(self, glob(os.path.join(parent._slc._SAFE,
                                                          'annotation',
                                                                base))[0])

        self.__bursts = [
            
            # Build Bursts.
            Burst(i, self) for i in range(self._num_bursts)
        
        ]
        
        self._geoloc = Geolocator(self)

    @property
    def band(self):
        "Returning the name of the band it belongs to."
        # Should it return the name or the band instance?
        return self._band
    
    @property
    def array(self):
        "# TODO Assemble the swath's array."

    def __getitem__(self, idx):
        "Return desired bursts of Band object."
        
        if not isinstance(idx, slice):
            idx = slice(idx, idx + 1 or None)
            
        returnable = self.__bursts[idx]
        # TODO
        # Should probably return a "BurstGroup" object regardless.
        # Likely rename the class to "Bursts".
        return returnable[0] if len(returnable) == 1 else BurstGroup(returnable)

    def __repr__(self):
        return f"<{type(self).__name__} {self.__i} object>"


from  .geolocation import Geolocator


class Burst:
    "Class representing a S1 TOPSAR burst and all its attributes."
    def __init__(self, i: int, parent: Swath):
        self.burst_info = parent.swathTiming.burstList[f'burst_{i}']
        self._swath     = parent
        self._band      = parent.band
        self._i         = i
        self._burstn    = i + 1
        self._path      = parent.measurement
        self._lpb       = parent._linespb
        self._spb       = parent._samplespb
        self._id        = self.burst_info.burstId
        self._line      = self._lpb * i
        self._fsample   = np.fromstring(self.burst_info.firstValidSample.text,
                                        sep=' ',
                                        dtype=int)
        self._lsample   = np.fromstring(self.burst_info.lastValidSample.text,
                                        sep=' ',
                                        dtype=int)
        
        self._dt        = parent._dt
        self._t         = float(self.burst_info.azimuthAnxTime.text)
        self._atimes    = self._t + np.arange(self._lpb) * self._dt
        
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
        
        # The array coordinates of the burst.
        self._arr_coords = (self._hstart,
                            self._vstart,
                            self._hend-self._hstart,
                            self._vend-self._vstart,)
    
    @property
    def array(self):
        "Fetch the burst out of the swath array."
        return gdal_array.LoadFile(self._path, *self._arr_coords)
    
    def __geotransform(self):
        "Derive from Swath geocoordinates. # TODO"
        pass
    
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
        return f"<{type(self).__name__} {self._burstn} object>"


from  .assembly import Deburster, SwathMerger


class BurstGroup(Burst):
    def __init__(self, bursts: List[Burst]):
        self._bursts   = []
        self._i        = []
        self._burstn   = []
        self._width    = 0
        self._height   = 0
        
        for burst in bursts:
            self._bursts.append(burst)
            self._i.append(burst._i)
            self._burstn.append(burst._burstn)
            
            # Sum up all the valid burst heights ignoring
            # overlaps. Subtract overlaps later.
            self._height += burst._arr_coords[-1]
            
            # Update to the minimum width.
            if burst._arr_coords[2] < self._width or not self._width:
                # Minimum width ensures no blank space.
                # in range direction.
                self._width = burst._arr_coords[2]
        
        self.__deburst   = Deburster(bursts)
        self.__shape     = (
            
            # Sum of burst heights minus overlaps.
            self._height - sum(self.__deburst.overlaps),
            
            # Minimum of all burst widths.
            # Temporary solution, to be revisited.
            self._width
        )
        
    @property
    def array(self):
        array = np.zeros((self.__shape), dtype=np.complex64)
        
        for p, barray in self.__deburst:      
            # Index destination array.
            array[# Starting from 0,
                  # Number of burst lines minus overlap
                  # with previous burst.
                  p:barray.shape[0] + p,
                  # Enforce width.
                   :self.__shape[1]] = barray[:,
                                              # Enforce width.
                                              :self.__shape[1]]
        return array
    
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
        return f"<{type(self).__name__} {self._burstn} object>"


class GRD(S1SARImage):
    pass
