from .structures import Swath
from  osgeo      import gdal


class Geolocator:
    """
    Class that handles the geolocation of a single Swath.
    """
    def __init__(self, swath: Swath) -> None:
        self.swath = swath
        self.ds    = swath.ds
        self.GCPsReference  = self.ds and self.ds.GetGCPSpatialRef()
        self.GCPsProjection = self.ds and self.ds.GetGCPProjection()
        self.GCPs           = self.ds and self.ds.GetGCPs()
        
        print('Reference:')
        print(self.GCPsReference)
        print('Projection:')
        print(self.GCPsProjection)
        print('Bursts:')
        print(swath._bursts[0]._vstart, swath._bursts[0]._vend)
        print(swath._bursts[1]._vstart, swath._bursts[1]._vend)
    
    def move_line(self, GCP: gdal.GCP):
        return GCP
    
    def move_pixel(self, GCP: gdal.GCP):
        return GCP