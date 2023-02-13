from .structures import Swath, Burst
from  osgeo      import gdal
from  typing     import List, Tuple


class SimpleGeoreferencer:
    """
    Class that handles the geolocation of a single Swath.
    
        Currently makes use of the geolocation grid provided
        with the SLC product. Not accurate.
    """
    def __init__(self, swath: Swath) -> None:
        self._swath = swath
        self._ds    = swath.ds
        self._GCPsReference        = self._ds and self._ds.GetGCPSpatialRef()
        self._GCPsProjection       = self._ds and self._ds.GetGCPProjection()
        self._GCPs: List[gdal.GCP] = self._ds and self._ds.GetGCPs()
        self._dt = swath._azdt
        self._OL = [self._overlap(swath._bursts[i-1],
                                  swath._bursts[i])
                    for i in range(1, len(swath))]

        for i in range(1, len(swath)):
            # "Remove overlap from Ground Control Points."
            if not self._GCPs: break; 
            self._move_line(i, self._GCPs[
                                i*21:
                                # Throw in the final row of GCPs
                                # at last burst (21 extra for IW).
                                (i+1)*21 + ((i == len(swath) - 1) and 21)])
        del self._ds
    
    def _overlap(self, burst_a: Burst, burst_b: Burst):
        return (burst_a._atimes[-1] - burst_b._atimes[0]) // self._dt + 1
    
    def _move_line(self, i: int, GCPs: List[gdal.GCP]):
        "Get overlap and remove it from gcps."
        offset = sum(self._OL[:i])
        for GCP in GCPs:
            GCP.GCPLine -= offset
            
    def _move_origin(self, geo_transformation: Tuple[float]):
        """
        # [1]: X_geo = GT(0) + X_pixel * GT(1) + Y_line * GT(2)
        # [2]: Y_geo = GT(3) + X_pixel * GT(4) + Y_line * GT(5)
        # topleft_x, pixel_width, x_rotation, topleft_y, y_rotation, pixel_height
        # Source: https://gdal.org/tutorials/geotransforms_tut.html
        """
        x_offset = self._swath._bursts[0]._arr_coords[0]
        y_offset = self._swath._bursts[0]._arr_coords[1]
        geo_transformation = (
            geo_transformation[0] + x_offset * geo_transformation[1] + y_offset * geo_transformation[2],
            geo_transformation[1],
            geo_transformation[2],
            geo_transformation[3] + x_offset * geo_transformation[4] + y_offset * geo_transformation[5],
            geo_transformation[4],
            geo_transformation[5]
        )
        return geo_transformation
    
    def GetGeoTransform(self):
        return self._move_origin(gdal.GCPsToGeoTransform(self._GCPs))
    
    def GetProjection(self):
        return self._GCPsProjection
    
    def GetSpatialReference(self):
        return self._GCPsReference


class DEMGeoreferencer:
    "Georeferencer implementation making use of high-resolution DEM."
    # TODO
