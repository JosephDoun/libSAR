from  typing     import List
from .structures import Burst, Swath


class Deburster:
    "Class to handle debursting."
    def __init__(self, bursts: List[Burst]):
        self._bursts   = bursts
        self._pointerA = 0
        self._trackerW = 0
        self._overlaps = [self.__overlap(i, bursts[i-1], bursts[i])
                          for i in range(len(bursts))]
    
    def __len__(self):
        return len(self._bursts)
    
    def __call__(self, burst_idx: int):
        overlap = self._overlaps[burst_idx]
        self._pointerA += self._trackerW - overlap
        self._trackerW  = self._bursts[burst_idx]._Burst__src_coords[3]
        return self._pointerA, self._bursts[burst_idx].array
    
    def __iter__(self):
        self._iter     = -1
        self._pointerA =  0
        return self
    
    def __next__(self):
        self._iter += 1
        while self._iter < len(self):
            return self.__call__(self._iter)
        raise StopIteration
    
    def __overlap(self, idx, x: Burst, y: Burst):
        """
        Compare azimuth times of last and first valid samples
        of consecutive bursts. Convert the azimuth difference
        to line overlap.
        
        Returns:
            The overlap in number of lines that two consecutive
            bursts have.
        """
        if not idx: return 0
        
        overlap = (
            # Azimuth time of last valid azimuth sample of previous burst.
            x._Burst__atimes[x._Burst__src_coords[1] +
                             x._Burst__src_coords[3] -
                             x._Burst__line] -
            # Azimuth time of first valid azimuth sample of current burst.
            y._Burst__atimes[y._Burst__src_coords[1] -
                             y._Burst__line]
        # This difference is expected to te positive.
        )
        # Following lines might have to change
        # from floor division to rounding.
        overlap //= y._Burst__dt
        overlap  += 1
        return int(overlap)
    
    @property
    def overlaps(self):
        return self._overlaps
    

class SwathMerger:
    pass
