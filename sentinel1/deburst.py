from  typing     import List
from .structures import Burst


class Deburster:
    "Class to handle debursting."
    def __init__(self, bursts: List[Burst]):
        self._bursts   = bursts
        self._pointerA = 0
        self._trackerW = 0
    
    def __len__(self):
        return len(self._bursts)
    
    def __call__(self, burst_idx: int):
        overlap = self.__overlap(burst_idx,
                                 self._bursts[burst_idx-1],
                                 self._bursts[burst_idx])
        self._pointerA += self._trackerW - overlap
        self._trackerW  = self._bursts[burst_idx]._src_coords[3]
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
            x._atimes[x._src_coords[1] + x._src_coords[3] - x._line] -
            # Azimuth time of first valid azimuth sample of current burst.
            y._atimes[y._src_coords[1] - y._line]
        # This difference is expected to te positive.
        )
        # Following lines might have to change
        # from floor division to rounding.
        overlap //= y._dt
        overlap  += 1
        return int(overlap)
        