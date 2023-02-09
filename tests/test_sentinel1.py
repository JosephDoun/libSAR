from libSAR.sentinel1.structures import *
from libSAR.sentinel1.assembly   import *

import os
import pytest

os.chdir('tests')


DATASLC = {
    
    "descending": "test_data/S1A_IW_SLC__1SDV_20230111T060136_20230111T060203_046732_059A26_F5B0.SAFE",
    "ascending" : ""
    
}

DATAGRD = {
    
    "descending": None,
    "ascending" : None
    
}

DATADUMMY = {
    "1": "test_data/S1B_EW_SLC__1SDH_20190515T888888_20190517T999999_023444_555B11_K9H2.SAFE"
}


class TestSAFEDirectory:
    "Run & Test SAFEDirectory class."
    def test_safe(self):
        SAFEDirectory(DATASLC['descending'])


class TestS1SARImage:
    "Test class for S1SARImage class."    
    def test_name_parsing_slc(self):
        S1SARSLC: S1SARImage = S1SARImage(DATASLC['descending'])
        assert S1SARSLC.platform == 'S1A'
        assert S1SARSLC.mode == 'IW'
        assert S1SARSLC.product == 'SLC'
        assert S1SARSLC.p_level == '1'
        assert S1SARSLC.p_class == 'S'
        assert S1SARSLC.start_time == '20230111T060136'
        assert S1SARSLC.stop_time == '20230111T060203'
        assert S1SARSLC.abs_orbit == '046732'
        assert S1SARSLC.mission_data_take_id == '059A26'
        assert S1SARSLC.unique_id == 'F5B0'
        assert S1SARSLC.bands == ['VV', 'VH']
        
    def test_name_parsing_grd(self):
        pass
    
    def test_exists(self):
        "Test for non existent directory."
        pytest.raises(FileNotFoundError, S1SARImage, DATADUMMY['1'])
    

class TestSLC:
    "Test class for SLC class."
    
    desc_safe = "test_data/S1A_IW_SLC__1SDV_20230111T060136_20230111T060203_046732_059A26_F5B0.SAFE"
    asce_safe = "test_data/..."
    
    slc_desc: SLC = SLC(desc_safe)
    # slc_asc:  SLC = SLC(asce_safe)
    
    def test_s1sar_name_parsing(self):
        """
        Test for name parsing.
        """
        assert isinstance(self.slc_desc, SLC)


class TestDeburster:
    "Test class for Deburster class."
    
    def test_1(self):
        pass

