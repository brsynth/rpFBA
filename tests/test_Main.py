from unittest import TestCase
from os       import path as os_path
from rpfba    import rp_fba, rp_fraction, rp_pfba
from brs_libs import rpSBML

class TestMain(TestCase):

    """
    @classmethod
    def setUpClass(self):
    """

    def test_fba(self):
        rpsbml = rpSBML(os_path.join('data', 'merged.xml'))
        obj_value, rpsbml = rp_fba(rpsbml, 'RP1_sink')
        self.assertTrue(rpsbml)
        self.assertAlmostEqual(obj_value, 9.230769230769237)
        # make sure that the results are written to the file
        all_json = rpsbml.genJSON()
        self.assertAlmostEqual(all_json['pathway']['brsynth']['fba_obj_RP1_sink']['value'], 9.230769230769237)

    def test_fraction(self):
        rpsbml = rpSBML(os_path.join('data', 'merged.xml'))
        obj_value, rpsbml = rp_fraction(rpsbml, 'biomass', 1.0, 'RP1_sink', 1.0)
        self.assertTrue(rpsbml)
        self.assertAlmostEqual(obj_value, 2.3076923076923888)
        # make sure that the results are written to the file
        all_json = rpsbml.genJSON()
        self.assertAlmostEqual(all_json['pathway']['brsynth']['fba_obj_RP1_sink__restricted_biomass']['value'], 2.3076923076923888)
        self.assertAlmostEqual(all_json['pathway']['brsynth']['fba_obj_biomass']['value'], 3.6794124272706443)

    def test_pfba(self):
        rpsbml = rpSBML(os_path.join('data', 'merged.xml'))
        obj_value, rpsbml = rp_pfba(rpsbml, 'RP1_sink')
        self.assertTrue(rpsbml)
        self.assertAlmostEqual(obj_value, 859.3846153846168)
        # make sure that the results are written to the file
        all_json = rpsbml.genJSON()
        self.assertAlmostEqual(all_json['pathway']['brsynth']['fba_obj_RP1_sink']['value'], 859.3846153846168)
