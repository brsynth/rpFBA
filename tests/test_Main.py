from unittest import TestCase
from os       import path as os_path
from rpfba    import runFBA, runFractionReaction, runParsimoniousFBA
from brs_libs import rpSBML

class TestMain(TestCase):

    """
    @classmethod
    def setUpClass(self):
    """

    def test_runFBA(self):
        rpsbml = rpSBML(os_path.join('data', 'merged.xml'))
        print(rpSBML)
        print(rpSBML.getModel)
        obj_value, status = runFBA(rpsbml, 'RP1_sink')
        self.assertAlmostEqual(obj_value, 9.230769230769237)
        self.assertTrue(status)
        # make sure that the results are written to the file
        all_json = rpsbml.genJSON()
        self.assertAlmostEqual(all_json['pathway']['brsynth']['fba_obj_RP1_sink']['value'], 9.230769230769237)

    def test_runFractionReaction(self):
        rpsbml = rpSBML(os_path.join('data', 'merged.xml'))
        obj_value, status = runFractionReaction(rpsbml, 'biomass', 1.0, 'RP1_sink', 1.0)
        self.assertAlmostEqual(obj_value, 2.3076923076923888)
        self.assertTrue(status)
        # make sure that the results are written to the file
        all_json = rpsbml.genJSON()
        self.assertAlmostEqual(all_json['pathway']['brsynth']['fba_obj_RP1_sink__restricted_biomass']['value'], 2.3076923076923888)
        self.assertAlmostEqual(all_json['pathway']['brsynth']['fba_obj_biomass']['value'], 3.6794124272706443)

    def test_runParsimoniousFBA(self):
        rpsbml = rpSBML(os_path.join('data', 'merged.xml'))
        obj_value, status = runParsimoniousFBA(rpsbml, 'RP1_sink')
        self.assertAlmostEqual(obj_value, 859.3846153846168)
        self.assertTrue(status)
        # make sure that the results are written to the file
        all_json = rpsbml.genJSON()
        self.assertAlmostEqual(all_json['pathway']['brsynth']['fba_obj_RP1_sink']['value'], 859.3846153846168)
