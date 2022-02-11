import unittest
from scripts.processOCV import processOCV


class Test(unittest.TestCase):

    def test_SOC_temp_to_OCV(self):

        model = processOCV()
        model.run()

        test_values = [
            (0.55,33,3.88),
            (0.95,5,4.07),
            (0.11,-9,3.54),
            (0.50,25,3.86),
            (1.00,50,4.15),
        ]
        for z,t,v in test_values:
            self.assertAlmostEqual(model.SOC_temp_to_OCV(z,t)[0], v, places=1, msg="Should be %s" % str(v))

if __name__ == '__main__':
    test_SOC_temp_to_OCV()
