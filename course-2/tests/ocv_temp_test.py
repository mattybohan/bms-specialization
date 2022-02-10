import unittest
from scripts.processOCV import SOC_temp_to_OCV


class Test(unittest.TestCase):

    def test_SOC_temp_to_OCV(self):

        test_values = [
            (0.55,33,3.88),
            (0.55,33,3.88),
            (0.55,33,3.88),
        ]
        for z,t,v in test_values:
            self.assertEqual(SOC_temp_to_OCV(z,t), v, "Should be %s" % str(v))

if __name__ == '__main__':
    test_SOC_temp_to_OCV()
