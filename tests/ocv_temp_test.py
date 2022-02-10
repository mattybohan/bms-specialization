import unittest
from scripts.processOCV import SOC_temp_to_OCV


class Test(unittest.TestCase):

    def test_SOC_temp_to_OCV(self):
        self.assertEqual(SOC_temp_to_OCV(0.55,33), 3.88, "Should be 3.88")


if __name__ == '__main__':
    test_SOC_temp_to_OCV()
